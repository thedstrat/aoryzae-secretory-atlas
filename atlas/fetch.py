"""Data fetching: cached HTTP client plus UniProt and KEGG calls.

One module. The caching class is only ever used by the fetchers below, so
splitting them across two files bought nothing.

Licences:
  UniProt  CC-BY 4.0   - derived data redistributable with attribution
  NCBI     public domain
  KEGG     RESTRICTED  - do NOT commit bulk dumps or redistribute pathway maps
                         (data/raw/ is gitignored for this reason)
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

log = logging.getLogger(__name__)

USER_AGENT = "aoryzae-atlas/0.1 (+https://github.com/YOURUSER/aoryzae-atlas)"


# ------------------------------------------------------------------ HTTP cache

class CachedClient:
    """GET with disk caching and polite rate limiting.

    Every response is cached keyed on a hash of the URL, with a provenance
    sidecar recording the fetch timestamp. Two reasons this matters more than
    it looks: you will iterate on parsing dozens of times and should not
    re-hit UniProt for 12,000 genes each round; and when someone asks where a
    row came from, you still have the exact bytes and the date.
    """

    def __init__(self, cache_dir, rate_limit_seconds: float = 0.0, timeout: int = 60):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.rate_limit_seconds = rate_limit_seconds
        self.timeout = timeout
        self._last = 0.0
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": USER_AGENT})

    def _path(self, url: str, params: Optional[dict]) -> Path:
        key = url + json.dumps(params or {}, sort_keys=True)
        return self.cache_dir / f"{hashlib.sha256(key.encode()).hexdigest()[:20]}.txt"

    def get(self, url: str, params: Optional[dict] = None, use_cache: bool = True) -> str:
        cache_path = self._path(url, params)
        if use_cache and cache_path.exists():
            return cache_path.read_text(encoding="utf-8")

        if self.rate_limit_seconds > 0:
            elapsed = time.monotonic() - self._last
            if elapsed < self.rate_limit_seconds:
                time.sleep(self.rate_limit_seconds - elapsed)

        log.info("GET %s", url)
        resp = self._session.get(url, params=params, timeout=self.timeout)
        self._last = time.monotonic()
        resp.raise_for_status()

        cache_path.write_text(resp.text, encoding="utf-8")
        cache_path.with_suffix(".meta.json").write_text(
            json.dumps(
                {
                    "url": resp.url,
                    "status": resp.status_code,
                    "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                },
                indent=2,
            )
        )
        return resp.text


# --------------------------------------------------------------------- UniProt

def uniprot_proteome(
    client: CachedClient,
    base_url: str,
    fields: list[str],
    taxon_id: int,
    proteome: Optional[str] = None,
    page_size: int = 500,
) -> pd.DataFrame:
    """Pull the organism proteome.

    Prefers a pinned proteome accession. A taxonomy fallback is noisier -
    expect several accessions per gene. Check for fan-out before joining.
    """
    query = f"proteome:{proteome}" if proteome else f"taxonomy_id:{taxon_id}"
    if not proteome:
        log.warning("No proteome pinned; querying by taxonomy (noisier).")

    # The /search endpoint is paginated and a single request silently returns
    # only the first page. The TSV /stream endpoint returns the complete query,
    # which is essential before using coverage counts or joining annotations.
    stream_url = base_url.removesuffix("/search") + "/stream"
    text = client.get(
        stream_url,
        params={
            "query": query,
            "fields": ",".join(fields),
            "format": "tsv",
        },
    )
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
    log.info("UniProt: %d entries", len(df))
    return df


def normalise_uniprot(df: pd.DataFrame) -> pd.DataFrame:
    """Map UniProt's human-readable TSV headers onto schema names.

    These headers have changed before. If the warning fires, inspect the real
    columns and update this mapping rather than assuming it is current.
    """
    rename = {
        "Entry": "uniprot_accession",
        "Gene Names (primary)": "gene_name",
        "Gene Names": "gene_synonyms_raw",
        "Protein names": "function",
        "KEGG": "kegg_gene_id_raw",
        "GeneID": "ncbi_gene_id",
        "Gene Ontology IDs": "go_terms_raw",
        "Subcellular location [CC]": "compartment_raw",
        "Protein existence": "protein_existence",
        "Reviewed": "reviewed",
        "Sequence": "sequence",
    }
    missing = set(rename) - set(df.columns)
    if missing:
        log.warning("UniProt headers not found: %s | actual: %s",
                    sorted(missing), list(df.columns))
    return df.rename(columns={k: v for k, v in rename.items() if k in df.columns})


# ------------------------------------------------------------------------ KEGG

def verify_kegg_org(client: CachedClient, base_url: str, org: str) -> bool:
    """Confirm the org code exists before pulling anything. Cheap guard."""
    text = client.get(f"{base_url}/list/organism")
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1] == org:
            log.info("KEGG org %r verified", org)
            return True
    log.error("KEGG org %r not found - check https://rest.kegg.jp/list/organism", org)
    return False


def kegg_gene_list(client: CachedClient, base_url: str, org: str) -> pd.DataFrame:
    text = client.get(f"{base_url}/list/{org}")
    rows = [
        {"kegg_gene_id": p[0], "kegg_description": p[-1]}
        for p in (l.split("\t") for l in text.strip().splitlines())
        if len(p) >= 2
    ]
    log.info("KEGG: %d genes for %s", len(rows), org)
    return pd.DataFrame(rows)


def kegg_conv(client: CachedClient, base_url: str, org: str, target: str) -> pd.DataFrame:
    """KEGG gene ID -> outside database ID.

    target: "uniprot" | "ncbi-geneid" | "ncbi-proteinid"

    Note KEGG only covers genes assigned to orthology groups, so coverage is
    partial. Do not anchor the crosswalk on this.
    """
    text = client.get(f"{base_url}/conv/{target}/{org}")
    col = {"uniprot": "uniprot_accession", "ncbi-geneid": "ncbi_gene_id"}.get(
        target, target.replace("-", "_")
    )
    prefix = {"uniprot": "up:"}.get(target, f"{target}:")
    rows = []
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            rows.append({"kegg_gene_id": parts[0], col: parts[1].replace(prefix, "")})
    log.info("KEGG->%s: %d mappings", target, len(rows))
    return pd.DataFrame(rows)


def kegg_ko(client: CachedClient, base_url: str, org: str) -> pd.DataFrame:
    """Gene -> KEGG Orthology group.

    This is how you find gene families in A. oryzae. Symbol-prefix matching
    ("och", "mnn", "pmt") barely works because most genes have no symbol at
    all. KO numbers are cross-species orthology groups, so they find the
    OCH1 equivalent regardless of local naming.
    """
    text = client.get(f"{base_url}/link/ko/{org}")
    rows = []
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            rows.append({"kegg_gene_id": parts[0], "kegg_ko": parts[1]})
    log.info("KEGG KO links: %d", len(rows))
    return pd.DataFrame(rows)


def kegg_pathway_members(
    client: CachedClient, base_url: str, org: str,
    pathway_ids: list[str], pathway_group: str,
) -> pd.DataFrame:
    """Genes per pathway map. Long format, one row per gene-pathway pair."""
    rows = []
    for pid in pathway_ids:
        try:
            text = client.get(f"{base_url}/link/{org}/path:{org}{pid}")
        except Exception as exc:  # noqa: BLE001
            log.warning("pathway %s%s failed: %s", org, pid, exc)
            continue
        n = 0
        for line in text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                rows.append({
                    "kegg_gene_id": parts[1],
                    "kegg_pathway": f"{org}{pid}",
                    "pathway_group": pathway_group,
                })
                n += 1
        log.info("pathway %s%s (%s): %d genes", org, pid, pathway_group, n)
    return pd.DataFrame(rows)
