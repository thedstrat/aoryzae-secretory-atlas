"""Independent validation of Liu 2014 AO090 identifier mappings.

The audit deliberately does not use function or gene-name similarity to
establish identity. It compares locus identifiers and database cross-references
from a current NCBI RefSeq GFF, KEGG, and UniProt organism-wide exports.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

import pandas as pd


AO_RE = re.compile(r"^AO090\d{9}$")
EXACT_AO_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9])AO090\d{9}(?![-A-Za-z0-9])")
LEGACY_AO_TOKEN_RE = re.compile(r"AO090\d{9}-[A-Za-z0-9]+")

AUDIT_COLUMNS = [
    "original_id", "proposed_current_id", "verified_current_id",
    "mapping_method", "supporting_sources", "source_agreement",
    "alternative_candidates", "audit_status", "audit_notes",
    "manual_review_required",
]


def _split_semicolon(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [part.strip() for part in str(value).split(";") if part.strip()]


def parse_kegg(gene_list: Path, uniprot_conv: Path, ncbi_conv: Path) -> dict:
    genes: set[str] = set()
    for line in gene_list.read_text(encoding="utf-8").splitlines():
        key = line.split("\t", 1)[0]
        if key.startswith("aor:"):
            genes.add(key.removeprefix("aor:"))

    uniprot: dict[str, set[str]] = defaultdict(set)
    for line in uniprot_conv.read_text(encoding="utf-8").splitlines():
        left, right = line.split("\t")
        uniprot[left.removeprefix("aor:")].add(right.removeprefix("up:"))

    ncbi: dict[str, set[str]] = defaultdict(set)
    for line in ncbi_conv.read_text(encoding="utf-8").splitlines():
        left, right = line.split("\t")
        ncbi[left.removeprefix("aor:")].add(right.removeprefix("ncbi-geneid:"))

    return {"genes": genes, "uniprot": uniprot, "ncbi": ncbi}


def parse_uniprot(path: Path) -> dict:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    by_locus: dict[str, list[dict]] = defaultdict(list)
    legacy_by_base: dict[str, list[dict]] = defaultdict(list)

    for row in frame.to_dict("records"):
        exact = set(EXACT_AO_TOKEN_RE.findall(row.get("Gene Names", "")))
        exact.update(
            value.removeprefix("aor:")
            for value in _split_semicolon(row.get("KEGG"))
            if AO_RE.fullmatch(value.removeprefix("aor:"))
        )
        record = {
            "accession": row.get("Entry", ""),
            "loci": exact,
            "geneids": set(_split_semicolon(row.get("GeneID"))),
        }
        for locus in exact:
            by_locus[locus].append(record)
        for legacy in LEGACY_AO_TOKEN_RE.findall(row.get("Gene Names", "")):
            legacy_by_base[legacy.split("-", 1)[0]].append({**record, "legacy_name": legacy})

    return {"by_locus": by_locus, "legacy_by_base": legacy_by_base}


def _parse_gff_attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def parse_ncbi_gff(path: Path) -> dict[str, dict[str, set[str]]]:
    loci: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"geneids": set(), "transcripts": set(), "proteins": set()}
    )
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            attrs = _parse_gff_attributes(parts[8])
            locus = attrs.get("locus_tag", "")
            if not AO_RE.fullmatch(locus):
                continue
            for xref in attrs.get("Dbxref", "").split(","):
                if xref.startswith("GeneID:"):
                    loci[locus]["geneids"].add(xref.removeprefix("GeneID:"))
            if attrs.get("transcript_id"):
                loci[locus]["transcripts"].add(attrs["transcript_id"])
            if attrs.get("protein_id"):
                loci[locus]["proteins"].add(attrs["protein_id"])
    return dict(loci)


def classify_mapping(
    original_id: str,
    proposed_id: str | None,
    candidates_by_source: dict[str, set[str]],
    reverse_proposed_count: int = 1,
    relationship_hint: str | None = None,
) -> dict[str, object]:
    """Classify an identifier without using functional/name similarity."""
    available = {name: values for name, values in candidates_by_source.items() if values}
    union = set().union(*available.values()) if available else set()
    singleton_values = {next(iter(values)) for values in available.values() if len(values) == 1}
    internally_multiple = any(len(values) > 1 for values in available.values())

    if internally_multiple or len(union) > 1:
        if not internally_multiple and len(singleton_values) > 1:
            status, agreement = "conflicting", "sources disagree"
        else:
            status, agreement = "ambiguous", "multiple current candidates"
        verified = ""
        method = "split" if relationship_hint == "split" else "ambiguous"
        manual = True
    elif len(union) == 1:
        verified = next(iter(union))
        agreement = "all available sources agree" if len(available) > 1 else "single-source support"
        status = "confirmed" if len(available) > 1 else "probable"
        method = "unchanged/exact identifier" if verified == original_id else "authoritative database cross-reference"
        if reverse_proposed_count > 1 and verified != original_id:
            method = "merged"
        manual = status != "confirmed" or bool(proposed_id and proposed_id != verified)
    else:
        verified, method = "", "unresolved"
        status, agreement, manual = "unresolved", "no current-source match", True

    return {
        "verified_current_id": verified,
        "mapping_method": method,
        "source_agreement": agreement,
        "audit_status": status,
        "manual_review_required": manual,
        "candidates": union,
    }


def build_audit(
    machinery: pd.DataFrame,
    kegg: dict,
    uniprot: dict,
    ncbi: dict[str, dict[str, set[str]]],
) -> pd.DataFrame:
    proposed_counts = machinery["ao_locus_tag"].dropna().value_counts().to_dict()
    rows = []

    for source_row in machinery.to_dict("records"):
        original = str(source_row["liu_ao_locus_tag"])
        proposed = source_row.get("ao_locus_tag")
        proposed = "" if proposed is None or pd.isna(proposed) else str(proposed)

        ncbi_candidates = {original} if original in ncbi else set()
        kegg_candidates = {original} if original in kegg["genes"] else set()
        uni_records = uniprot["by_locus"].get(original, [])
        uniprot_candidates = set().union(*(r["loci"] for r in uni_records)) if uni_records else set()
        candidates = {"NCBI RefSeq": ncbi_candidates, "KEGG": kegg_candidates, "UniProt": uniprot_candidates}
        result = classify_mapping(original, proposed, candidates, proposed_counts.get(proposed, 1))

        sources = [name for name, values in candidates.items() if values]
        notes = []
        if original in ncbi:
            rec = ncbi[original]
            if rec["geneids"]:
                notes.append("NCBI GeneID:" + ",".join(sorted(rec["geneids"])))
            if rec["transcripts"]:
                notes.append("RefSeq transcript:" + ",".join(sorted(rec["transcripts"])))
            if rec["proteins"]:
                notes.append("RefSeq protein:" + ",".join(sorted(rec["proteins"])))
        if original in kegg["genes"]:
            notes.append("KEGG:aor:" + original)
            if kegg["ncbi"].get(original):
                notes.append("KEGG->GeneID:" + ",".join(sorted(kegg["ncbi"][original])))
            if kegg["uniprot"].get(original):
                notes.append("KEGG->UniProt:" + ",".join(sorted(kegg["uniprot"][original])))
        if uni_records:
            notes.append("UniProt:" + ",".join(sorted({r["accession"] for r in uni_records})))

        # A shared locus string is not sufficient if the databases attach
        # incompatible GeneIDs or UniProt accessions to that locus.
        ncbi_geneids = ncbi.get(original, {}).get("geneids", set())
        kegg_geneids = kegg["ncbi"].get(original, set())
        uniprot_geneids = set().union(*(r["geneids"] for r in uni_records)) if uni_records else set()
        geneid_sets = [values for values in (ncbi_geneids, kegg_geneids, uniprot_geneids) if values]
        kegg_uniprot = kegg["uniprot"].get(original, set())
        direct_uniprot = {r["accession"] for r in uni_records}
        id_disagreement = (
            any(values != geneid_sets[0] for values in geneid_sets[1:])
            or bool(kegg_uniprot and direct_uniprot and kegg_uniprot != direct_uniprot)
        )
        if id_disagreement:
            result["audit_status"] = "conflicting"
            result["source_agreement"] = "database identifiers disagree"
            result["manual_review_required"] = True
            notes.append("CONFLICT: current sources attach different database identifiers to the locus")

        alternatives = sorted(result.pop("candidates") - ({result["verified_current_id"]} if result["verified_current_id"] else set()))
        for record in uniprot["legacy_by_base"].get(original, []):
            alternatives.append(f"UniProt:{record['accession']} ({record['legacy_name']}; legacy suffixed name)")

        if proposed and result["verified_current_id"] and proposed != result["verified_current_id"]:
            notes.append(f"PROPOSED CORRECTION: {proposed} -> {result['verified_current_id']}")
        elif source_row.get("mapping_status") != result["audit_status"] and original == "AO090011000795":
            notes.append("PROPOSED CORRECTION: existing ambiguous classification -> confirmed unchanged locus")
        if not sources:
            notes.append("Identifier absent from fresh NCBI RefSeq, KEGG aor, and UniProt organism-wide records; no sequence-based replacement established")

        rows.append({
            "original_id": original,
            "proposed_current_id": proposed,
            **result,
            "supporting_sources": "; ".join(sources),
            "alternative_candidates": "; ".join(dict.fromkeys(alternatives)),
            "audit_notes": "; ".join(notes),
        })

    return pd.DataFrame(rows).reindex(columns=AUDIT_COLUMNS)


def run_audit(project_root: Path) -> pd.DataFrame:
    raw = project_root / "data/raw/crosswalk_audit"
    machinery = pd.read_csv(project_root / "data/processed/secretion_machinery_genes.csv", dtype=str)
    kegg = parse_kegg(raw / "kegg_aor_genes.tsv", raw / "kegg_uniprot.tsv", raw / "kegg_ncbi.tsv")
    uniprot = parse_uniprot(raw / "uniprot_aoryzae.tsv")
    ncbi = parse_ncbi_gff(raw / "GCF_000184455.2_ASM18445v3_genomic.gff.gz")
    return build_audit(machinery, kegg, uniprot, ncbi)
