"""Identifier resolution, with the mapping ladder in the right order.

THE LADDER (each step only runs on what the previous one left unresolved):

    1. exact          direct AO090... locus-tag match
    2. cross_reference via UniProt / KEGG / NCBI cross-references
    3. sequence       protein sequence match for changed gene models
    4. orthology      fresh re-derivation from the yeast ortholog
    5. unresolved     recorded, never dropped

Step 1 is the primary path, not a hopeful first try. NCBI Gene records for
RIB40 still use AO090... locus tags on the current RefSeq assembly, and KEGG
keys on aor:AO090.... AspGD going unmaintained did not kill the identifiers.

Anything resolved at step 4 is tagged record_origin=new_orthology_inference.
It is OUR 2026 inference and must stay distinguishable from Liu's canonical
assignments - never silently substituted for them.

THE INVARIANT: every input row produces at least one output row. A gene that
fails to resolve gets mapping_status=unresolved, not deletion. The count of
unresolved genes is a finding worth reporting.
"""

from __future__ import annotations

import logging
import re
import uuid

import pandas as pd

from .schema import GlycosylationRole, MappingStatus, RecordOrigin, SUBSYSTEM_ORDER

log = logging.getLogger(__name__)


_SUBSYSTEM_TEXT_RULES = {
    "dolichol_pathway": (r"\bdolichol\b",),
    "er_glycosylation": (r"oligosaccharyltransferase|n-glycan transfer|glycosylation site occupancy",),
    "gpi_biosynthesis": (r"gpi[- ]anchor biosynth|gpi mannosyltransferase",),
    "copii": (r"\bcopii\b|sec23|sec24|sec13.sec31",),
    "copi": (r"\bcopi\b|coatomer",),
    "tc": (r"translocon|signal recognition particle|signal peptidase|protein translocation",),
    "erad": (r"\berad\b|er-associated degradation",),
    "folding": (r"chaperone|protein folding|peptidyl.prolyl|disulfide isomerase|calnexin|calreticulin",),
    "golgi_processing": (r"golgi.*(?:mannosyl|glycosyl)|(?:mannosyl|glycosyl).*golgi",),
    "snare": (r"\bsnare\b",),
    "cpy_pathway": (r"\bcpy pathway\b|vacuolar protein sorting",),
    "alp_pathway": (r"\balp pathway\b",),
    "ldsv": (r"\bldsv\b",),
    "hdsv": (r"\bhdsv\b",),
}


def _text_candidates(text: object) -> list[str]:
    value = "" if pd.isna(text) else str(text)
    return sorted({name for name, patterns in _SUBSYSTEM_TEXT_RULES.items()
                   if any(re.search(pattern, value, re.I) for pattern in patterns)})


def assign_subsystems_from_evidence(
    df: pd.DataFrame, yeast_scaffold: pd.DataFrame,
) -> pd.DataFrame:
    """Fill blank subsystems from Liu descriptions, then the yeast scaffold,
    then current annotation text. Conflicting evidence is retained for review.

    ``yeast_scaffold`` must contain ``yeast_ortholog`` and ``subsystem``.
    Existing source labels are never overwritten.
    """
    out = df.copy()
    if "subsystem_source" not in out:
        out["subsystem_source"] = "unassigned"
    if "manual_review_required" not in out:
        out["manual_review_required"] = False
    out["manual_review_required"] = out["manual_review_required"].fillna(False).astype(bool)
    out["subsystem_candidates"] = ""
    out["subsystem_evidence"] = ""
    existing = out["subsystem"].notna() & out["subsystem"].astype(str).str.strip().ne("")
    out.loc[existing & out["subsystem_source"].isin(["", "unassigned"]), "subsystem_source"] = "liu2014"
    out.loc[existing, "subsystem_evidence"] = "Liu 2014 Table S1 subsystem"

    scaffold = yeast_scaffold.dropna(subset=["yeast_ortholog", "subsystem"]).copy()
    scaffold["_key"] = scaffold["yeast_ortholog"].astype(str).str.upper().str.replace(r"P$", "", regex=True)
    lookup = scaffold.groupby("_key")["subsystem"].agg(lambda x: sorted(set(x))).to_dict()

    blank = out["subsystem"].isna() | out["subsystem"].astype(str).str.strip().eq("")
    for idx in out.index[blank]:
        description = _text_candidates(out.at[idx, "Description"] if "Description" in out else "")
        ortholog_value = out.at[idx, "yeast_ortholog"]
        orthologs = re.split(r"[/,; ]+", "" if pd.isna(ortholog_value) else str(ortholog_value))
        yeast = sorted({candidate for token in orthologs for candidate in lookup.get(token.upper().removesuffix("P"), [])})
        current_text = " ".join(str(out.at[idx, c]) for c in ("annotation_function", "ncbi_annotation")
                                if c in out and pd.notna(out.at[idx, c]))
        current = _text_candidates(current_text)
        stages = (("liu2014_description", description), ("yeast_scaffold", yeast),
                  ("current_annotation", current))
        for source, choices in stages:
            if len(choices) == 1:
                out.at[idx, "subsystem"] = choices[0]
                out.at[idx, "subsystem_source"] = source
                out.at[idx, "subsystem_evidence"] = source.replace("_", " ")
                break
            if len(choices) > 1:
                out.at[idx, "subsystem_candidates"] = "; ".join(choices)
                out.at[idx, "subsystem_evidence"] = f"conflicting {source.replace('_', ' ')} evidence"
                out.at[idx, "manual_review_required"] = True
                break
    out["pathway_order"] = out["subsystem"].map(SUBSYSTEM_ORDER)
    return out


def assign_controlled_compartments(df: pd.DataFrame) -> pd.DataFrame:
    """Assign the eight user-facing compartments from direct location text,
    falling back to a conservative subsystem-level placement."""
    out = df.copy()
    out["compartment"] = "unknown"
    out["compartment_source"] = "unassigned"
    out["compartment_confidence"] = "low"
    patterns = [
        ("ER", r"endoplasmic reticulum|\ber membrane\b"), ("Golgi", r"golgi"),
        ("vacuole", r"vacuol|lysosom"), ("extracellular", r"secreted|extracellular|cell wall"),
        ("vesicle", r"vesicle|endosome"), ("cytosol", r"cytoplasm|cytosol"),
        ("membrane", r"plasma membrane|cell membrane|membrane"),
    ]
    fallback = {
        "tc": "ER", "dolichol_pathway": "ER", "er_glycosylation": "ER",
        "folding": "ER", "gpi_biosynthesis": "ER", "copii": "vesicle",
        "copi": "vesicle", "ldsv": "vesicle", "hdsv": "vesicle", "snare": "vesicle",
        "golgi_processing": "Golgi", "cpy_pathway": "vacuole", "alp_pathway": "vacuole",
        "translation": "cytosol", "septin": "cytosol",
        "beta_1_6_glucan_biosynthesis": "membrane",
    }
    for idx in out.index:
        raw_value = out.at[idx, "annotation_compartment_raw"] if "annotation_compartment_raw" in out else None
        raw = "" if pd.isna(raw_value) else str(raw_value)
        hits = []
        for label, pattern in patterns:
            if re.search(pattern, raw, re.I):
                hits.append(label)
                if label == "ER":
                    raw = re.sub(r"endoplasmic reticulum|\ber membrane\b", "", raw, flags=re.I)
        hits = sorted(set(hits))
        if len(hits) == 1:
            out.at[idx, "compartment"] = hits[0]
            out.at[idx, "compartment_source"] = "uniprot"
            out.at[idx, "compartment_confidence"] = "high"
        elif len(hits) > 1:
            out.at[idx, "compartment_source"] = "uniprot_conflict"
        elif out.at[idx, "subsystem"] in fallback:
            out.at[idx, "compartment"] = fallback[out.at[idx, "subsystem"]]
            out.at[idx, "compartment_source"] = "subsystem_inference"
            out.at[idx, "compartment_confidence"] = "medium"
    return out


def assign_glycosylation_roles(
    df: pd.DataFrame, glycosylation_ko: dict[str, list[str]] | None = None,
) -> pd.DataFrame:
    """Assign conservative glycosylation roles with source and confidence."""
    out = df.copy()
    if "glycosylation_role" not in out:
        out["glycosylation_role"] = "none"
    if "glycosylation_role_source" not in out:
        out["glycosylation_role_source"] = None
    if "glycosylation_role_confidence" not in out:
        out["glycosylation_role_confidence"] = None
    direct = {"dolichol_pathway": "n_glycan_assembly", "er_glycosylation": "n_glycan_transfer",
              "gpi_biosynthesis": "gpi_anchor"}
    for idx in out.index:
        if out.at[idx, "glycosylation_role"] != "none":
            continue
        subsystem = out.at[idx, "subsystem"]
        source = out.at[idx, "subsystem_source"]
        if subsystem in direct:
            out.at[idx, "glycosylation_role"] = direct[subsystem]
            out.at[idx, "glycosylation_role_source"] = source
            out.at[idx, "glycosylation_role_confidence"] = "high" if source == "liu2014" else "medium"
            continue
        text = " ".join(str(out.at[idx, c]) for c in ("Description", "annotation_function", "yeast_ortholog")
                        if c in out and pd.notna(out.at[idx, c]))
        rules = [
            ("o_glycosylation", r"protein o-mannosyltransferase|\bpmt\d*\b"),
            ("gpi_anchor", r"gpi[- ]anchor biosynth|gpi mannosyltransferase"),
            ("n_glycan_transfer", r"oligosaccharyltransferase"),
            ("n_glycan_trimming", r"glucosidase i{1,2}|alpha-mannosidase|calnexin|calreticulin"),
            ("golgi_mannosylation", r"\b(?:och|hoc|mnn|ktr)\d*\b|golgi.*mannosyltransferase"),
            ("n_glycan_assembly", r"\b(?:alg|dpm)\d*\b|dolichol.*(?:glycosyl|mannosyl|glucosyl)"),
        ]
        hits = [role for role, pattern in rules if re.search(pattern, text, re.I)]
        if len(set(hits)) == 1:
            out.at[idx, "glycosylation_role"] = hits[0]
            out.at[idx, "glycosylation_role_source"] = "liu_description_or_current_annotation"
            out.at[idx, "glycosylation_role_confidence"] = "medium"
    if glycosylation_ko:
        ko_to_role = {ko: role for role, kos in glycosylation_ko.items() for ko in kos}
        invalid = set(ko_to_role.values()) - {role.value for role in GlycosylationRole}
        if invalid:
            raise ValueError(f"unknown glycosylation roles in KO map: {sorted(invalid)}")
        blank = out["glycosylation_role"].eq("none")
        ko_roles = out["kegg_ko"].map(ko_to_role) if "kegg_ko" in out else pd.Series(index=out.index, dtype=object)
        fill = blank & ko_roles.notna()
        out.loc[fill, "glycosylation_role"] = ko_roles[fill]
        out.loc[fill, "glycosylation_role_source"] = "kegg_ko"
        out.loc[fill, "glycosylation_role_confidence"] = "medium"
    return out


def make_record_ids(df: pd.DataFrame, prefix: str = "AOR") -> pd.DataFrame:
    """Assign a stable surrogate key to every row.

    Deterministic on row position so reruns are reproducible. If you later
    need IDs stable across source revisions, hash the natural key instead.
    """
    out = df.copy()
    out["record_id"] = [f"{prefix}{i:05d}" for i in range(1, len(out) + 1)]
    return out


def build_crosswalk(
    kegg_genes: pd.DataFrame,
    kegg_uniprot: pd.DataFrame,
    kegg_ncbi: pd.DataFrame,
    kegg_ko: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble the identifier lookup, anchored on the locus tag.

    Anchored on KEGG's gene list rather than UniProt because the locus tag is
    our primary biological identifier and KEGG exposes it directly in the
    gene ID. Verify that assumption in notebook 02 before trusting this.
    """
    cw = kegg_genes.copy()

    # 'aor:AO090011000215' -> 'AO090011000215'
    cw["ao_locus_tag"] = cw["kegg_gene_id"].str.split(":").str[-1]

    for other in (kegg_uniprot, kegg_ncbi, kegg_ko):
        if not other.empty:
            cw = cw.merge(other, on="kegg_gene_id", how="left")

    n_before = len(kegg_genes)
    if len(cw) != n_before:
        log.warning(
            "Crosswalk fan-out: %d -> %d rows. A gene maps to multiple "
            "accessions; this is expected but must be handled explicitly.",
            n_before, len(cw),
        )

    log.info(
        "Crosswalk: %d rows | %d with UniProt | %d with KO",
        len(cw),
        cw["uniprot_accession"].notna().sum() if "uniprot_accession" in cw else 0,
        cw["kegg_ko"].notna().sum() if "kegg_ko" in cw else 0,
    )
    return cw


def resolve_by_locus_tag(
    seed: pd.DataFrame,
    crosswalk: pd.DataFrame,
    seed_id_col: str = "liu_ao_locus_tag",
) -> pd.DataFrame:
    """Step 1 of the ladder: direct locus-tag match.

    Returns every seed row with mapping_status set. Nothing is dropped.
    Ambiguous hits (one tag -> several current genes) are flagged rather than
    arbitrarily resolved to the first match.
    """
    s = seed.copy()
    s["_tag"] = s[seed_id_col].astype(str).str.strip().str.upper()

    lookup = crosswalk[["ao_locus_tag", "kegg_gene_id"]].copy()
    if "uniprot_accession" in crosswalk:
        lookup["uniprot_accession"] = crosswalk["uniprot_accession"]
    if "ncbi_gene_id" in crosswalk:
        lookup["ncbi_gene_id"] = crosswalk["ncbi_gene_id"]
    if "kegg_ko" in crosswalk:
        lookup["kegg_ko"] = crosswalk["kegg_ko"]
    lookup["_tag"] = lookup["ao_locus_tag"].astype(str).str.strip().str.upper()

    merged = s.merge(lookup, on="_tag", how="left", indicator=True)

    merged["mapping_status"] = MappingStatus.UNRESOLVED.value
    merged["mapping_method"] = None

    hit = merged["_merge"] == "both"
    merged.loc[hit, "mapping_status"] = MappingStatus.EXACT.value
    merged.loc[hit, "mapping_method"] = "direct_locus_tag"

    # One seed row fanning out to several crosswalk rows. Two real causes:
    # a split gene model, or several protein products (isoforms) for one gene.
    # Either way it needs a human, not an arbitrary pick of the first match.
    # Counted per record_id, not per tag, so duplicate tags in the SOURCE
    # don't produce false positives.
    n_matches = merged.groupby("record_id")["record_id"].transform("size")
    ambiguous = hit & (n_matches > 1)
    merged.loc[ambiguous, "mapping_status"] = MappingStatus.AMBIGUOUS.value
    merged.loc[ambiguous, "mapping_method"] = "ambiguous"

    merged["manual_review_required"] = ambiguous
    merged.loc[~hit, "mapping_method"] = "unresolved"
    merged = merged.drop(columns=["_merge"])

    total = merged["_tag"].nunique()
    unresolved = merged.loc[
        merged["mapping_status"] == MappingStatus.UNRESOLVED.value, "_tag"
    ].nunique()
    log.info(
        "Step 1 (locus tag): %d/%d unresolved (%.1f%%)",
        unresolved, total, 100 * unresolved / total if total else 0,
    )
    return merged


def resolve_by_cross_reference(
    unresolved: pd.DataFrame,
    crosswalk: pd.DataFrame,
    via_col: str,
) -> pd.DataFrame:
    """Step 2: resolve via an existing database cross-reference.

    Only run this on rows step 1 could not place. Set mapping_status to
    cross_reference for anything it recovers.
    """
    raise NotImplementedError(
        "Step 2. Build only if step 1 leaves a meaningful unresolved count - "
        "check the number first rather than writing the whole ladder up front."
    )


def report_unresolved(df: pd.DataFrame) -> pd.DataFrame:
    """Rows needing attention. A deliverable in its own right.

    "Of the N genes in the canonical 2014 list, X no longer map cleanly" is
    genuinely useful to a biologist and surfaces only if someone does the
    plumbing carefully.
    """
    flagged = [
        MappingStatus.UNRESOLVED.value,
        MappingStatus.AMBIGUOUS.value,
        MappingStatus.SPLIT.value,
        MappingStatus.MERGED.value,
    ]
    return df[df["mapping_status"].isin(flagged)].copy()


def assign_missing_subsystems(
    df: pd.DataFrame,
    pathway_members: pd.DataFrame,
    pathway_map: dict[str, str],
) -> pd.DataFrame:
    """Fill blank subsystem labels from KEGG pathway membership."""
    invalid = set(pathway_map.values()) - set(SUBSYSTEM_ORDER)
    if invalid:
        raise ValueError(f"pathway_map contains unknown subsystems: {sorted(invalid)}")

    out = df.copy()
    blank = out["subsystem"].isna() | out["subsystem"].astype(str).str.strip().eq("")
    before = int((~blank).sum())
    if "subsystem_source" not in out:
        out["subsystem_source"] = "unassigned"
    out["subsystem_source"] = out["subsystem_source"].fillna("unassigned")
    out.loc[~blank & out["subsystem_source"].eq("unassigned"), "subsystem_source"] = "liu2014"
    if "manual_review_required" not in out:
        out["manual_review_required"] = False
    out["manual_review_required"] = out["manual_review_required"].fillna(False).astype(bool)

    members = pathway_members[["kegg_gene_id", "kegg_pathway"]].copy()
    members["subsystem_candidate"] = members["kegg_pathway"].map(pathway_map)
    members = members.dropna(subset=["subsystem_candidate"])
    candidates = (
        members.groupby("kegg_gene_id")["subsystem_candidate"]
        .agg(lambda values: sorted(set(values)))
        .to_dict()
    )

    for idx in out.index[blank]:
        choices = candidates.get(out.at[idx, "kegg_gene_id"], [])
        if len(choices) == 1:
            out.at[idx, "subsystem"] = choices[0]
            out.at[idx, "subsystem_source"] = "kegg_pathway"
        elif len(choices) > 1:
            out.at[idx, "subsystem"] = None
            out.at[idx, "subsystem_source"] = "unassigned"
            out.at[idx, "manual_review_required"] = True

    out["pathway_order"] = out["subsystem"].map(SUBSYSTEM_ORDER)
    after = int(out["subsystem"].notna().sum())
    log.info("subsystem coverage: %d/%d -> %d/%d", before, len(out), after, len(out))
    return out


def apply_curated_subsystem_resolutions(
    df: pd.DataFrame, resolutions: dict[str, dict[str, str]],
    id_col: str = "liu_ao_locus_tag",
) -> pd.DataFrame:
    """Apply evidence-reviewed primary subsystems without replacing Liu labels.

    Each resolution requires ``subsystem``, ``source``, ``confidence``, and
    ``rationale``. Rows with a verified Liu subsystem are immutable.
    """
    out = df.copy()
    for field in ("subsystem_confidence", "subsystem_rationale"):
        if field not in out:
            out[field] = None
    for tag, resolution in resolutions.items():
        required = {"subsystem", "source", "confidence", "rationale"}
        missing = required - set(resolution)
        if missing:
            raise ValueError(f"{tag}: missing resolution fields {sorted(missing)}")
        if resolution["subsystem"] not in SUBSYSTEM_ORDER:
            raise ValueError(f"{tag}: unknown subsystem {resolution['subsystem']}")
        mask = out[id_col].eq(tag)
        preserve_liu = not bool(resolution.get("override_verified_liu", False))
        verified_liu = mask & out["subsystem_source"].eq("liu2014") & preserve_liu
        target = mask & ~verified_liu
        out.loc[target, "subsystem"] = resolution["subsystem"]
        out.loc[target, "subsystem_source"] = resolution["source"]
        out.loc[target, "subsystem_confidence"] = resolution["confidence"]
        out.loc[target, "subsystem_rationale"] = resolution["rationale"]
        out.loc[target, "manual_review_required"] = resolution.get("review", False)
    out["pathway_order"] = out["subsystem"].map(SUBSYSTEM_ORDER)
    return out


def assert_no_row_loss(before: pd.DataFrame, after: pd.DataFrame, key: str) -> None:
    """Guard against the failure mode that matters: silent drops.

    A bad join does not raise - it just returns fewer rows, and you get a
    clean-looking result that is quietly missing data.
    """
    lost = set(before[key].astype(str)) - set(after[key].astype(str))
    if lost:
        raise AssertionError(
            f"{len(lost)} rows lost in join. Examples: {sorted(lost)[:5]}"
        )
