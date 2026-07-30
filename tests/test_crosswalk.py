"""Crosswalk tests - all failure paths.

The happy path is obvious. What bites you is a locus tag silently matching
two current genes, or unresolved rows vanishing in a join.
"""

import pandas as pd
import pytest

from atlas.crosswalk import (
    assert_no_row_loss, assign_controlled_compartments, assign_glycosylation_roles,
    apply_curated_subsystem_resolutions, assign_missing_subsystems,
    assign_subsystems_from_evidence, build_crosswalk,
    make_record_ids, report_unresolved, resolve_by_locus_tag,
)
from atlas.schema import MappingStatus


@pytest.fixture
def crosswalk_df():
    return pd.DataFrame({
        "kegg_gene_id": ["aor:AO090001000653", "aor:AO090011000215"],
        "ao_locus_tag": ["AO090001000653", "AO090011000215"],
        "uniprot_accession": ["Q2UNU5", "Q2UNU6"],
        "ncbi_gene_id": ["5993412", "5993413"],
        "kegg_ko": ["ko:K09027", "ko:K05528"],
    })


def test_record_ids_are_unique_and_stable():
    df = pd.DataFrame({"x": range(5)})
    out = make_record_ids(df)
    assert out["record_id"].nunique() == 5
    assert make_record_ids(df)["record_id"].tolist() == out["record_id"].tolist()


def test_unresolved_rows_survive(crosswalk_df):
    seed = make_record_ids(pd.DataFrame({
        "liu_ao_locus_tag": ["AO090001000653", "AO999999999999"]
    }))
    out = resolve_by_locus_tag(seed, crosswalk_df)

    assert len(out) == 2, "unresolved rows must not be dropped"
    assert MappingStatus.UNRESOLVED.value in set(out["mapping_status"])
    assert MappingStatus.EXACT.value in set(out["mapping_status"])


def test_ambiguous_mapping_flagged_not_silently_resolved():
    """One seed row fanning out to several records needs a human.

    Real causes: a split gene model, or several protein products for one gene.
    Either way, arbitrarily taking the first match loses information.
    """
    seed = make_record_ids(pd.DataFrame({"liu_ao_locus_tag": ["AO090001000653"]}))
    cw = pd.DataFrame({
        "kegg_gene_id": ["aor:AO090001000653", "aor:AO090001000653"],
        "ao_locus_tag": ["AO090001000653", "AO090001000653"],
        "uniprot_accession": ["Q1", "Q2"],   # two products, one gene
    })
    out = resolve_by_locus_tag(seed, cw)

    assert all(out["mapping_status"] == MappingStatus.AMBIGUOUS.value)
    assert all(out["manual_review_required"])
    assert len(out) == 2, "both candidate records must be retained"


def test_source_duplicates_do_not_false_positive_as_ambiguous():
    """Two seed rows with the same tag is a source-data issue, not ambiguity."""
    seed = make_record_ids(pd.DataFrame({
        "liu_ao_locus_tag": ["AO090001000653", "AO090001000653"]
    }))
    cw = pd.DataFrame({
        "kegg_gene_id": ["aor:AO090001000653"],
        "ao_locus_tag": ["AO090001000653"],
        "uniprot_accession": ["Q1"],
    })
    out = resolve_by_locus_tag(seed, cw)

    assert all(out["mapping_status"] == MappingStatus.EXACT.value)
    assert not any(out["manual_review_required"])


def test_case_and_whitespace_insensitive(crosswalk_df):
    seed = make_record_ids(pd.DataFrame({
        "liu_ao_locus_tag": ["  ao090001000653  "]
    }))
    out = resolve_by_locus_tag(seed, crosswalk_df)
    assert out["mapping_status"].iloc[0] == MappingStatus.EXACT.value


def test_report_unresolved_includes_ambiguous():
    df = pd.DataFrame({"mapping_status": [
        MappingStatus.EXACT.value,
        MappingStatus.UNRESOLVED.value,
        MappingStatus.AMBIGUOUS.value,
        MappingStatus.SPLIT.value,
    ]})
    assert len(report_unresolved(df)) == 3


def test_assert_no_row_loss_raises_on_drop():
    before = pd.DataFrame({"record_id": ["A", "B", "C"]})
    after = pd.DataFrame({"record_id": ["A", "B"]})
    with pytest.raises(AssertionError, match="rows lost"):
        assert_no_row_loss(before, after, "record_id")


def test_build_crosswalk_extracts_locus_tag():
    genes = pd.DataFrame({
        "kegg_gene_id": ["aor:AO090001000653"], "kegg_description": ["x"]
    })
    empty = pd.DataFrame()
    out = build_crosswalk(genes, empty, empty, empty)
    assert out["ao_locus_tag"].iloc[0] == "AO090001000653"


def test_assign_missing_subsystems_never_overwrites_liu_label():
    genes = pd.DataFrame({
        "kegg_gene_id": ["aor:A"], "subsystem": ["folding"],
        "manual_review_required": [False],
    })
    members = pd.DataFrame({
        "kegg_gene_id": ["aor:A"], "kegg_pathway": ["aor03060"],
    })
    out = assign_missing_subsystems(genes, members, {"aor03060": "tc"})
    assert out.loc[0, "subsystem"] == "folding"
    assert out.loc[0, "subsystem_source"] == "liu2014"


def test_assign_missing_subsystems_flags_conflicting_pathways():
    genes = pd.DataFrame({
        "kegg_gene_id": ["aor:A"], "subsystem": [None],
        "manual_review_required": [False],
    })
    members = pd.DataFrame({
        "kegg_gene_id": ["aor:A", "aor:A"],
        "kegg_pathway": ["aor03060", "aor04130"],
    })
    out = assign_missing_subsystems(
        genes, members, {"aor03060": "tc", "aor04130": "snare"}
    )
    assert pd.isna(out.loc[0, "subsystem"])
    assert out.loc[0, "subsystem_source"] == "unassigned"
    assert bool(out.loc[0, "manual_review_required"])


def test_assign_missing_subsystems_retains_unmapped_gene():
    genes = pd.DataFrame({
        "kegg_gene_id": ["aor:A"], "subsystem": [None],
        "manual_review_required": [False],
    })
    members = pd.DataFrame(columns=["kegg_gene_id", "kegg_pathway"])
    out = assign_missing_subsystems(genes, members, {"aor03060": "tc"})
    assert len(out) == 1
    assert pd.isna(out.loc[0, "subsystem"])
    assert out.loc[0, "subsystem_source"] == "unassigned"


def test_evidence_assignment_preserves_liu_and_flags_description_conflict():
    genes = pd.DataFrame({
        "subsystem": ["folding", None], "subsystem_source": ["unassigned", "unassigned"],
        "Description": ["COPII", "CPY pathway, HDSV"],
        "yeast_ortholog": ["X", "Y"], "annotation_function": [None, None],
        "manual_review_required": [False, False],
    })
    scaffold = pd.DataFrame({"yeast_ortholog": ["X"], "subsystem": ["copii"]})
    out = assign_subsystems_from_evidence(genes, scaffold)
    assert out.loc[0, "subsystem"] == "folding"
    assert out.loc[0, "subsystem_source"] == "liu2014"
    assert pd.isna(out.loc[1, "subsystem"])
    assert bool(out.loc[1, "manual_review_required"])


def test_controlled_compartment_uses_direct_text_before_subsystem():
    genes = pd.DataFrame({"annotation_compartment_raw": ["SUBCELLULAR LOCATION: Golgi apparatus."],
                          "subsystem": ["folding"], "manual_review_required": [False]})
    out = assign_controlled_compartments(genes)
    assert out.loc[0, "compartment"] == "Golgi"
    assert out.loc[0, "compartment_source"] == "uniprot"


def test_glycosylation_role_has_source_and_confidence():
    genes = pd.DataFrame({"subsystem": ["gpi_biosynthesis"], "subsystem_source": ["yeast_scaffold"],
                          "Description": [""], "annotation_function": [""], "yeast_ortholog": ["GPI1"]})
    out = assign_glycosylation_roles(genes)
    assert out.loc[0, "glycosylation_role"] == "gpi_anchor"
    assert out.loc[0, "glycosylation_role_source"] == "yeast_scaffold"
    assert out.loc[0, "glycosylation_role_confidence"] == "medium"


def test_kegg_ko_fills_blank_role_but_never_overwrites_liu():
    genes = pd.DataFrame({
        "subsystem": [None, None], "subsystem_source": ["unassigned", "unassigned"],
        "Description": ["", ""], "annotation_function": ["", ""],
        "yeast_ortholog": ["X", "Y"], "kegg_ko": ["ko:K05528", "ko:K05528"],
        "glycosylation_role": ["none", "n_glycan_transfer"],
        "glycosylation_role_source": [None, "liu2014"],
        "glycosylation_role_confidence": [None, "high"],
    })
    out = assign_glycosylation_roles(genes, {"golgi_mannosylation": ["ko:K05528"]})
    assert out.loc[0, "glycosylation_role"] == "golgi_mannosylation"
    assert out.loc[0, "glycosylation_role_source"] == "kegg_ko"
    assert out.loc[0, "glycosylation_role_confidence"] == "medium"
    assert out.loc[1, "glycosylation_role"] == "n_glycan_transfer"
    assert out.loc[1, "glycosylation_role_source"] == "liu2014"


def test_multi_pathway_membership_does_not_null_existing_primary():
    genes = pd.DataFrame({"kegg_gene_id": ["aor:A"], "subsystem": ["er_glycosylation"],
                          "subsystem_source": ["curated_evidence"],
                          "manual_review_required": [False]})
    members = pd.DataFrame({"kegg_gene_id": ["aor:A", "aor:A"],
                            "kegg_pathway": ["aor04141", "aor00510"]})
    out = assign_missing_subsystems(
        genes, members, {"aor04141": "folding", "aor00510": "er_glycosylation"}
    )
    assert out.loc[0, "subsystem"] == "er_glycosylation"
    assert out.loc[0, "subsystem_source"] == "curated_evidence"


def test_curated_primary_requires_and_records_provenance():
    genes = pd.DataFrame({"liu_ao_locus_tag": ["A"], "subsystem": [None],
                          "subsystem_source": ["unassigned"],
                          "manual_review_required": [True]})
    resolution = {"A": {"subsystem": "erad", "source": "direct_annotation",
                         "confidence": "high", "rationale": "Specific ERAD function."}}
    out = apply_curated_subsystem_resolutions(genes, resolution)
    assert out.loc[0, "subsystem_source"] == "direct_annotation"
    assert out.loc[0, "subsystem_confidence"] == "high"
    assert out.loc[0, "subsystem_rationale"] == "Specific ERAD function."


def test_curated_resolution_never_overwrites_verified_liu():
    genes = pd.DataFrame({"liu_ao_locus_tag": ["A"], "subsystem": ["folding"],
                          "subsystem_source": ["liu2014"],
                          "manual_review_required": [False]})
    resolution = {"A": {"subsystem": "erad", "source": "direct_annotation",
                         "confidence": "high", "rationale": "candidate"}}
    out = apply_curated_subsystem_resolutions(genes, resolution)
    assert out.loc[0, "subsystem"] == "folding"
    assert out.loc[0, "subsystem_source"] == "liu2014"
