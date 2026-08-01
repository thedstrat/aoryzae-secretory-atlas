from pathlib import Path

import pandas as pd

from atlas.crosswalk_audit import AUDIT_COLUMNS, classify_mapping


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_known_exact_mapping_example():
    result = classify_mapping(
        "AO090011000795", "AO090011000795",
        {"NCBI RefSeq": {"AO090011000795"}, "KEGG": {"AO090011000795"}},
    )
    assert result["audit_status"] == "confirmed"
    assert result["mapping_method"] == "unchanged/exact identifier"


def test_known_ambiguous_mapping_pattern():
    result = classify_mapping("OLD", None, {"source": {"NEW1", "NEW2"}})
    assert result["audit_status"] == "ambiguous"
    assert result["mapping_method"] == "ambiguous"


def test_known_split_mapping_pattern():
    result = classify_mapping(
        "OLD", None, {"authoritative annotation": {"NEW1", "NEW2"}},
        relationship_hint="split",
    )
    assert result["audit_status"] == "ambiguous"
    assert result["mapping_method"] == "split"


def test_known_merged_mapping_pattern():
    result = classify_mapping(
        "OLD1", "NEW", {"NCBI RefSeq": {"NEW"}, "KEGG": {"NEW"}},
        reverse_proposed_count=2,
    )
    assert result["audit_status"] == "confirmed"
    assert result["mapping_method"] == "merged"


def test_known_unresolved_mapping_example():
    result = classify_mapping("AO090005001666", None, {})
    assert result["audit_status"] == "unresolved"
    assert result["mapping_method"] == "unresolved"


def test_conflicting_sources_are_not_silently_resolved():
    result = classify_mapping("OLD", "NEW1", {"NCBI RefSeq": {"NEW1"}, "KEGG": {"NEW2"}})
    assert result["audit_status"] == "conflicting"
    assert result["manual_review_required"] is True


def test_committed_crosswalk_audit_contract():
    audit = pd.read_csv(PROJECT_ROOT / "data/processed/crosswalk_audit.csv")
    assert audit.columns.tolist() == AUDIT_COLUMNS
    assert len(audit) == 369
    assert audit.original_id.is_unique
    assert set(audit.audit_status) <= {"confirmed", "probable", "ambiguous", "conflicting", "unresolved"}

    exact = audit.loc[audit.original_id.eq("AO090011000795")].iloc[0]
    assert exact.audit_status == "confirmed"
    assert exact.verified_current_id == "AO090011000795"

    unresolved = audit.loc[audit.original_id.eq("AO090005001666")].iloc[0]
    assert unresolved.audit_status == "unresolved"
    assert bool(unresolved.manual_review_required)
