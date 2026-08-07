"""Schema contract tests.

The schema is a contract. These assert it stays coherent as fields get added.
"""

from pathlib import Path

import pandas as pd

from atlas.schema import (
    COLUMN_ORDER, DATA_DICTIONARY, REQUIRED_FIELDS, SUBSYSTEM_ORDER,
    GeneRecord, GlycosylationRole, MappingStatus, RecordOrigin, RecordType,
)


def test_every_column_has_a_definition():
    """A field with no data dictionary entry is a field nobody can interpret."""
    undocumented = set(COLUMN_ORDER) - set(DATA_DICTIONARY)
    assert not undocumented, f"undocumented fields: {sorted(undocumented)}"


def test_no_orphan_definitions():
    orphans = set(DATA_DICTIONARY) - set(COLUMN_ORDER)
    assert not orphans, f"documented but not in COLUMN_ORDER: {sorted(orphans)}"


def test_data_dictionary_has_structured_provenance():
    required = {"definition", "allowed_values", "value_meanings", "blank_means", "source"}
    for name, entry in DATA_DICTIONARY.items():
        assert set(entry) == required, f"{name}: incomplete dictionary metadata"
        assert entry["source"] in {"Liu 2014", "fetched", "generated"}
        assert all(str(entry[key]).strip() for key in {"definition", "allowed_values", "blank_means", "source"})


def test_dictionary_explains_confidence_mapping_and_blank_roles():
    assert "high =" in DATA_DICTIONARY["subsystem_confidence"]["value_meanings"]
    assert "medium =" in DATA_DICTIONARY["glycosylation_role_confidence"]["value_meanings"]
    assert "glycosylation_role is none" in DATA_DICTIONARY["glycosylation_role_source"]["blank_means"]
    assert "exact =" in DATA_DICTIONARY["mapping_status"]["value_meanings"]


def test_exported_column_descriptions_match_schema():
    root = Path(__file__).resolve().parents[1]
    exported = pd.read_csv(root / "data/processed/column_descriptions.csv").fillna("")
    expected = pd.DataFrame([{"field": field, **entry} for field, entry in DATA_DICTIONARY.items()]).fillna("")
    pd.testing.assert_frame_equal(exported, expected, check_dtype=False)


def test_record_serialises_to_every_column():
    r = GeneRecord(record_id="AOR00001")
    row = r.to_row()
    missing = set(COLUMN_ORDER) - set(row)
    assert not missing, f"COLUMN_ORDER references non-fields: {sorted(missing)}"


def test_required_fields_are_real_columns():
    assert set(REQUIRED_FIELDS).issubset(COLUMN_ORDER)


def test_enums_serialise_to_strings():
    r = GeneRecord(
        record_id="AOR00001",
        record_type=RecordType.MACHINERY,
        record_origin=RecordOrigin.NEW_ORTHOLOGY,
        mapping_status=MappingStatus.AMBIGUOUS,
        glycosylation_role=GlycosylationRole.GOLGI_MANNOSYLATION,
    )
    row = r.to_row()
    assert row["record_type"] == "machinery_component"
    assert row["record_origin"] == "new_orthology_inference"
    assert row["mapping_status"] == "ambiguous"
    assert row["glycosylation_role"] == "golgi_mannosylation"


def test_subsystem_order_is_monotonic_and_unique():
    """pathway_order drives left-to-right diagram layout - ties would collide."""
    values = list(SUBSYSTEM_ORDER.values())
    assert len(values) == len(set(values)), "duplicate pathway_order values"
    assert values == sorted(values), "SUBSYSTEM_ORDER should read in pathway order"
