"""Schema contract tests.

The schema is a contract. These assert it stays coherent as fields get added.
"""

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
        record_type=RecordType.CLIENT,
        record_origin=RecordOrigin.NEW_ORTHOLOGY,
        mapping_status=MappingStatus.AMBIGUOUS,
        glycosylation_role=GlycosylationRole.GOLGI_MANNOSYLATION,
    )
    row = r.to_row()
    assert row["record_type"] == "secretory_client"
    assert row["record_origin"] == "new_orthology_inference"
    assert row["mapping_status"] == "ambiguous"
    assert row["glycosylation_role"] == "golgi_mannosylation"


def test_subsystem_order_is_monotonic_and_unique():
    """pathway_order drives left-to-right diagram layout - ties would collide."""
    values = list(SUBSYSTEM_ORDER.values())
    assert len(values) == len(set(values)), "duplicate pathway_order values"
    assert values == sorted(values), "SUBSYSTEM_ORDER should read in pathway order"


def test_machinery_and_clients_are_distinct_values():
    assert RecordType.MACHINERY.value != RecordType.CLIENT.value
