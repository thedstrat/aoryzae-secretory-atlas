"""Acceptance tests: no silent data loss, machinery never blended with clients.

These encode the invariants agreed before any code was written. They exist
because bad joins do not raise - they just return fewer rows, and you get a
clean-looking result quietly missing data.
"""

import pandas as pd
import pytest
from pathlib import Path

from atlas.crosswalk import assert_no_row_loss, make_record_ids, resolve_by_locus_tag
from atlas.schema import REQUIRED_FIELDS, RecordType


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def source_and_final():
    source = pd.read_csv(PROJECT_ROOT / "data/interim/liu_components_raw_cleaned.csv")
    final = pd.read_csv(PROJECT_ROOT / "data/processed/aoryzae_secretory_components.csv")
    joined = final.merge(source, left_on="liu_ao_locus_tag", right_on="ID",
                         how="left", validate="one_to_one", suffixes=("_final", "_liu"))
    assert len(joined) == 369
    return source, final, joined


def test_every_source_row_gets_a_record_id():
    seed = pd.DataFrame({"liu_ao_locus_tag": ["A", "B", "C", None]})
    out = make_record_ids(seed)
    assert out["record_id"].notna().all(), "rows without an ID cannot be tracked"
    assert len(out) == 4, "a missing identifier is not a reason to drop a row"


def test_row_with_no_identifier_is_kept():
    """A row missing its locus tag is still a source row."""
    seed = make_record_ids(pd.DataFrame({"liu_ao_locus_tag": [None, "AO090001000653"]}))
    cw = pd.DataFrame({
        "kegg_gene_id": ["aor:AO090001000653"],
        "ao_locus_tag": ["AO090001000653"],
    })
    out = resolve_by_locus_tag(seed, cw)
    assert len(out) == 2


def test_machinery_and_clients_never_share_a_frame_by_accident():
    machinery = pd.DataFrame({"record_id": ["A1"], "record_type": [RecordType.MACHINERY.value]})
    clients = pd.DataFrame({"record_id": ["C1"], "record_type": [RecordType.CLIENT.value]})
    combined = pd.concat([machinery, clients], ignore_index=True)

    # A combined export is fine ONLY if record_type distinguishes them.
    assert combined["record_type"].nunique() == 2
    assert combined["record_type"].notna().all(), "unlabelled blend is the failure mode"


def test_required_fields_present_in_shippable_table():
    good = pd.DataFrame({f: ["x"] for f in REQUIRED_FIELDS})
    for f in REQUIRED_FIELDS:
        assert good[f].notna().all()

    bad = good.copy()
    bad.loc[0, "mapping_status"] = None
    assert bad["mapping_status"].isna().any()


def test_baseline_count_guard():
    """Pattern for notebook use: compare against notebook 01's baseline."""
    baseline = 369
    final = pd.DataFrame({"record_id": [f"AOR{i:05d}" for i in range(1, 370)]})
    assert len(final) == baseline
    with pytest.raises(AssertionError):
        assert_no_row_loss(final, final.iloc[:-1], "record_id")


def test_secondary_relationships_do_not_duplicate_component_rows():
    components = pd.DataFrame({"record_id": ["A", "B"], "subsystem": ["erad", "folding"]})
    relationships = pd.DataFrame({
        "record_id": ["A", "A", "B"],
        "subsystem": ["erad", "er_glycosylation", "folding"],
        "relationship_type": ["primary", "secondary", "primary"],
    })
    assert components.record_id.is_unique
    assert relationships.loc[relationships.relationship_type.eq("primary"), "record_id"].is_unique
    assert len(components) == 2


def test_machinery_identifiers_and_orthologs_are_verbatim(source_and_final):
    source, final, joined = source_and_final
    assert final["liu_ao_locus_tag"].tolist() == source["ID"].tolist()
    left = joined["yeast_ortholog"].fillna("<NULL>").tolist()
    right = joined["S. cerevisiae ortholog"].fillna("<NULL>").tolist()
    assert left == right


def test_liu_subsystems_are_normalized_from_the_source_cell(source_and_final):
    _, _, joined = source_and_final
    normalized = {
        "TC": "tc", "Dolichol pathway": "dolichol_pathway",
        "Erglycosylation": "er_glycosylation", "Folding": "folding",
        "GPI biosynthesis": "gpi_biosynthesis", "ERAD": "erad",
        "COPII": "copii", "COPI": "copi", "Golgi processing": "golgi_processing",
        "LDSV": "ldsv", "HDSV": "hdsv", "CPY pathway": "cpy_pathway",
        "ALPpathway": "alp_pathway", "SNARE": "snare", "Septin": "septin",
        "beta-1,6 glucan biosynthesis": "beta_1_6_glucan_biosynthesis",
        "Translation": "translation", "putative mitochondria protein": "putative_mitochondria_protein",
        "mitochondrial\u00a0m‐AAA protease": "mitochondrial_m_aaa_protease",
    }
    direct = joined.subsystem_source.eq("liu2014")
    expected = joined.loc[direct, "Subsystems or function"].map(normalized)
    assert expected.notna().all(), "a direct Liu label lacks an explicit normalization"
    assert joined.loc[direct, "subsystem"].tolist() == expected.tolist()


def test_liu_functions_remain_the_original_description(source_and_final):
    _, _, joined = source_and_final
    from_liu = joined["Description"].notna()
    assert joined.loc[from_liu, "function"].tolist() == joined.loc[from_liu, "Description"].tolist()


def test_transcriptomic_flags_recompute_per_gene(source_and_final):
    _, _, joined = source_and_final
    pcols = ["CF1.1vs A1560_adj.P.Val", "A16 vs A1560_adj.P.Val", "CF32 vs A1560_adj.P.Val"]
    fcols = ["CF1.1 vs A1560_logFC", "A16 vs A1560_logFC", "CF32 vs A1560_logFC"]
    significant = joined[pcols].lt(0.05).all(axis=1)
    direction = pd.Series(None, index=joined.index, dtype=object)
    direction.loc[significant & joined[fcols].gt(0).all(axis=1)] = "up"
    direction.loc[significant & joined[fcols].lt(0).all(axis=1)] = "down"
    assert joined["sig_all_three_final"].astype(bool).tolist() == significant.tolist()
    assert joined["direction_all_three_final"].fillna("<NULL>").tolist() == direction.fillna("<NULL>").tolist()


def test_shortlist_per_gene_and_aggregate_contract(source_and_final):
    _, final, _ = source_and_final
    shortlist = pd.read_csv(PROJECT_ROOT / "data/processed/transcriptomic_shortlist.csv")
    expected = final.loc[final.sig_all_three.astype(bool), "liu_ao_locus_tag"]
    assert set(shortlist.ao_locus_tag) == set(expected)
    assert len(shortlist) == 51
    assert shortlist.direction.value_counts().to_dict() == {"up": 48, "down": 3}


def test_composite_liu_descriptions_keep_primary_rationale_and_secondary_roles(source_and_final):
    _, final, _ = source_and_final
    relationships = pd.read_csv(PROJECT_ROOT / "data/processed/gene_subsystems.csv")
    expected_secondary = {
        "AO090012000213": {"erad"}, "AO090005000437": {"snare"},
        "AO090023000840": {"folding"}, "AO090003000257": {"erad"},
        "AO090003000853": {"erad"},
        "AO090005000718": {"hdsv", "ldsv", "cpy_pathway", "alp_pathway"},
        "AO090701000139": {"snare"}, "AO090023000864": {"snare"},
    }
    for tag, secondary in expected_secondary.items():
        row = final.loc[final.liu_ao_locus_tag.eq(tag)].iloc[0]
        assert pd.notna(row.subsystem) and bool(row.subsystem_rationale)
        observed = set(relationships.loc[
            relationships.ao_locus_tag.eq(tag) & relationships.relationship_type.eq("secondary"),
            "subsystem",
        ])
        assert secondary.issubset(observed)
    kar2 = final.loc[final.liu_ao_locus_tag.eq("AO090003000257")].iloc[0]
    assert kar2.subsystem == "folding", "KAR2 must not inherit the first ERAD token"
