import json
from pathlib import Path

from src import config
from src.data import select_complexes as sc

FIXTURES = Path(__file__).parent / "fixtures"


def load_fixture(name):
    with open(FIXTURES / name) as f:
        return json.load(f)


# --- query builder -----------------------------------------------------


def test_build_search_query_uses_config_values_not_literals():
    query = sc.build_search_query()
    nodes = query["query"]["nodes"]
    params_by_attr = {n["parameters"]["attribute"]: n["parameters"] for n in nodes}

    assert (
        params_by_attr["rcsb_entry_info.resolution_combined"]["value"]
        == config.RESOLUTION_CUTOFF_ANGSTROM
    )
    assert (
        params_by_attr["rcsb_accession_info.initial_release_date"]["value"]
        == config.PDB_RELEASE_DATE_CUTOFF
    )
    protein_count_params = params_by_attr["rcsb_entry_info.polymer_entity_count_protein"]
    assert protein_count_params["value"]["from"] == 2
    assert protein_count_params["value"]["to"] == config.MAX_PROTEIN_ENTITIES
    assert (
        params_by_attr["rcsb_entry_info.polymer_entity_count_nucleic_acid"]["value"] == 0
    )
    assert (
        params_by_attr["rcsb_entry_info.polymer_entity_count_nucleic_acid_hybrid"]["value"]
        == 0
    )
    assert params_by_attr["exptl.method"]["value"] == "X-RAY DIFFRACTION"


def test_build_search_query_honors_max_protein_entities_override():
    query = sc.build_search_query(max_protein_entities=4)
    nodes = query["query"]["nodes"]
    params_by_attr = {n["parameters"]["attribute"]: n["parameters"] for n in nodes}
    assert (
        params_by_attr["rcsb_entry_info.polymer_entity_count_protein"]["value"]["to"] == 4
    )


# --- parsing on fixture JSON (no network) -------------------------------


def test_parse_entry_summary():
    entry_data = load_fixture("rcsb_entry.json")
    summary = sc.parse_entry_summary("1ABC", entry_data)

    assert summary["pdb_id"] == "1ABC"
    assert summary["release_date"] == "2019-06-12"
    assert summary["resolution"] == 1.95
    assert summary["polymer_entity_count_protein"] == 2
    assert summary["polymer_entity_ids"] == ["1", "2"]


def test_parse_polymer_entity_normal_produces_one_row_per_chain():
    entity_data = load_fixture("rcsb_polymer_entity_normal.json")
    rows = sc.parse_polymer_entity("1ABC", "1", entity_data)

    assert len(rows) == 2
    chain_ids = {r["chain_id"] for r in rows}
    assert chain_ids == {"A", "C"}
    for row in rows:
        assert row["polymer_type"] == "Protein"
        assert row["seq_length"] == 141
        assert row["uniprot_ids"] == ["P69905"]


# --- SIFTS mapping -------------------------------------------------------


def test_load_sifts_mapping_dedupes_and_uppercases_pdb_id():
    mapping = sc.load_sifts_mapping(FIXTURES / "sifts_sample.tsv")
    assert mapping[("1ABC", "A")] == ("P69905",)
    assert mapping[("1ABC", "B")] == ("P00001",)
    assert mapping[("2DEF", "A")] == ("P99999",)
    assert ("1ABC", "Z") not in mapping


def test_cross_check_uniprot_agreement_and_disagreement():
    mapping = sc.load_sifts_mapping(FIXTURES / "sifts_sample.tsv")

    sifts_ids, agrees = sc.cross_check_uniprot("1ABC", "A", ["P69905"], mapping)
    assert agrees is True
    assert sifts_ids == ("P69905",)

    # RCSB says chimera P11111/P22222, SIFTS says a single different
    # accession P99999 for 2DEF chain A -> a genuine disagreement.
    sifts_ids, agrees = sc.cross_check_uniprot(
        "2DEF", "A", ["P11111", "P22222"], mapping
    )
    assert agrees is False
    assert sifts_ids == ("P99999",)


# --- chain-level filters -------------------------------------------------


def _row(pdb_id="1ABC", entity_id="1", chain_id="A", polymer_type="Protein",
         seq_length=100, uniprot_ids=None):
    return {
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "chain_id": chain_id,
        "polymer_type": polymer_type,
        "seq_length": seq_length,
        "uniprot_ids": uniprot_ids if uniprot_ids is not None else ["P00001"],
    }


def test_filter_non_protein_entities_drops_and_logs_reason():
    rows = [_row(chain_id="A", polymer_type="Protein"), _row(chain_id="D", polymer_type="Other")]
    kept, dropped = sc.filter_non_protein_entities(rows)
    assert [r["chain_id"] for r in kept] == ["A"]
    assert len(dropped) == 1
    assert dropped[0]["chain_id"] == "D"
    assert "non_protein_entity_type" in dropped[0]["drop_reason"]


def test_filter_min_chain_length_drops_short_chain():
    rows = [_row(chain_id="A", seq_length=141), _row(chain_id="B", seq_length=12)]
    kept, dropped = sc.filter_min_chain_length(rows, config.MIN_CHAIN_LENGTH)
    assert [r["chain_id"] for r in kept] == ["A"]
    assert len(dropped) == 1
    assert dropped[0]["chain_id"] == "B"
    assert "chain_length" in dropped[0]["drop_reason"]


def test_filter_uniprot_mapping_drops_unmapped_chain():
    rows = [_row(chain_id="A", uniprot_ids=["P69905"]), _row(chain_id="D", uniprot_ids=[])]
    kept, dropped = sc.filter_uniprot_mapping(rows)
    assert [r["chain_id"] for r in kept] == ["A"]
    assert dropped[0]["drop_reason"] == "no_uniprot_mapping"


def test_filter_chimera_drops_multi_uniprot_chain():
    rows = [
        _row(chain_id="A", uniprot_ids=["P69905"]),
        _row(chain_id="X", uniprot_ids=["P11111", "P22222"]),
    ]
    kept, dropped = sc.filter_chimera(rows)
    assert [r["chain_id"] for r in kept] == ["A"]
    assert len(dropped) == 1
    assert dropped[0]["chain_id"] == "X"
    assert dropped[0]["drop_reason"] == "chimera:P11111,P22222"


def test_run_filter_pipeline_combines_all_filters_and_tracks_attrition():
    rows = [
        _row(pdb_id="1ABC", chain_id="A", polymer_type="Protein", seq_length=141, uniprot_ids=["P69905"]),
        _row(pdb_id="1ABC", chain_id="D", polymer_type="Other", seq_length=50, uniprot_ids=[]),
        _row(pdb_id="1ABC", chain_id="B", polymer_type="Protein", seq_length=12, uniprot_ids=["P00001"]),
        _row(pdb_id="3GHI", chain_id="A", polymer_type="Protein", seq_length=180, uniprot_ids=[]),
        _row(pdb_id="2DEF", chain_id="A", polymer_type="Protein", seq_length=220, uniprot_ids=["P11111", "P22222"]),
    ]

    survivors, attrition = sc.run_filter_pipeline(rows, config.MIN_CHAIN_LENGTH)

    assert [r["chain_id"] for r in survivors] == ["A"]
    assert survivors[0]["pdb_id"] == "1ABC"

    stage_names = [a["stage"] for a in attrition]
    assert stage_names[0] == "initial_protein_chains"
    assert "non_protein_entity_type" in stage_names
    assert any(s.startswith("min_chain_length") for s in stage_names)
    assert "uniprot_mapping" in stage_names
    assert "chimera" in stage_names

    # every stage's n_remaining should equal the next stage's n_before
    for i in range(len(attrition) - 1):
        assert attrition[i]["n_remaining"] == attrition[i + 1]["n_before"]

    final = attrition[-1]
    assert final["n_remaining"] == 1
