import json
from pathlib import Path

import pytest

from src import config
from src.data import select_complexes as sc

FIXTURES = Path(__file__).parent / "fixtures"


def load_fixture(name):
    with open(FIXTURES / name) as f:
        return json.load(f)


class FakeResponse:
    def __init__(self, json_data, status_code=200):
        self._json_data = json_data
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"fake HTTP error {self.status_code}")

    def json(self):
        return self._json_data


class FakeSession:
    """Minimal requests.Session stand-in that records calls and fails loudly
    on any call not explicitly wired up, so resumability tests can assert
    "no network happened" by simply not wiring anything."""

    def __init__(self, post_response=None, get_responses=None):
        self._post_response = post_response
        self._get_responses = get_responses or {}
        self.post_calls = []
        self.get_calls = []

    def post(self, url, json=None, headers=None):
        self.post_calls.append((url, json))
        if self._post_response is None:
            raise AssertionError(f"unexpected POST to {url} (test wired no response)")
        return FakeResponse(self._post_response)

    def get(self, url, stream=False):
        self.get_calls.append(url)
        if url not in self._get_responses:
            raise AssertionError(f"unexpected GET to {url} (test wired no response)")
        return FakeResponse(self._get_responses[url])


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


# --- deposition-group collapse (pure function, no I/O) -------------------


def test_collapse_deposition_groups_keeps_best_resolution_member():
    group_info = {
        "5QY1": {"resolution": 1.8, "group_id": "G_1"},
        "5QY2": {"resolution": 1.2, "group_id": "G_1"},
        "5QY3": {"resolution": 2.1, "group_id": "G_1"},
    }
    kept, dropped = sc.collapse_deposition_groups(group_info)

    assert kept == ["5QY2"]
    assert len(dropped) == 2
    kept_reps = {d["kept_representative"] for d in dropped}
    assert kept_reps == {"5QY2"}
    dropped_ids = {d["pdb_id"] for d in dropped}
    assert dropped_ids == {"5QY1", "5QY3"}


def test_collapse_deposition_groups_keeps_ungrouped_entries_as_singletons():
    group_info = {
        "4HHB": {"resolution": 1.74, "group_id": None},
        "10CY": {"resolution": 2.49, "group_id": None},
    }
    kept, dropped = sc.collapse_deposition_groups(group_info)

    assert set(kept) == {"4HHB", "10CY"}
    assert dropped == []


def test_collapse_deposition_groups_breaks_resolution_tie_with_entry_id():
    group_info = {
        "5QY9": {"resolution": 1.5, "group_id": "G_2"},
        "5QY2": {"resolution": 1.5, "group_id": "G_2"},
    }
    kept, dropped = sc.collapse_deposition_groups(group_info)
    assert kept == ["5QY2"]
    assert dropped[0]["pdb_id"] == "5QY9"


def test_collapse_deposition_groups_missing_resolution_sorts_last():
    group_info = {
        "5QY1": {"resolution": None, "group_id": "G_3"},
        "5QY2": {"resolution": 2.4, "group_id": "G_3"},
    }
    kept, dropped = sc.collapse_deposition_groups(group_info)
    assert kept == ["5QY2"]
    assert dropped[0]["pdb_id"] == "5QY1"


def test_collapse_deposition_groups_mixed_grouped_and_singleton():
    group_info = {
        "5QY1": {"resolution": 1.5, "group_id": "G_1"},
        "5QY2": {"resolution": 1.2, "group_id": "G_1"},
        "4HHB": {"resolution": 1.74, "group_id": None},
    }
    kept, dropped = sc.collapse_deposition_groups(group_info)
    assert set(kept) == {"5QY2", "4HHB"}
    assert len(dropped) == 1
    assert dropped[0]["pdb_id"] == "5QY1"


# --- group-membership batched fetch ---------------------------------------


def test_fetch_group_membership_batched_parses_resolution_and_filters_aggregation_method(monkeypatch, tmp_path):
    monkeypatch.setattr(sc, "GRAPHQL_GROUP_BATCH_CACHE_DIR", tmp_path / "graphql_group_batches")
    canned = {
        "data": {
            "entries": [
                {
                    "rcsb_id": "5QY1",
                    "rcsb_entry_info": {"resolution_combined": [1.8]},
                    "rcsb_entry_group_membership": [
                        {"group_id": "G_1002115", "aggregation_method": "matching_deposit_group_id"}
                    ],
                },
                {
                    "rcsb_id": "4HHB",
                    "rcsb_entry_info": {"resolution_combined": [1.74]},
                    # A sequence-similarity group, NOT a true deposit group --
                    # must be ignored, leaving group_id=None.
                    "rcsb_entry_group_membership": [
                        {"group_id": "G_seq_9", "aggregation_method": "sequence_identity"}
                    ],
                },
                {
                    "rcsb_id": "10CY",
                    "rcsb_entry_info": {"resolution_combined": [2.49]},
                    "rcsb_entry_group_membership": None,
                },
            ]
        }
    }
    session = FakeSession(post_response=canned)

    info = sc.fetch_group_membership_batched(session, ["5QY1", "4HHB", "10CY"])

    assert info["5QY1"] == {"resolution": 1.8, "group_id": "G_1002115"}
    assert info["4HHB"] == {"resolution": 1.74, "group_id": None}
    assert info["10CY"] == {"resolution": 2.49, "group_id": None}


# --- GraphQL batch response splitting (pure function, no I/O) ------------


def test_split_graphql_entry_matches_rest_fixture_shape_for_parsing():
    batch = load_fixture("rcsb_graphql_entry_batch.json")
    entry = batch["data"]["entries"][0]

    entry_json, entities = sc.split_graphql_entry(entry)

    assert "rcsb_id" not in entry_json
    assert "polymer_entities" not in entry_json
    assert set(entities.keys()) == {"1"}
    assert "rcsb_id" not in entities["1"]

    # The split output must be directly consumable by the same parsers the
    # REST-fixture tests above exercise.
    summary = sc.parse_entry_summary("1ABC", entry_json)
    assert summary["pdb_id"] == "1ABC"
    assert summary["release_date"] == "2019-06-12"
    assert summary["resolution"] == 1.95
    assert summary["polymer_entity_ids"] == ["1"]

    rows = sc.parse_polymer_entity("1ABC", "1", entities["1"])
    assert len(rows) == 2
    assert {r["chain_id"] for r in rows} == {"A", "C"}
    assert rows[0]["uniprot_ids"] == ["P69905"]


# --- warm_entry_cache_via_graphql (resumability + fallback) --------------


def test_warm_entry_cache_via_graphql_skips_entries_already_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(sc, "ENTRY_CACHE_DIR", tmp_path / "entries")
    monkeypatch.setattr(sc, "POLYMER_ENTITY_CACHE_DIR", tmp_path / "polymer_entities")
    monkeypatch.setattr(sc, "GRAPHQL_ENTRY_BATCH_CACHE_DIR", tmp_path / "graphql_entry_batches")
    (tmp_path / "entries").mkdir(parents=True)
    (tmp_path / "entries" / "1ABC.json").write_text('{"already": "cached"}')

    session = FakeSession()  # wired with nothing -> raises if any network call is made
    sc.warm_entry_cache_via_graphql(session, ["1ABC"])

    assert session.post_calls == []


def test_warm_entry_cache_via_graphql_writes_entry_and_entity_cache_files(tmp_path, monkeypatch):
    monkeypatch.setattr(sc, "ENTRY_CACHE_DIR", tmp_path / "entries")
    monkeypatch.setattr(sc, "POLYMER_ENTITY_CACHE_DIR", tmp_path / "polymer_entities")
    monkeypatch.setattr(sc, "GRAPHQL_ENTRY_BATCH_CACHE_DIR", tmp_path / "graphql_entry_batches")

    batch_response = load_fixture("rcsb_graphql_entry_batch.json")["data"]
    session = FakeSession(post_response={"data": batch_response})

    sc.warm_entry_cache_via_graphql(session, ["1ABC"])

    assert len(session.post_calls) == 1
    entry_path = tmp_path / "entries" / "1ABC.json"
    entity_path = tmp_path / "polymer_entities" / "1ABC_1.json"
    assert entry_path.exists()
    assert entity_path.exists()
    with open(entry_path) as f:
        entry_json = json.load(f)
    assert entry_json["rcsb_accession_info"]["initial_release_date"] == "2019-06-12T00:00:00Z"


def test_warm_entry_cache_via_graphql_falls_back_to_rest_for_ids_graphql_omits(tmp_path, monkeypatch):
    monkeypatch.setattr(sc, "ENTRY_CACHE_DIR", tmp_path / "entries")
    monkeypatch.setattr(sc, "POLYMER_ENTITY_CACHE_DIR", tmp_path / "polymer_entities")
    monkeypatch.setattr(sc, "GRAPHQL_ENTRY_BATCH_CACHE_DIR", tmp_path / "graphql_entry_batches")

    # GraphQL response omits "9ZZZ" entirely (as RCSB does for an unknown/
    # invalid id) -- the REST fallback must be used instead of silently
    # dropping it.
    batch_response = {"entries": []}
    rest_entry = load_fixture("rcsb_entry.json")  # polymer_entity_ids: ["1", "2"]
    rest_entity = load_fixture("rcsb_polymer_entity_normal.json")
    get_responses = {
        sc.ENTRY_URL_TEMPLATE.format(entry_id="9ZZZ"): rest_entry,
        sc.POLYMER_ENTITY_URL_TEMPLATE.format(entry_id="9ZZZ", entity_id="1"): rest_entity,
        sc.POLYMER_ENTITY_URL_TEMPLATE.format(entry_id="9ZZZ", entity_id="2"): rest_entity,
    }
    session = FakeSession(post_response={"data": batch_response}, get_responses=get_responses)

    sc.warm_entry_cache_via_graphql(session, ["9ZZZ"])

    assert (tmp_path / "entries" / "9ZZZ.json").exists()
    assert (tmp_path / "polymer_entities" / "9ZZZ_1.json").exists()
    assert (tmp_path / "polymer_entities" / "9ZZZ_2.json").exists()
    assert sc.ENTRY_URL_TEMPLATE.format(entry_id="9ZZZ") in session.get_calls
