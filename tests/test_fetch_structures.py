import csv
import json
from pathlib import Path

import pytest
import requests

from src import config
from src.data import fetch_structures as fs

FIXTURES = Path(__file__).parent / "fixtures"


def load_fixture(name):
    with open(FIXTURES / name) as f:
        return json.load(f)


class FakeResponse:
    def __init__(self, status_code=200, json_data=None, text_data="", content=b""):
        self.status_code = status_code
        self._json_data = json_data
        self.text = text_data
        self.content = content

    def json(self):
        return self._json_data


class FakeSession:
    """Minimal requests.Session stand-in: exact-URL-keyed responses, fails
    loudly (AssertionError) on any call not explicitly wired up, so a "no
    network happened" assertion is just "nothing needed to be wired up"."""

    def __init__(self, responses: dict | None = None):
        self._responses = responses or {}
        self.calls: list[str] = []

    def get(self, url, timeout=None, **kwargs):
        self.calls.append(url)
        if url not in self._responses:
            raise AssertionError(f"unexpected GET to {url} (test wired no response)")
        outcome = self._responses[url]
        if isinstance(outcome, Exception):
            raise outcome
        return outcome


@pytest.fixture(autouse=True)
def _isolate_cache_dirs(tmp_path, monkeypatch):
    """Redirect every cache directory to tmp_path and remove the polite
    delay, so tests never touch the real data/raw/ tree and run fast."""
    monkeypatch.setattr(fs, "PDB_UPDATED_DIR", tmp_path / "pdb")
    monkeypatch.setattr(fs, "ALPHAFOLD_DIR", tmp_path / "alphafold")
    monkeypatch.setattr(fs, "ALPHAFOLD_API_CACHE_DIR", tmp_path / "alphafold" / "api")
    monkeypatch.setattr(fs, "UNIPROT_CACHE_DIR", tmp_path / "uniprot")
    monkeypatch.setattr(config, "PHASE3_REQUEST_DELAY_SECONDS", 0)
    yield tmp_path


# --- pure URL builders ---------------------------------------------------


def test_pdbe_updated_cif_url_lowercases_pdb_id():
    assert fs.pdbe_updated_cif_url("7YF2") == (
        "https://www.ebi.ac.uk/pdbe/entry-files/download/7yf2_updated.cif"
    )


def test_alphafold_api_url_uses_accession_verbatim():
    assert fs.alphafold_api_url("Q9H1A3") == "https://alphafold.ebi.ac.uk/api/prediction/Q9H1A3"


def test_uniprot_fasta_url_uses_accession_verbatim():
    assert fs.uniprot_fasta_url("Q9H1A3") == "https://rest.uniprot.org/uniprotkb/Q9H1A3.fasta"


# --- pure response classification, on saved fixtures ---------------------


def test_classify_pdb_http_status_success():
    assert fs.classify_pdb_http_status(200) is None


def test_classify_pdb_http_status_failure_is_distinct_reason():
    assert fs.classify_pdb_http_status(404) == "pdb_download_error:404"
    assert fs.classify_pdb_http_status(500) == "pdb_download_error:500"


def test_classify_alphafold_api_response_success_matches_exact_accession():
    body = load_fixture("alphafold_api_success.json")
    failure, entry = fs.classify_alphafold_api_response(200, body, "Q9H1A3")
    assert failure is None
    assert entry["uniprotAccession"] == "Q9H1A3"


def test_classify_alphafold_api_response_accession_absent_on_404():
    failure, entry = fs.classify_alphafold_api_response(404, None, "Q9H1A3")
    assert failure == fs.FAILURE_ACCESSION_ABSENT
    assert entry is None


def test_classify_alphafold_api_response_accession_absent_on_empty_body():
    failure, entry = fs.classify_alphafold_api_response(200, [], "Q9H1A3")
    assert failure == fs.FAILURE_ACCESSION_ABSENT


def test_classify_alphafold_api_response_isoforms_only_is_accession_absent():
    # AFDB returned isoform entries (e.g. "Q9UPN9-2") but none matching the
    # queried base accession "Q9UPN9" exactly -- not a match.
    body = load_fixture("alphafold_api_isoforms_only.json")
    failure, entry = fs.classify_alphafold_api_response(200, body, "Q9UPN9")
    assert failure == fs.FAILURE_ACCESSION_ABSENT
    assert entry is None


def test_classify_alphafold_api_response_multiple_exact_matches_is_fragmented():
    body = load_fixture("alphafold_api_fragmented.json")
    failure, entry = fs.classify_alphafold_api_response(200, body, "Q9Y4C0")
    assert failure == fs.FAILURE_FRAGMENTED
    assert entry is None


def test_parse_fasta_sequence_strips_header_and_joins_lines():
    text = ">sp|Q9H1A3|METL9_HUMAN Protein\nMRLLAGWLCL\nSLASVWLARR\n"
    assert fs.parse_fasta_sequence(text) == "MRLLAGWLCLSLASVWLARR"


def test_sequences_match_true_for_identical_nonempty():
    assert fs.sequences_match("MKV", "MKV") is True


def test_sequences_match_false_for_different_or_empty():
    assert fs.sequences_match("MKV", "MKA") is False
    assert fs.sequences_match("", "MKV") is False


# --- priority ordering -----------------------------------------------------


def test_leakage_priority_rank_orders_homolog_then_exact_train_then_rest():
    primary = {"eligible_homolog": True, "eligible_exact_train": True, "cluster_id": "9ZZZ_A"}
    sensitivity_a = {"eligible_homolog": False, "eligible_exact_train": True, "cluster_id": "1AAA_A"}
    rest = {"eligible_homolog": False, "eligible_exact_train": False, "cluster_id": "1AAA_A"}
    assert fs.leakage_priority_rank(primary)[0] == 0
    assert fs.leakage_priority_rank(sensitivity_a)[0] == 1
    assert fs.leakage_priority_rank(rest)[0] == 2


def test_sort_by_priority_puts_primary_set_first_deterministically():
    rows = [
        {"eligible_homolog": False, "eligible_exact_train": False, "cluster_id": "3CCC_A"},
        {"eligible_homolog": True, "eligible_exact_train": True, "cluster_id": "2BBB_A"},
        {"eligible_homolog": True, "eligible_exact_train": True, "cluster_id": "1AAA_A"},
    ]
    ordered = fs.sort_by_priority(rows)
    assert [r["cluster_id"] for r in ordered] == ["1AAA_A", "2BBB_A", "3CCC_A"]


def test_unique_ids_in_order_preserves_first_seen_order():
    rows = [{"pdb_id": "7YF2"}, {"pdb_id": "6YR7"}, {"pdb_id": "7YF2"}]
    assert fs.unique_ids_in_order(rows, "pdb_id") == ["7YF2", "6YR7"]


# --- fetch report row ------------------------------------------------------


def _rep(**overrides):
    rep = {
        "cluster_id": "7YF2_A",
        "pdb_id": "7YF2",
        "chain_id": "A",
        "uniprot_acc": "Q9H1A3",
        "release_date": "2022-01-01",
        "eligible_homolog": True,
        "eligible_exact_train": True,
        "eligible_none": True,
    }
    rep.update(overrides)
    return rep


def _pdb_result(status="success", failure_reason=""):
    return {"path": "/tmp/x.cif.gz" if status == "success" else "", "status": status, "failure_reason": failure_reason, "bytes": 0, "network_call": False}


def _af_result(status="success", failure_reason=""):
    return {
        "cif_path": "/tmp/y.cif.gz" if status == "success" else "",
        "status": status,
        "failure_reason": failure_reason,
        "model_version": 6,
        "uniprot_start": 1,
        "uniprot_end": 10,
        "bytes": 0,
        "network_call": False,
    }


def test_build_fetch_report_row_schema_matches_declared_fields():
    row = fs.build_fetch_report_row(_rep(), _pdb_result(), _af_result())
    assert set(row.keys()) == set(fs.FETCH_REPORT_FIELDS)


def test_build_fetch_report_row_overall_status_complete_partial_failed():
    both_ok = fs.build_fetch_report_row(_rep(), _pdb_result("success"), _af_result("success"))
    assert both_ok["overall_status"] == "complete"

    pdb_only = fs.build_fetch_report_row(
        _rep(), _pdb_result("success"), _af_result("failed", fs.FAILURE_ACCESSION_ABSENT)
    )
    assert pdb_only["overall_status"] == "partial"

    neither = fs.build_fetch_report_row(
        _rep(), _pdb_result("failed", fs.FAILURE_PDB_DOWNLOAD_ERROR), _af_result("failed", fs.FAILURE_NETWORK_ERROR)
    )
    assert neither["overall_status"] == "failed"


# --- load_representatives / candidates_dedup.csv parsing -------------------


def test_load_representatives_parses_booleans_and_renames_uniprot_column(tmp_path):
    path = tmp_path / "candidates_dedup.csv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "cluster_id", "pdb_id", "chain_id", "release_date", "uniprot_ids",
                "eligible_homolog", "eligible_exact_train", "eligible_none",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "cluster_id": "7YF2_A", "pdb_id": "7YF2", "chain_id": "A",
                "release_date": "2022-01-01", "uniprot_ids": "Q9H1A3",
                "eligible_homolog": True, "eligible_exact_train": True, "eligible_none": True,
            }
        )
    rows = fs.load_representatives(path)
    assert len(rows) == 1
    assert rows[0]["uniprot_acc"] == "Q9H1A3"
    assert rows[0]["eligible_homolog"] is True
    assert isinstance(rows[0]["eligible_homolog"], bool)


# --- compute_attrition -------------------------------------------------


def test_compute_attrition_counts_per_leakage_mode():
    report_rows = [
        {"eligible_homolog": True, "eligible_exact_train": True, "eligible_none": True,
         "pdb_status": "success", "alphafold_status": "success"},
        {"eligible_homolog": False, "eligible_exact_train": True, "eligible_none": True,
         "pdb_status": "success", "alphafold_status": "failed"},
        {"eligible_homolog": False, "eligible_exact_train": False, "eligible_none": True,
         "pdb_status": "failed", "alphafold_status": "failed"},
    ]
    attrition = fs.compute_attrition(report_rows)
    by_mode = {r["leakage_filter_mode"]: r for r in attrition}

    assert by_mode["homolog"]["n_eligible_representatives"] == 1
    assert by_mode["homolog"]["n_both_success"] == 1

    assert by_mode["exact_train"]["n_eligible_representatives"] == 2
    assert by_mode["exact_train"]["n_both_success"] == 1
    assert by_mode["exact_train"]["n_pdb_success"] == 2
    assert by_mode["exact_train"]["n_alphafold_success"] == 1

    assert by_mode["none"]["n_eligible_representatives"] == 3
    assert by_mode["none"]["n_both_success"] == 1


# --- fetch_pdb_structure: cache hit makes zero network calls ---------------


def test_fetch_pdb_structure_cache_hit_makes_no_network_call(tmp_path):
    cache_path = fs.PDB_UPDATED_DIR / "7YF2_updated.cif.gz"
    cache_path.parent.mkdir(parents=True)
    cache_path.write_bytes(b"fake-gzip-bytes")

    session = FakeSession()  # no responses wired -> any .get() call fails the test
    result = fs.fetch_pdb_structure(session, "7YF2")

    assert result["status"] == "success"
    assert result["path"] == str(cache_path)
    assert session.calls == []


def test_fetch_pdb_structure_success_downloads_and_gzips():
    url = fs.pdbe_updated_cif_url("7YF2")
    session = FakeSession({url: FakeResponse(200, content=b"data_7YF2\n#\n")})
    result = fs.fetch_pdb_structure(session, "7YF2")

    assert result["status"] == "success"
    assert Path(result["path"]).exists()
    assert session.calls == [url]


def test_fetch_pdb_structure_http_error_is_pdb_download_error():
    url = fs.pdbe_updated_cif_url("9ZZZ")
    session = FakeSession({url: FakeResponse(404)})
    result = fs.fetch_pdb_structure(session, "9ZZZ")

    assert result["status"] == "failed"
    assert result["failure_reason"] == "pdb_download_error:404"


def test_fetch_pdb_structure_network_error_is_network_error():
    url = fs.pdbe_updated_cif_url("9ZZZ")
    session = FakeSession({url: requests.ConnectionError("boom")})
    result = fs.fetch_pdb_structure(session, "9ZZZ")

    assert result["status"] == "failed"
    assert result["failure_reason"] == fs.FAILURE_NETWORK_ERROR


# --- fetch_alphafold_model: full success/failure paths, no network on cache hit --


def test_fetch_alphafold_model_fully_cached_makes_no_network_call(tmp_path):
    acc = "Q9H1A3"
    api_path = fs.ALPHAFOLD_API_CACHE_DIR / f"{acc}.json"
    api_path.parent.mkdir(parents=True)
    api_path.write_text(json.dumps({"status_code": 200, "body": [
        {"uniprotAccession": acc, "uniprotSequence": "MKV", "uniprotStart": 1, "uniprotEnd": 3,
         "latestVersion": 6, "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-Q9H1A3-F1-model_v6.cif"}
    ]}))
    uniprot_path = fs.UNIPROT_CACHE_DIR / f"{acc}.fasta"
    uniprot_path.parent.mkdir(parents=True)
    uniprot_path.write_text(">sp|Q9H1A3|X\nMKV\n")
    cif_path = fs.ALPHAFOLD_DIR / f"{acc}.cif.gz"
    cif_path.write_bytes(b"fake-gzip-bytes")

    session = FakeSession()
    result = fs.fetch_alphafold_model(session, acc)

    assert result["status"] == "success"
    assert session.calls == []


def test_fetch_alphafold_model_success_downloads_cif_and_records_metadata():
    acc = "Q9H1A3"
    api_url = fs.alphafold_api_url(acc)
    uniprot_url = fs.uniprot_fasta_url(acc)
    cif_url = "https://alphafold.ebi.ac.uk/files/AF-Q9H1A3-F1-model_v6.cif"
    entry = {
        "uniprotAccession": acc, "uniprotSequence": "MKV", "uniprotStart": 1, "uniprotEnd": 3,
        "latestVersion": 6, "cifUrl": cif_url,
    }
    session = FakeSession(
        {
            api_url: FakeResponse(200, json_data=[entry]),
            uniprot_url: FakeResponse(200, text_data=">sp|Q9H1A3|X\nMKV\n"),
            cif_url: FakeResponse(200, content=b"data_AF\n#\n"),
        }
    )
    result = fs.fetch_alphafold_model(session, acc)

    assert result["status"] == "success"
    assert result["model_version"] == 6
    assert result["uniprot_start"] == 1
    assert result["uniprot_end"] == 3
    assert Path(result["cif_path"]).exists()


def test_fetch_alphafold_model_accession_absent():
    acc = "Q0000000"
    api_url = fs.alphafold_api_url(acc)
    session = FakeSession({api_url: FakeResponse(404)})
    result = fs.fetch_alphafold_model(session, acc)

    assert result["status"] == "failed"
    assert result["failure_reason"] == fs.FAILURE_ACCESSION_ABSENT


def test_fetch_alphafold_model_fragmented():
    acc = "Q9Y4C0"
    api_url = fs.alphafold_api_url(acc)
    body = [
        {"uniprotAccession": acc, "entryId": "AF-Q9Y4C0-F1"},
        {"uniprotAccession": acc, "entryId": "AF-Q9Y4C0-F2"},
    ]
    session = FakeSession({api_url: FakeResponse(200, json_data=body)})
    result = fs.fetch_alphafold_model(session, acc)

    assert result["status"] == "failed"
    assert result["failure_reason"] == fs.FAILURE_FRAGMENTED


def test_fetch_alphafold_model_sequence_mismatch_never_downloads_cif():
    acc = "Q9H1A3"
    api_url = fs.alphafold_api_url(acc)
    uniprot_url = fs.uniprot_fasta_url(acc)
    cif_url = "https://alphafold.ebi.ac.uk/files/AF-Q9H1A3-F1-model_v6.cif"
    entry = {
        "uniprotAccession": acc, "uniprotSequence": "MKV", "uniprotStart": 1, "uniprotEnd": 3,
        "latestVersion": 6, "cifUrl": cif_url,
    }
    session = FakeSession(
        {
            api_url: FakeResponse(200, json_data=[entry]),
            # Current UniProt sequence differs from the AF model's.
            uniprot_url: FakeResponse(200, text_data=">sp|Q9H1A3|X\nMKA\n"),
        }
    )
    result = fs.fetch_alphafold_model(session, acc)

    assert result["status"] == "failed"
    assert result["failure_reason"] == fs.FAILURE_SEQUENCE_MISMATCH
    assert cif_url not in session.calls


def test_fetch_alphafold_model_network_error_on_api_call():
    acc = "Q9H1A3"
    api_url = fs.alphafold_api_url(acc)
    session = FakeSession({api_url: requests.ConnectionError("boom")})
    result = fs.fetch_alphafold_model(session, acc)

    assert result["status"] == "failed"
    assert result["failure_reason"] == fs.FAILURE_NETWORK_ERROR
