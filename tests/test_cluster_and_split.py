import gzip
import tempfile
from pathlib import Path

import pytest

from src import config
from src.data import cluster_and_split as cas

FIXTURES = Path(__file__).parent / "fixtures"


# --- chain id formatting -----------------------------------------------


def test_chain_id_str_uppercases_pdb_id_only():
    assert cas.chain_id_str("1abc", "A") == "1ABC_A"
    assert cas.chain_id_str("1ABC", "a") == "1ABC_a"


# --- deterministic representative selection -----------------------------


def test_choose_representative_prefers_best_resolution():
    index = {
        "1ABC_A": {"resolution": "2.5", "seq_length": "100"},
        "2XYZ_B": {"resolution": "1.2", "seq_length": "80"},
    }
    assert cas.choose_representative(["1ABC_A", "2XYZ_B"], index) == "2XYZ_B"


def test_choose_representative_breaks_resolution_tie_with_longest_chain():
    index = {
        "1ABC_A": {"resolution": "2.0", "seq_length": "100"},
        "2XYZ_B": {"resolution": "2.0", "seq_length": "250"},
    }
    assert cas.choose_representative(["1ABC_A", "2XYZ_B"], index) == "2XYZ_B"


def test_choose_representative_breaks_full_tie_with_id_sort_order():
    index = {
        "2XYZ_B": {"resolution": "2.0", "seq_length": "100"},
        "1ABC_A": {"resolution": "2.0", "seq_length": "100"},
    }
    assert cas.choose_representative(["2XYZ_B", "1ABC_A"], index) == "1ABC_A"


def test_choose_representative_is_deterministic_regardless_of_input_order():
    index = {
        "1ABC_A": {"resolution": "1.8", "seq_length": "200"},
        "2XYZ_B": {"resolution": "1.8", "seq_length": "200"},
        "3DEF_C": {"resolution": "2.4", "seq_length": "300"},
    }
    members = ["1ABC_A", "2XYZ_B", "3DEF_C"]
    assert cas.choose_representative(members, index) == cas.choose_representative(
        list(reversed(members)), index
    )


# --- MMseqs2 output parsing (fixture files, no mmseqs invocation) -------


def test_parse_cluster_tsv_groups_by_mmseqs_representative():
    groups = cas.parse_cluster_tsv(FIXTURES / "mmseqs_cluster_sample.tsv")
    assert groups == {
        "1ABC_A": ["1ABC_A", "2XYZ_B"],
        "3DEF_C": ["3DEF_C"],
    }


def test_parse_search_tsv_reads_hits_and_casts_fident():
    hits = cas.parse_search_tsv(FIXTURES / "mmseqs_search_sample.tsv")
    assert set(hits.keys()) == {"1ABC_A", "3DEF_C"}
    assert len(hits["1ABC_A"]) == 2
    targets = {h["target"]: h["fident"] for h in hits["1ABC_A"]}
    assert targets["4UF6_C"] == pytest.approx(0.95)
    assert targets["4GBS_A"] == pytest.approx(0.40)
    assert isinstance(hits["1ABC_A"][0]["fident"], float)


# --- leakage flag computation --------------------------------------------


def test_compute_leakage_flags_picks_best_hit_per_split_and_flags_overlap():
    split_membership = {
        "train": {"4UF6_C"},
        "test": {"4GBS_A"},
        "validation": set(),
    }
    search_hits = cas.parse_search_tsv(FIXTURES / "mmseqs_search_sample.tsv")

    flags = cas.compute_leakage_flags(["1ABC_A", "3DEF_C"], search_hits, split_membership)

    rep = flags["1ABC_A"]
    assert rep["max_identity_train"] == pytest.approx(0.95)
    assert rep["best_hit_train_id"] == "4UF6_C"
    assert rep["max_identity_test"] == pytest.approx(0.40)
    assert rep["best_hit_test_id"] == "4GBS_A"
    assert rep["max_identity_validation"] == 0.0
    assert rep["max_identity_any"] == pytest.approx(0.95)
    assert rep["pesto_homolog_overlap"] is True  # 0.95 >= SEQUENCE_IDENTITY_CUTOFF
    assert rep["pesto_exact_train_overlap"] is False  # "1ABC_A" not literally in train set

    other = flags["3DEF_C"]
    # 3DEF_C's only hit (3U22_A, 0.25) isn't in any split set -> no overlap at all.
    assert other["max_identity_any"] == 0.0
    assert other["pesto_homolog_overlap"] is False


def test_compute_leakage_flags_exact_train_overlap_is_literal_id_membership():
    split_membership = {"train": {"1ABC_A"}, "test": set(), "validation": set()}
    flags = cas.compute_leakage_flags(["1ABC_A"], {}, split_membership)
    assert flags["1ABC_A"]["pesto_exact_train_overlap"] is True
    assert flags["1ABC_A"]["pesto_homolog_overlap"] is False  # no search hits at all


def test_leakage_mode_eligible_homolog_excludes_overlap():
    overlap_row = {"pesto_homolog_overlap": True, "pesto_exact_train_overlap": False}
    clean_row = {"pesto_homolog_overlap": False, "pesto_exact_train_overlap": False}
    assert cas.leakage_mode_eligible("homolog", overlap_row) is False
    assert cas.leakage_mode_eligible("homolog", clean_row) is True


def test_leakage_mode_eligible_exact_train_only_looks_at_exact_flag():
    row = {"pesto_homolog_overlap": True, "pesto_exact_train_overlap": False}
    # Homolog-overlapping but not an exact training-set match -> still eligible
    # under the narrower exact_train mode.
    assert cas.leakage_mode_eligible("exact_train", row) is True
    row["pesto_exact_train_overlap"] = True
    assert cas.leakage_mode_eligible("exact_train", row) is False


def test_leakage_mode_eligible_none_is_always_true():
    row = {"pesto_homolog_overlap": True, "pesto_exact_train_overlap": True}
    assert cas.leakage_mode_eligible("none", row) is True


def test_leakage_mode_eligible_rejects_unknown_mode():
    with pytest.raises(ValueError):
        cas.leakage_mode_eligible("bogus", {"pesto_homolog_overlap": False, "pesto_exact_train_overlap": False})


def test_leakage_eligibility_fields_match_configured_modes():
    assert cas.LEAKAGE_ELIGIBILITY_FIELDS == [f"eligible_{m}" for m in config.LEAKAGE_FILTER_MODES]
    for field in cas.LEAKAGE_ELIGIBILITY_FIELDS:
        assert field in cas.DEDUP_FIELDS


# --- threshold sweep -------------------------------------------------------


def _synthetic_flags(n_below_030, n_between_030_050, n_above_095, n_exact_only):
    """Build a synthetic {rep_id: flags} table with a known survivor profile."""
    flags = {}
    i = 0
    for _ in range(n_below_030):
        flags[f"r{i}"] = {"max_identity_any": 0.10, "pesto_exact_train_overlap": False}
        i += 1
    for _ in range(n_between_030_050):
        flags[f"r{i}"] = {"max_identity_any": 0.35, "pesto_exact_train_overlap": False}
        i += 1
    for _ in range(n_above_095):
        flags[f"r{i}"] = {"max_identity_any": 0.99, "pesto_exact_train_overlap": False}
        i += 1
    for _ in range(n_exact_only):
        flags[f"r{i}"] = {"max_identity_any": 0.0, "pesto_exact_train_overlap": True}
        i += 1
    return flags


def test_sweep_thresholds_counts_survivors_per_threshold():
    flags = _synthetic_flags(n_below_030=120, n_between_030_050=10, n_above_095=5, n_exact_only=3)
    rows = cas.sweep_thresholds(flags)
    by_mode = {r["mode"]: r for r in rows}

    # identity<0.30: only the 120 below-0.30 rows + the 3 exact-only rows (max_identity_any=0.0) survive.
    assert by_mode["identity<0.30"]["n_survive"] == 123
    # identity<0.50: below-0.30 + between-0.30-0.50 are all < 0.50 -> excluded (0.35 < 0.50 survives)
    assert by_mode["identity<0.50"]["n_survive"] == 120 + 10 + 3
    # identity<0.95: everything except the 5 at 0.99 survives
    assert by_mode["identity<0.95"]["n_survive"] == 120 + 10 + 3
    # identity<0.30 threshold and MIN_TEST_CHAINS_TARGET/FLOOR from config
    assert by_mode["identity<0.30"]["meets_target"] == (123 >= config.MIN_TEST_CHAINS_TARGET)
    assert by_mode["identity<0.30"]["meets_floor"] == (123 >= config.MIN_TEST_CHAINS_FLOOR)


def test_sweep_thresholds_exact_id_mode_counts_non_exact_overlap_only():
    flags = _synthetic_flags(n_below_030=40, n_between_030_050=0, n_above_095=0, n_exact_only=5)
    rows = cas.sweep_thresholds(flags)
    by_mode = {r["mode"]: r for r in rows}
    # exact_id_train_only survives everything except the 5 flagged pesto_exact_train_overlap=True.
    assert by_mode["exact_id_train_only"]["n_survive"] == 40
    assert by_mode["exact_id_train_only"]["meets_floor"] is False  # 40 < MIN_TEST_CHAINS_FLOOR (50)


def test_sweep_thresholds_covers_every_configured_threshold_plus_exact_id():
    flags = _synthetic_flags(1, 0, 0, 0)
    rows = cas.sweep_thresholds(flags)
    modes = {r["mode"] for r in rows}
    for threshold in config.LEAKAGE_SWEEP_IDENTITY_THRESHOLDS:
        assert f"identity<{threshold:.2f}" in modes
    assert "exact_id_train_only" in modes


# --- PeSTo split membership loading --------------------------------------


def test_load_pesto_split_membership_preserves_chain_id_case(tmp_path, monkeypatch):
    # Chain IDs are case-sensitive in the PDB (e.g. "C" and "c" can be
    # distinct chains in the same entry) -- only the PDB ID half should be
    # normalized to uppercase, not the whole "PDBID_CHAINID" line.
    train_path = tmp_path / "train.txt"
    train_path.write_text("4uf6_C\n4uf6_c\n1abc_A\n")
    test_path = tmp_path / "test.txt"
    test_path.write_text("\n")
    val_path = tmp_path / "val.txt"
    val_path.write_text("\n")

    monkeypatch.setitem(cas.PESTO_SPLIT_FILES, "train", train_path)
    monkeypatch.setitem(cas.PESTO_SPLIT_FILES, "test", test_path)
    monkeypatch.setitem(cas.PESTO_SPLIT_FILES, "validation", val_path)

    membership = cas.load_pesto_split_membership()

    assert membership["train"] == {"4UF6_C", "4UF6_c", "1ABC_A"}


# --- pdb_seqres bulk lookup (synthetic gzip fixture, no network) --------


def test_lookup_pdb_seqres_sequences_finds_wanted_and_reports_missing():
    content = (
        ">100d_A mol:na length:10  some description\n"
        "CCGGCGCCGG\n"
        ">4uf6_C mol:protein length:5  a protein\n"
        "MKVLA\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        gz_path = Path(tmp) / "pdb_seqres_sample.txt.gz"
        with gzip.open(gz_path, "wt") as f:
            f.write(content)

        found, missing = cas.lookup_pdb_seqres_sequences(gz_path, {"4UF6_C", "9ZZZ_Z"})

    assert found == {"4UF6_C": "MKVLA"}
    assert missing == {"9ZZZ_Z"}


# --- real MMseqs2 integration (skipped if the binary isn't installed) ---


MMSEQS_AVAILABLE = cas.MMSEQS_BIN.exists()


@pytest.mark.skipif(not MMSEQS_AVAILABLE, reason="MMseqs2 binary not installed at external/mmseqs/")
def test_mmseqs_cluster_and_search_roundtrip_on_tiny_synthetic_sequences(tmp_path):
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(
        ">A_1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEE\n"
        ">B_1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEQ\n"
        ">C_1\nGSHMASMTGGQQMGRGSEFEEQLNQWLQERKQMLEELNTIKQQVQELTAKLSEKSSELAAAHAELVKLQAENQQLKAQ\n"
    )

    cluster_prefix = tmp_path / "out" / "cluster"
    cluster_tmp = tmp_path / "cluster_tmp"
    cluster_tsv = cas.run_mmseqs_cluster(query_fasta, cluster_prefix, cluster_tmp)
    groups = cas.parse_cluster_tsv(cluster_tsv)
    # A_1 and B_1 differ by a single residue (>99% identity) -> same cluster;
    # C_1 is unrelated -> its own cluster.
    all_members = {m for members in groups.values() for m in members}
    assert all_members == {"A_1", "B_1", "C_1"}
    cluster_of = {m: rep for rep, members in groups.items() for m in members}
    assert cluster_of["A_1"] == cluster_of["B_1"]
    assert cluster_of["C_1"] != cluster_of["A_1"]

    target_fasta = tmp_path / "target.fasta"
    target_fasta.write_text(
        ">TARGET_X\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEE\n"
    )
    search_out = tmp_path / "search.tsv"
    search_tmp = tmp_path / "search_tmp"
    result_path = cas.run_mmseqs_search(query_fasta, target_fasta, search_out, search_tmp)
    hits = cas.parse_search_tsv(result_path)
    assert "A_1" in hits
    assert hits["A_1"][0]["target"] == "TARGET_X"
    assert hits["A_1"][0]["fident"] == pytest.approx(1.0)
