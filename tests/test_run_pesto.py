import csv
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

from src import config
from src.models import run_pesto as rp


def _write_labels_report(path: Path, rows: list[dict]) -> None:
    fields = [
        "cluster_id", "pdb_id", "chain_id", "uniprot_acc", "status",
        "eligible_homolog", "eligible_exact_train", "eligible_none",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_fetch_report(path: Path, rows: list[dict]) -> None:
    fields = ["pdb_id", "chain_id", "uniprot_acc", "pdb_cif_path", "alphafold_cif_path"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_residue_map(path: Path, keys: list[tuple[int, str]], uniprot_start: int = 100) -> None:
    rows = [
        {"auth_seq_id": k[0], "auth_ins_code": k[1], "uniprot_resnum": uniprot_start + i}
        for i, k in enumerate(keys)
    ]
    pd.DataFrame(rows).to_parquet(path, index=False)


@pytest.fixture()
def _isolated_paths(tmp_path, monkeypatch):
    labels_path = tmp_path / "labels_report.csv"
    fetch_path = tmp_path / "fetch_report.csv"
    residue_maps_dir = tmp_path / "residue_mappings"
    inputs_dir = tmp_path / "pesto_inputs"
    predictions_dir = tmp_path / "predictions"
    residue_maps_dir.mkdir()

    monkeypatch.setattr(rp, "LABELS_REPORT_PATH", labels_path)
    monkeypatch.setattr(rp, "FETCH_REPORT_PATH", fetch_path)
    monkeypatch.setattr(rp, "RESIDUE_MAPS_DIR", residue_maps_dir)
    monkeypatch.setattr(rp, "PESTO_INPUTS_DIR", inputs_dir)
    monkeypatch.setattr(rp, "PREDICTIONS_DIR", predictions_dir)

    return {
        "labels_path": labels_path, "fetch_path": fetch_path,
        "residue_maps_dir": residue_maps_dir, "inputs_dir": inputs_dir, "predictions_dir": predictions_dir,
    }


# --- --subset / --inputs select the right jobs -----------------------------


def test_load_representatives_filters_by_subset(_isolated_paths):
    _write_labels_report(
        _isolated_paths["labels_path"],
        [
            {"cluster_id": "A", "pdb_id": "AAAA", "chain_id": "A", "uniprot_acc": "P1", "status": "labeled",
             "eligible_homolog": True, "eligible_exact_train": True, "eligible_none": True},
            {"cluster_id": "B", "pdb_id": "BBBB", "chain_id": "A", "uniprot_acc": "P2", "status": "labeled",
             "eligible_homolog": False, "eligible_exact_train": True, "eligible_none": True},
            {"cluster_id": "C", "pdb_id": "CCCC", "chain_id": "A", "uniprot_acc": "P3", "status": "labeled",
             "eligible_homolog": False, "eligible_exact_train": False, "eligible_none": True},
            {"cluster_id": "D", "pdb_id": "DDDD", "chain_id": "A", "uniprot_acc": "P4", "status": "excluded",
             "eligible_homolog": True, "eligible_exact_train": True, "eligible_none": True},
        ],
    )
    _write_fetch_report(
        _isolated_paths["fetch_path"],
        [{"pdb_id": p, "chain_id": "A", "uniprot_acc": u, "pdb_cif_path": "x", "alphafold_cif_path": "y"}
         for p, u in [("AAAA", "P1"), ("BBBB", "P2"), ("CCCC", "P3"), ("DDDD", "P4")]],
    )

    primary = rp.load_representatives("primary")
    assert {r["pdb_id"] for r in primary} == {"AAAA"}

    exact_train = rp.load_representatives("exact_train")
    assert {r["pdb_id"] for r in exact_train} == {"AAAA", "BBBB"}

    none_mode = rp.load_representatives("none")
    assert {r["pdb_id"] for r in none_mode} == {"AAAA", "BBBB", "CCCC"}  # DDDD excluded by Phase 5, not this filter


def _fake_prepare_experimental_input(pdb_cif_path, chain_id, phase4_keys, out_pdb, out_keys):
    rp._write_keys([[k[0], k[1]] for k in phase4_keys], ["auth_seq_id", "auth_ins_code"], out_keys)
    return 100, {k: 0 for k in phase4_keys}


def _fake_prepare_alphafold_input(af_cif_path, uniprot_range, out_pdb, out_keys):
    rp._write_keys([[100], [101]], ["uniprot_resnum"], out_keys)
    return 100


def test_build_jobs_respects_requested_inputs(monkeypatch, _isolated_paths):
    _write_residue_map(_isolated_paths["residue_maps_dir"] / "AAAA_A.parquet", [(1, ""), (2, "")])
    monkeypatch.setattr(rp, "prepare_experimental_input", _fake_prepare_experimental_input)
    monkeypatch.setattr(rp, "prepare_alphafold_input", _fake_prepare_alphafold_input)

    rep = {"cluster_id": "A", "pdb_id": "AAAA", "chain_id": "A", "uniprot_acc": "P1",
           "pdb_cif_path": "x", "alphafold_cif_path": "y"}

    jobs = rp.build_jobs([rep], ("exp",))
    assert [j.variant for j in jobs] == ["exp"]

    jobs = rp.build_jobs([rep], ("exp", "af_trimmed"))
    assert sorted(j.variant for j in jobs) == ["af_trimmed", "exp"]


# --- oversized inputs skipped; the other input still runs ------------------


def test_oversized_input_skipped_other_input_still_runs(monkeypatch, _isolated_paths):
    monkeypatch.setattr(config, "PESTO_MAX_ATOMS", 500)
    _write_residue_map(_isolated_paths["residue_maps_dir"] / "AAAA_A.parquet", [(1, ""), (2, "")])

    # exp is oversized, af_trimmed is not.
    def oversized_exp(pdb_cif_path, chain_id, phase4_keys, out_pdb, out_keys):
        rp._write_keys([[k[0], k[1]] for k in phase4_keys], ["auth_seq_id", "auth_ins_code"], out_keys)
        return 999, {k: 0 for k in phase4_keys}

    monkeypatch.setattr(rp, "prepare_experimental_input", oversized_exp)
    monkeypatch.setattr(rp, "prepare_alphafold_input", _fake_prepare_alphafold_input)

    rep = {"cluster_id": "A", "pdb_id": "AAAA", "chain_id": "A", "uniprot_acc": "P1",
           "pdb_cif_path": "x", "alphafold_cif_path": "y"}
    jobs = rp.build_jobs([rep], ("exp", "af_trimmed"))
    by_variant = {j.variant: j for j in jobs}

    assert by_variant["exp"].status == "skipped_too_large"
    assert by_variant["exp"].reason == rp.REASON_EXCEEDS_MEMORY
    assert by_variant["af_trimmed"].status == "pending"


# --- resumability: existing valid outputs are skipped -----------------------


def test_existing_valid_output_is_skipped_without_reprocessing(monkeypatch, _isolated_paths):
    out_path = rp.output_path_for("exp", "AAAA", "A")
    out_path.parent.mkdir(parents=True)
    pd.DataFrame({"auth_seq_id": [1], "auth_ins_code": [""], "pesto_interface_prob": [0.5]}).to_parquet(out_path)

    _write_residue_map(_isolated_paths["residue_maps_dir"] / "AAAA_A.parquet", [(1, "")])

    prep_called = MagicMock()
    monkeypatch.setattr(rp, "prepare_experimental_input", prep_called)

    rep = {"cluster_id": "A", "pdb_id": "AAAA", "chain_id": "A", "uniprot_acc": "P1",
           "pdb_cif_path": "x", "alphafold_cif_path": "y"}
    jobs = rp.build_jobs([rep], ("exp",))

    assert jobs[0].status == "skipped_resumed"
    prep_called.assert_not_called()


def test_is_valid_output_rejects_missing_or_empty_file(tmp_path):
    missing = tmp_path / "missing.parquet"
    assert rp.is_valid_output(missing) is False

    empty = tmp_path / "empty.parquet"
    pd.DataFrame({"auth_seq_id": [], "pesto_interface_prob": []}).to_parquet(empty)
    assert rp.is_valid_output(empty) is False

    valid = tmp_path / "valid.parquet"
    pd.DataFrame({"auth_seq_id": [1], "pesto_interface_prob": [0.1]}).to_parquet(valid)
    assert rp.is_valid_output(valid) is True


# --- a crashed/OOM-killed worker is logged and the run continues -----------


def test_crashed_worker_logged_and_other_jobs_still_run(monkeypatch, _isolated_paths):
    def fake_run(cmd, capture_output, text, timeout):
        pdb_path = Path(cmd[4])
        output_path = Path(cmd[6])
        if "CRASH" in pdb_path.name:
            return MagicMock(returncode=-9, stderr="")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"auth_seq_id": [1], "auth_ins_code": [""], "pesto_interface_prob": [0.5]}).to_parquet(output_path)
        return MagicMock(returncode=0, stderr="")

    monkeypatch.setattr(rp.subprocess, "run", fake_run)

    ok_job = rp.Job(
        cluster_id="A", pdb_id="AAAA", chain_id="A", uniprot_acc="P1", variant="exp",
        pdb_path=Path("/tmp/OK.pdb"), keys_path=Path("/tmp/OK.keys.csv"),
        output_path=_isolated_paths["predictions_dir"] / "exp" / "AAAA_A.parquet",
    )
    crash_job = rp.Job(
        cluster_id="B", pdb_id="BBBB", chain_id="A", uniprot_acc="P2", variant="exp",
        pdb_path=Path("/tmp/CRASH.pdb"), keys_path=Path("/tmp/CRASH.keys.csv"),
        output_path=_isolated_paths["predictions_dir"] / "exp" / "BBBB_A.parquet",
    )

    rp.run_tiered([ok_job, crash_job])

    assert crash_job.status == "failed"
    assert crash_job.reason.startswith(rp.REASON_WORKER_CRASHED)
    assert ok_job.status == "success"  # the crash didn't stop the other job
