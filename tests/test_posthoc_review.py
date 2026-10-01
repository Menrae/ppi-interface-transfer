"""Tests for src.analysis.posthoc_review. Synthetic data only, no network."""

import numpy as np
import pandas as pd
import pytest

from src import config
from src.analysis import posthoc_review as ph


# --- template window ----------------------------------------------------------


def test_model_created_date_ignores_isoforms_and_takes_latest():
    body = [
        {"uniprotAccession": "P1", "modelCreatedDate": "2022-06-01T00:00:00Z"},
        {"uniprotAccession": "P1", "modelCreatedDate": "2024-01-01T00:00:00Z"},
        {"uniprotAccession": "P1-2", "modelCreatedDate": "2025-01-01T00:00:00Z"},
    ]
    assert ph.model_created_date(body, "P1") == "2024-01-01"
    assert ph.model_created_date(body, "Q9") is None
    assert ph.model_created_date(None, "P1") is None


def test_classify_template_window():
    assert ph.classify_template_window("2021-05-01", "2022-06-01") == "could_be_template"
    assert ph.classify_template_window("2022-06-01", "2022-06-01") == "after_model"
    assert ph.classify_template_window("2023-01-01", "2022-06-01") == "after_model"
    assert ph.classify_template_window("2023-01-01", None) == "unknown"


def test_template_window_analysis_counts_every_chain():
    chains = pd.DataFrame({
        "pdb_id": ["A", "B", "C", "D"], "chain_id": ["X"] * 4,
        "release_date": ["2021-01-01", "2021-02-01", "2023-01-01", "2023-02-01"],
        "model_created_date": ["2022-06-01"] * 4,
    })
    drops = pd.DataFrame({"pdb_id": ["A", "B", "C"], "chain_id": ["X"] * 3, "drop": [0.1, 0.2, -0.1]})
    summary, per_chain = ph.template_window_analysis(chains, drops)
    s = summary.set_index("template_window")
    assert s.loc["could_be_template", "n_chains"] == 2
    assert s.loc["after_model", "n_chains"] == 2
    assert s.loc["after_model", "n_scoreable"] == 1  # D has no defined drop
    assert len(per_chain) == 4


# --- signed error -------------------------------------------------------------


def test_add_signed_error_sign_convention():
    df = pd.DataFrame({"is_interface_contact": [True, False],
                       "exp_prob": [0.9, 0.1], "af_trimmed_prob": [0.6, 0.05]})
    out = ph.add_signed_error(df)
    # residue 0: AF further from y=1 (0.4 vs 0.1) -> positive
    assert out["error_delta"].iloc[0] == pytest.approx(0.3)
    # residue 1: AF closer to y=0 (0.05 vs 0.1) -> negative
    assert out["error_delta"].iloc[1] == pytest.approx(-0.05)


# --- paired F1 ----------------------------------------------------------------


def test_paired_f1_stats_known_values():
    chains = [{
        "labels": np.array([True, True, False, False]),
        "exp": np.array([0.9, 0.9, 0.1, 0.1]),       # perfect
        "af": np.array([0.9, 0.1, 0.9, 0.1]),        # tp=1, fp=1, fn=1 -> F1 0.5
        "plddt": np.array([95, 40, 40, 95]),
    }]
    stats = ph.paired_f1_stats(chains, lo_cutoff=0, hi_cutoff=90, threshold=0.5)
    assert stats["f1_exp_lo"] == pytest.approx(1.0)
    assert stats["f1_af_lo"] == pytest.approx(0.5)
    # cutoff 90 keeps residues 0 and 3, where AF is also correct
    assert stats["f1_af_hi"] == pytest.approx(1.0)
    assert stats["af_gain"] == pytest.approx(0.5)
    assert stats["gap_change"] == pytest.approx(-0.5)


def test_paired_f1_bootstrap_resamples_chains_jointly():
    rng = np.random.default_rng(0)
    rows = []
    for c in range(6):
        n = 20
        labels = rng.random(n) < 0.4
        rows.append(pd.DataFrame({
            "pdb_id": f"P{c}", "chain_id": "A", "is_interface_contact": labels,
            "exp_prob": np.where(labels, 0.8, 0.2), "af_trimmed_prob": rng.random(n),
            "plddt": rng.uniform(30, 100, n),
        }))
    features = pd.concat(rows, ignore_index=True)
    out = ph.paired_f1_bootstrap(features, n_resamples=50, seed=0).set_index("statistic")
    for k in ("af_gain", "gap_change"):
        assert out.loc[k, "ci_lo"] <= out.loc[k, "ci_hi"]
    # exp input is perfect at every cutoff, so its change is exactly 0 in every draw
    assert out.loc["exp_change", "value"] == pytest.approx(0.0)
    assert out.loc["exp_change", "ci_lo"] == pytest.approx(0.0)


# --- pooled vs median ---------------------------------------------------------


def test_within_chain_rank_removes_chain_offsets():
    df = pd.DataFrame({"pdb_id": ["A", "A", "B", "B"], "chain_id": ["X"] * 4,
                       "p": [0.1, 0.2, 0.8, 0.9]})
    ranks = ph.within_chain_rank(df, "p")
    assert list(ranks) == [0.5, 1.0, 0.5, 1.0]


def test_top_drop_chains_picks_largest():
    drops = pd.DataFrame({"pdb_id": list("ABCDEFGHIJ"), "chain_id": ["X"] * 10,
                          "drop": [0.0, 0.5, 0.1, 0.9, 0.2, 0.3, 0.05, 0.4, 0.6, 0.7]})
    assert ph.top_drop_chains(drops, 0.2) == {("D", "X"), ("J", "X")}


# --- flag audit ---------------------------------------------------------------


def test_flag_mapping_audit_groups():
    chains = pd.DataFrame({
        "geometrically_flagged": [True, True, False, False, False],
        "mapping_method": ["sifts", "fallback", "sifts", "sifts", "fallback"],
        "mapping_identity": [1.0, 0.95, 1.0, 0.92, 1.0],
    })
    out = ph.flag_mapping_audit(chains).set_index("group")
    assert out.loc["flagged", "n_chains"] == 2
    assert out.loc["flagged", "n_fallback"] == 1
    assert out.loc["unflagged", "n_identity_below_1"] == 1


def test_extreme_local_rmsd_audit_counts_per_chain():
    features = pd.DataFrame({"pdb_id": ["A", "A", "A", "B"], "chain_id": ["X"] * 4,
                             "local_rmsd": [12.0, 15.0, 1.0, 0.5]})
    info = pd.DataFrame({"pdb_id": ["A", "B"], "chain_id": ["X", "X"], "geometrically_flagged": [True, False]})
    out = ph.extreme_local_rmsd_audit(features, info, threshold=10.0)
    assert len(out) == 1
    row = out.iloc[0]
    assert row["n_extreme"] == 2 and row["n_residues"] == 3 and row["max_local_rmsd"] == 15.0
