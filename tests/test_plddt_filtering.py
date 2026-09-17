"""Tests for src.analysis.plddt_filtering. No network access."""

import numpy as np
import pandas as pd
import pytest

from src.analysis import plddt_filtering as pf


# --- compute_prf -------------------------------------------------------------


def test_compute_prf_known_values():
    labels = np.array([True, True, False, False, True])
    preds = np.array([True, False, False, True, True])
    # tp=2 (idx0,4), fp=1 (idx3), fn=1 (idx1)
    result = pf.compute_prf(labels, preds)
    assert result["precision"] == pytest.approx(2 / 3)
    assert result["recall"] == pytest.approx(2 / 3)
    assert result["f1"] == pytest.approx(2 / 3)


def test_compute_prf_undefined_when_no_predicted_positives():
    labels = np.array([True, False, True])
    preds = np.array([False, False, False])
    result = pf.compute_prf(labels, preds)
    assert np.isnan(result["precision"])  # 0/0
    assert result["recall"] == 0.0  # tp=0, fn=2 -> defined, just 0
    assert np.isnan(result["f1"])


def test_compute_prf_undefined_when_no_true_positives():
    labels = np.array([False, False, False])
    preds = np.array([True, False, True])
    result = pf.compute_prf(labels, preds)
    assert result["precision"] == 0.0  # tp=0, fp=2 -> defined, just 0
    assert np.isnan(result["recall"])  # 0/0
    assert np.isnan(result["f1"])


# --- compute_filtering_curve --------------------------------------------------


def _make_features(n_chains=10, residues_per_chain=20, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for c in range(n_chains):
        for _ in range(residues_per_chain):
            rows.append({
                "pdb_id": f"P{c:03d}", "chain_id": "A",
                "plddt": rng.uniform(0, 100),
                "is_interface_contact": bool(rng.integers(0, 2)),
                "exp_prob": rng.random(), "af_trimmed_prob": rng.random(),
            })
    return pd.DataFrame(rows)


def test_cutoff_zero_keeps_every_residue():
    features = _make_features()
    curve = pf.compute_filtering_curve(features, cutoffs=(0,), n_resamples=50)
    assert (curve["n_residues"] == len(features)).all()
    assert np.allclose(curve["frac_residues_kept"], 1.0)
    assert np.allclose(curve["frac_positives_kept"], 1.0)


def test_cutoff_above_max_plddt_keeps_zero_and_is_reported():
    features = _make_features()
    features["plddt"] = 50.0  # fixed, well below an unreachable cutoff
    curve = pf.compute_filtering_curve(features, cutoffs=(200,), n_resamples=50)
    assert len(curve) == 2  # both inputs still get a row, not silently dropped
    assert (curve["n_residues"] == 0).all()
    assert (curve["frac_residues_kept"] == 0).all()
    assert curve["precision"].isna().all()
    assert curve["recall"].isna().all()
    assert curve["f1"].isna().all()


def test_higher_cutoff_never_increases_kept_fraction():
    features = _make_features(n_chains=15, residues_per_chain=30)
    curve = pf.compute_filtering_curve(features, cutoffs=(0, 50, 70, 90), n_resamples=50)
    for input_name in pf.INPUTS:
        sub = curve[curve["input"] == input_name].sort_values("cutoff")
        fracs = sub["frac_residues_kept"].to_numpy()
        assert all(fracs[i] >= fracs[i + 1] for i in range(len(fracs) - 1))


def test_perfect_predictor_scores_perfectly_regardless_of_cutoff():
    features = _make_features(n_chains=8, residues_per_chain=25)
    # af_trimmed_prob perfectly encodes the true label at a 0.5 threshold.
    features["af_trimmed_prob"] = np.where(features["is_interface_contact"], 0.9, 0.1)
    curve = pf.compute_filtering_curve(features, cutoffs=(0, 50), n_resamples=50)
    af = curve[curve["input"] == "af_trimmed"]
    assert (af["precision"] == 1.0).all()
    assert (af["recall"] == 1.0).all()
    assert (af["f1"] == 1.0).all()


# --- bootstrap resamples chains, not residues --------------------------------


def test_bootstrap_prf_ci_resamples_at_chain_granularity():
    # Two chains: one all-interface (guarantees recall well-defined but
    # precision/recall trivially 1 when isolated), one all-non-interface
    # with some false positives. With only 2 chains, resampling can only
    # ever produce {chain0 x2, chain1 x2, chain0+chain1} as pooled sets.
    chain0_labels = np.array([True] * 10)
    chain0_probs = np.full(10, 0.9)  # all predicted positive, perfect precision/recall on its own
    chain1_labels = np.array([False] * 10)
    chain1_probs = np.array([0.9] * 3 + [0.1] * 7)  # 3 false positives, precision=0 on its own

    per_chain = [(chain0_labels, chain0_probs), (chain1_labels, chain1_probs)]
    draws, n_undefined = pf.bootstrap_prf_ci(per_chain, threshold=0.5, n_resamples=500, seed=0)

    # chain1 x2 alone has zero true positives -> recall undefined for that outcome;
    # with only 2 chains this should happen a substantial (~25%) fraction of draws,
    # not a rare/near-zero fraction (which would indicate residue-level resampling).
    assert n_undefined["recall"] / 500 == pytest.approx(0.25, abs=0.15)
    assert all(0.0 <= v <= 1.0 for v in draws["precision"])
    assert all(0.0 <= v <= 1.0 for v in draws["recall"])


def test_bootstrap_prf_ci_empty_per_chain_returns_no_draws():
    draws, n_undefined = pf.bootstrap_prf_ci([], threshold=0.5, n_resamples=100, seed=0)
    assert draws == {"precision": [], "recall": [], "f1": []}
    assert n_undefined == {"precision": 0, "recall": 0, "f1": 0}
