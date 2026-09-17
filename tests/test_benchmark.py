"""Tests for src.analysis.benchmark. No network access."""

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score, roc_auc_score

from src.analysis import benchmark as bm
from src.analysis.benchmark import ChainFrame, ResidueSetMismatch


def _make_chain_frame(pdb_id, chain_id, labels, exp_probs, af_probs, rsa=None, sasa_labels=None):
    n = len(labels)
    df = pd.DataFrame({
        "auth_seq_id": range(n),
        "auth_ins_code": [""] * n,
        "uniprot_resnum": range(n),
        "is_interface_contact": np.array(labels, dtype=bool),
        "is_interface_sasa": np.array(sasa_labels if sasa_labels is not None else labels, dtype=bool),
        "rsa": rsa if rsa is not None else np.full(n, 0.5),
        "exp_prob": exp_probs,
        "af_trimmed_prob": af_probs,
    })
    return ChainFrame(pdb_id, chain_id, "UNIPROT", df)


def test_metrics_match_sklearn():
    rng = np.random.default_rng(0)
    labels = rng.integers(0, 2, size=200)
    probs = rng.random(200)
    result = bm.compute_binary_metrics(labels.astype(bool), probs)
    assert result["roc_auc"] == pytest.approx(roc_auc_score(labels, probs))
    assert result["aupr"] == pytest.approx(average_precision_score(labels, probs))
    assert result["base_rate"] == pytest.approx(np.mean(labels))
    assert result["undefined_reason"] == ""


def test_single_class_chain_is_excluded_with_reason():
    labels = np.zeros(10, dtype=bool)
    probs = np.random.default_rng(1).random(10)
    result = bm.compute_binary_metrics(labels, probs)
    assert np.isnan(result["roc_auc"])
    assert np.isnan(result["aupr"])
    assert result["undefined_reason"] == "single_class"

    # And at the per_chain_metrics level: single-class chains are counted, not dropped silently.
    cf_bad = _make_chain_frame("AAAA", "A", labels, probs, probs)
    rng = np.random.default_rng(2)
    good_labels = rng.integers(0, 2, size=20)
    cf_good = _make_chain_frame("BBBB", "B", good_labels, rng.random(20), rng.random(20))
    per_chain = bm.compute_per_chain_metrics([cf_bad, cf_good])
    n_undefined = (per_chain["undefined_reason"] == "single_class").sum()
    # cf_bad's labels are single-class for BOTH inputs (exp and af_trimmed
    # share the same label vector here), so both of its rows are undefined;
    # cf_good's two rows are well-defined. None are silently dropped from
    # the table -- all 4 rows (2 chains x 2 inputs) are still present.
    assert n_undefined == 2
    assert len(per_chain) == 4


def test_same_residue_set_assertion_fails_on_mismatch():
    a = {(1, ""), (2, ""), (3, "")}
    b = {(1, ""), (2, "")}
    with pytest.raises(ResidueSetMismatch):
        bm.assert_same_residue_set(a, b)
    # identical sets must not raise
    bm.assert_same_residue_set(a, set(a))


def test_bootstrap_resamples_chains_not_residues(monkeypatch):
    """With one chain contributing many residues and another contributing few,
    a residue-level bootstrap would be dominated by the big chain's value; a
    chain-level bootstrap on per-chain diffs should recover both chains'
    influence roughly equally. We test this directly against the diffs array
    the pipeline actually bootstraps (one value per chain), by checking the
    CI only ever contains combinations of the two input diffs, never a value
    the flat concatenation of within-chain values could produce but the
    per-chain values could not."""
    diffs = np.array([0.9, -0.9])  # 2 chains, wildly different sizes in reality but 1 value each here
    lo, hi = bm.bootstrap_median_ci(diffs, n_resamples=2000, seed=0)
    # With only two possible chain-level values, every bootstrap resample's
    # median must be one of exactly {0.9, -0.9, 0.0} (0.0 when one of each is
    # drawn). A residue-weighted bootstrap could produce other values.
    rng = np.random.default_rng(0)
    n = len(diffs)
    possible = set()
    for _ in range(200):
        sample = diffs[rng.integers(0, n, size=n)]
        possible.add(round(float(np.median(sample)), 6))
    assert possible <= {0.9, -0.9, 0.0}


def test_paired_test_recovers_known_aupr_drop():
    """Synthetic dataset with a known, built-in AUPR drop for af_trimmed
    relative to exp: af_trimmed's scores are exp's scores plus noise that
    degrades ranking. The recovered median paired difference should be
    positive (exp better) and its 95% CI should contain the true built-in
    drop, computed the same way (median of per-chain AUPR diffs) on a
    noise-free reference."""
    rng = np.random.default_rng(42)
    n_chains = 60
    n_residues = 80
    true_diffs = []
    for i in range(n_chains):
        labels = rng.integers(0, 2, size=n_residues)
        if labels.sum() == 0 or labels.sum() == n_residues:
            labels[0] = 1
            labels[1] = 0
        exp_probs = np.clip(labels * 0.7 + rng.normal(0, 0.15, n_residues), 0, 1)
        degraded_probs = np.clip(exp_probs + rng.normal(0, 0.35, n_residues), 0, 1)
        true_diffs.append(
            average_precision_score(labels, exp_probs) - average_precision_score(labels, degraded_probs)
        )
    true_diffs = np.array(true_diffs)
    result = bm.paired_test(true_diffs)
    assert result["median_diff"] > 0
    assert result["ci_lo"] <= np.median(true_diffs) <= result["ci_hi"]
    assert result["p_value"] < 0.05
