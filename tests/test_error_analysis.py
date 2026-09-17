"""Tests for src.analysis.error_analysis. No network access."""

import numpy as np
import pandas as pd
import pytest

from src.analysis import error_analysis as ea


def _make_features(n_chains=6, residues_per_chain=20, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for c in range(n_chains):
        pdb_id, chain_id = f"P{c:03d}", "A"
        for r in range(residues_per_chain):
            rows.append({
                "pdb_id": pdb_id, "chain_id": chain_id,
                "is_interface_contact": bool(rng.integers(0, 2)),
                "exp_prob": rng.random(), "af_trimmed_prob": rng.random(),
                "plddt": rng.uniform(0, 100), "local_rmsd": rng.uniform(0, 15),
                "rsa": rng.uniform(0, 1), "secondary_structure": rng.choice(["a", "b", "c"]),
                "global_ca_distance": rng.uniform(0, 5), "prob_shift": rng.random(),
            })
    return pd.DataFrame(rows)


# --- band assignment / reporting --------------------------------------------


def test_compute_band_metrics_reports_all_bands_and_flags_undersized():
    features = _make_features(n_chains=8, residues_per_chain=30)
    band_metrics = ea.compute_band_metrics(features)

    assert set(band_metrics["band"]) == set(ea.BAND_ORDER)
    assert set(band_metrics["input"]) == {"exp", "af_trimmed"}
    # every band/input combination present, none silently dropped
    assert len(band_metrics) == len(ea.BAND_ORDER) * 2
    for _, row in band_metrics.iterrows():
        if row["n_residues"] > 0:
            assert row["undersized"] == (row["n_residues"] < ea.config.PLDDT_BAND_MIN_RESIDUES)


def test_compute_band_metrics_empty_band_is_reported_not_skipped():
    # All pLDDT values pinned above 90 -- every other band is empty.
    features = _make_features(n_chains=5, residues_per_chain=10)
    features["plddt"] = 95.0
    band_metrics = ea.compute_band_metrics(features)
    empty_bands = band_metrics[band_metrics["band"] != ">=90"]
    assert (empty_bands["n_residues"] == 0).all()
    assert (empty_bands["undersized"]).all()
    # not silently dropped from the output table
    assert len(band_metrics) == len(ea.BAND_ORDER) * 2


# --- bootstrap resamples chains, not residues -------------------------------


def test_bootstrap_chain_scalar_ci_resamples_at_chain_granularity():
    # Two "chains" with scalar values 1.0 and 100.0. A chain-level bootstrap
    # of the mean can only ever land on {1.0, 50.5, 100.0} (the three
    # possible with-replacement pairs of 2 chosen from 2) -- a residue-level
    # (or otherwise finer-grained) resample would produce a continuum
    # instead. The 2.5th/97.5th percentiles of 1,000 draws should sit
    # exactly at the extremes (each has 25% probability mass).
    values = np.array([1.0, 100.0])
    lo, hi = ea.bootstrap_chain_scalar_ci(values, np.mean, n_resamples=2000, seed=0)
    assert lo == pytest.approx(1.0)
    assert hi == pytest.approx(100.0)


def test_bootstrap_band_metrics_resamples_whole_chains():
    # One "clean" chain (all interface) and one chain with a clear AUC
    # signal; with only 2 chains, chain-level resampling only ever produces
    # {chain0 x2, chain1 x2, chain0+chain1} as the pooled label/prob sets.
    # If it resampled residues instead, undefined (single-class) draws from
    # the always-all-positive chain0 would be far rarer than what chain-level
    # resampling produces (chain0 x2 => still single-class, guaranteed).
    n = 10
    chain0 = pd.DataFrame({
        "pdb_id": ["A"] * n, "chain_id": ["1"] * n,
        "is_interface_contact": [True] * n,
        "exp_prob": np.linspace(0, 1, n), "af_trimmed_prob": np.linspace(0, 1, n),
    })
    rng = np.random.default_rng(0)
    chain1 = pd.DataFrame({
        "pdb_id": ["B"] * n, "chain_id": ["1"] * n,
        "is_interface_contact": rng.integers(0, 2, n).astype(bool),
        "exp_prob": rng.random(n), "af_trimmed_prob": rng.random(n),
    })
    band_df = pd.concat([chain0, chain1], ignore_index=True)

    roc_draws, aupr_draws, n_undefined = ea.bootstrap_band_metrics(band_df, n_resamples=200, seed=0)
    # chain0 alone (or chain0 x2) is always single-class -> at least 1/4 of
    # draws (the "chain0 x2" outcome) must be undefined; since only 2
    # chains exist, this fraction should be roughly 25%, not near 0 (which
    # would indicate residue-level rather than chain-level resampling).
    assert n_undefined > 0
    assert n_undefined / 200 == pytest.approx(0.25, abs=0.15)
    # Defined draws still produced valid metrics in [0, 1].
    assert all(0.0 <= v <= 1.0 for v in roc_draws["exp"])
    assert all(0.0 <= v <= 1.0 for v in aupr_draws["exp"])


# --- regression recovers known coefficients on synthetic data --------------


def _synthetic_regression_df(n_chains=40, residues_per_chain=15, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for c in range(n_chains):
        pdb_id, chain_id = f"P{c:03d}", "A"
        for _ in range(residues_per_chain):
            plddt = rng.uniform(30, 100)
            local_rmsd = rng.uniform(0, 10)
            rsa = rng.uniform(0, 1)
            sse = rng.choice(["a", "b", "c"])
            rows.append({
                "pdb_id": pdb_id, "chain_id": chain_id,
                "plddt": plddt, "local_rmsd": local_rmsd, "rsa": rsa,
                "secondary_structure": sse,
            })
    df = pd.DataFrame(rows)

    plddt_z = (df["plddt"] - df["plddt"].mean()) / df["plddt"].std(ddof=0)
    rmsd_z = (df["local_rmsd"] - df["local_rmsd"].mean()) / df["local_rmsd"].std(ddof=0)
    rsa_z = (df["rsa"] - df["rsa"].mean()) / df["rsa"].std(ddof=0)

    true_coefs = {"const": 0.2, "plddt_z": -0.10, "local_rmsd_z": 0.15, "rsa_z": 0.05}
    rng2 = np.random.default_rng(seed + 1)
    noise = rng2.normal(scale=0.02, size=len(df))
    df["prob_shift"] = (
        true_coefs["const"] + true_coefs["plddt_z"] * plddt_z
        + true_coefs["local_rmsd_z"] * rmsd_z + true_coefs["rsa_z"] * rsa_z + noise
    )
    return df, true_coefs


def test_regression_recovers_known_coefficients():
    df, true_coefs = _synthetic_regression_df()
    result = ea.fit_shift_regression(df)
    result = result.set_index("term")

    assert result.loc["const", "coef"] == pytest.approx(true_coefs["const"], abs=0.03)
    assert result.loc["plddt_z", "coef"] == pytest.approx(true_coefs["plddt_z"], abs=0.03)
    assert result.loc["local_rmsd_z", "coef"] == pytest.approx(true_coefs["local_rmsd_z"], abs=0.03)
    assert result.loc["rsa_z", "coef"] == pytest.approx(true_coefs["rsa_z"], abs=0.03)
    assert result["se_method"].eq("cluster-robust (CRV1), clustered by chain").all()


def test_build_regression_frame_drops_missing_and_logs_count():
    features = _make_features(n_chains=3, residues_per_chain=10)
    features.loc[0, "plddt"] = np.nan
    features.loc[1, "local_rmsd"] = np.nan
    df, info = ea.build_regression_frame(features)
    assert info["n_total"] == 30
    assert info["n_dropped_missing"] == 2
    assert len(df) == 28
    assert df[["plddt", "local_rmsd", "rsa", "prob_shift"]].isna().sum().sum() == 0


# --- VIF computed correctly --------------------------------------------------


def test_vif_flags_collinear_predictors():
    rng = np.random.default_rng(2)
    n_chains, residues_per_chain = 30, 15
    rows = []
    for c in range(n_chains):
        for _ in range(residues_per_chain):
            plddt = rng.uniform(30, 100)
            # local_rmsd almost perfectly determined by plddt -> high VIF for both.
            local_rmsd = -0.1 * plddt + rng.normal(scale=0.05)
            rsa = rng.uniform(0, 1)  # independent
            rows.append({
                "pdb_id": f"P{c:03d}", "chain_id": "A",
                "plddt": plddt, "local_rmsd": local_rmsd, "rsa": rsa,
                "secondary_structure": rng.choice(["a", "b", "c"]),
                "prob_shift": rng.random(),
            })
    df = pd.DataFrame(rows)
    result = ea.fit_shift_regression(df).set_index("term")

    assert result.loc["plddt_z", "vif"] > ea.config.VIF_COLLINEARITY_THRESHOLD
    assert result.loc["local_rmsd_z", "vif"] > ea.config.VIF_COLLINEARITY_THRESHOLD
    assert bool(result.loc["plddt_z", "vif_flagged"]) is True
    assert bool(result.loc["local_rmsd_z", "vif_flagged"]) is True
    # rsa is independent -- should not be flagged.
    assert result.loc["rsa_z", "vif"] < ea.config.VIF_COLLINEARITY_THRESHOLD
    assert bool(result.loc["rsa_z", "vif_flagged"]) is False
