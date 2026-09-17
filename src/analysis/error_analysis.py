"""Phase 8 (lean): error analysis on the primary set (exp vs. af_trimmed).

Implements the pre-registered analysis plan in PLAN.md ("Phase 8 analysis
plan (pre-registered 2026-09-17...)"), committed to before this module was
run against real predictions. Do not add or change endpoints here without
updating that pre-registration first. All results are associations, not
causal claims (pLDDT and local RMSD are expected to correlate -- see VIF
in the regression output).

Reuses Phase 7's chain-level join (`src.analysis.benchmark.load_all_chain_frames`,
`build_chain_frame`) rather than re-deriving it, then augments each chain
with Phase 8's own structural features (`src.analysis.structural_metrics`).

Runnable as: .venv/bin/python -m src.analysis.error_analysis
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.outliers_influence import variance_inflation_factor

from src import config
from src.analysis import benchmark as bm
from src.analysis import structural_metrics as sme

logger = logging.getLogger(__name__)

# --- Paths -----------------------------------------------------------------

PDB_UPDATED_DIR = config.RAW_DATA_DIR / "pdb"
ALPHAFOLD_DIR = config.RAW_DATA_DIR / "alphafold"
LABELS_REPORT_PATH = config.INTERIM_DATA_DIR / "labels_report.csv"

ERROR_ANALYSIS_DIR = config.RESULTS_DIR / "error_analysis"
FIGURES_DIR = config.RESULTS_DIR / "figures"
PER_RESIDUE_FEATURES_PATH = ERROR_ANALYSIS_DIR / "per_residue_features.parquet"
BAND_METRICS_PATH = ERROR_ANALYSIS_DIR / "band_metrics.csv"
REGRESSION_PATH = ERROR_ANALYSIS_DIR / "regression.csv"
LENGTH_ANALYSIS_PATH = ERROR_ANALYSIS_DIR / "length_analysis.csv"
FLAGGED_CHAINS_SUMMARY_PATH = ERROR_ANALYSIS_DIR / "flagged_chains_summary.csv"
STRUCTURAL_EXCLUSIONS_PATH = ERROR_ANALYSIS_DIR / "structural_feature_exclusions.csv"

BAND_ORDER = sme.plddt_band_labels()
REGRESSION_PREDICTORS = ("plddt", "local_rmsd", "rsa")
SSE_REFERENCE_CATEGORY = "c"  # coil


# --- Building per-chain structural features --------------------------------


def compute_chain_structural_features(cf: bm.ChainFrame) -> tuple[pd.DataFrame | None, str]:
    """Augments a Phase-7 ChainFrame's residue table with pLDDT, local RMSD,
    global Calpha distance, and secondary structure. Returns (df, "") or
    (None, reason) if the chain can't be processed at all."""
    exp_path = PDB_UPDATED_DIR / f"{cf.pdb_id}_updated.cif.gz"
    af_path = ALPHAFOLD_DIR / f"{cf.uniprot_acc}.cif.gz"
    if not (exp_path.exists() and af_path.exists()):
        return None, "missing_structure_file"

    df = cf.df.copy()
    auth_keys = list(zip(df["auth_seq_id"].astype(int), df["auth_ins_code"]))
    uniprot_keys = list(df["uniprot_resnum"].astype(int))

    exp_coords_map = sme.read_ca_coords_by_auth(exp_path, cf.chain_id, auth_keys)
    af_coords_map = sme.read_ca_coords_by_uniprot(af_path, uniprot_keys)

    keep = np.array([a in exp_coords_map and u in af_coords_map for a, u in zip(auth_keys, uniprot_keys)])
    n_dropped = int((~keep).sum())
    if n_dropped:
        logger.warning(
            "%s_%s: dropping %d/%d residues missing a resolvable Calpha in exp and/or af_trimmed",
            cf.pdb_id, cf.chain_id, n_dropped, len(df),
        )
    if keep.sum() < config.LOCAL_RMSD_MIN_NEIGHBORS:
        return None, "too_few_ca_resolved_residues"

    df = df.loc[keep].reset_index(drop=True)
    auth_keys = [k for k, ok in zip(auth_keys, keep) if ok]
    uniprot_keys = [k for k, ok in zip(uniprot_keys, keep) if ok]

    exp_coords = np.array([exp_coords_map[k] for k in auth_keys])
    af_coords = np.array([af_coords_map[k] for k in uniprot_keys])

    local_rmsd, n_neighbors, reasons = sme.local_rmsd_per_residue(
        exp_coords, af_coords, config.LOCAL_RMSD_RADIUS_ANGSTROM, config.LOCAL_RMSD_MIN_NEIGHBORS
    )
    global_dist = sme.global_ca_distance_per_residue(exp_coords, af_coords)
    plddt_map = sme.read_plddt_by_uniprot(af_path, uniprot_keys)
    sse_map, sse_method = sme.read_secondary_structure_by_auth(exp_path, cf.chain_id)

    df["plddt"] = [plddt_map.get(k, np.nan) for k in uniprot_keys]
    df["local_rmsd"] = local_rmsd
    df["local_rmsd_n_neighbors"] = n_neighbors
    df["local_rmsd_reason"] = reasons
    df["global_ca_distance"] = global_dist
    df["secondary_structure"] = [sse_map.get(k, "") for k in auth_keys]
    df["secondary_structure_method"] = sse_method
    df["pdb_id"] = cf.pdb_id
    df["chain_id"] = cf.chain_id
    df["uniprot_acc"] = cf.uniprot_acc
    df["prob_shift"] = (df["af_trimmed_prob"] - df["exp_prob"]).abs()

    n_missing_plddt = int(df["plddt"].isna().sum())
    if n_missing_plddt:
        logger.warning(
            "%s_%s: %d residues missing plddt (uniprot key not resolved in the AlphaFold model)",
            cf.pdb_id, cf.chain_id, n_missing_plddt,
        )
    n_unassigned_sse = int((df["secondary_structure"] == "").sum())
    if n_unassigned_sse:
        logger.info(
            "%s_%s: %d/%d residues have no P-SEA secondary-structure assignment",
            cf.pdb_id, cf.chain_id, n_unassigned_sse, len(df),
        )

    return df, ""


def load_all_structural_features() -> tuple[pd.DataFrame, list[dict]]:
    chain_frames, join_exclusions = bm.load_all_chain_frames()
    logger.info(
        "Loaded %d chain frames from Phase 7's join (%d already excluded there)",
        len(chain_frames), len(join_exclusions),
    )
    exclusions = [{**exc, "stage": "phase7_join"} for exc in join_exclusions]

    feature_dfs = []
    for cf in chain_frames:
        df, reason = compute_chain_structural_features(cf)
        if df is None:
            exclusions.append({"pdb_id": cf.pdb_id, "chain_id": cf.chain_id, "reason": reason, "stage": "phase8_structural"})
        else:
            feature_dfs.append(df)

    if not feature_dfs:
        raise RuntimeError("No chains produced structural features -- nothing to analyze")
    combined = pd.concat(feature_dfs, ignore_index=True)
    return combined, exclusions


def write_exclusions(exclusions: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pdb_id", "chain_id", "reason", "stage"])
        writer.writeheader()
        writer.writerows(exclusions)


# --- Chain-resampled bootstrap helpers --------------------------------------


def bootstrap_chain_scalar_ci(values: np.ndarray, statistic, n_resamples: int, seed: int) -> tuple[float, float]:
    """95% percentile bootstrap CI for `statistic(values)`, resampling
    values (one scalar per chain) with replacement -- cheap, same
    methodology as Phase 7's bootstrap_median_ci but for an arbitrary
    statistic."""
    rng = np.random.default_rng(seed)
    n = len(values)
    stats = np.empty(n_resamples)
    for i in range(n_resamples):
        sample = values[rng.integers(0, n, size=n)]
        stats[i] = statistic(sample)
    return float(np.percentile(stats, 2.5)), float(np.percentile(stats, 97.5))


def bootstrap_paired_stat_ci(x: np.ndarray, y: np.ndarray, statistic, n_resamples: int, seed: int) -> tuple[float, float]:
    """95% percentile bootstrap CI for `statistic(x, y)`, resampling
    (x[i], y[i]) pairs jointly with replacement (one pair per chain)."""
    rng = np.random.default_rng(seed)
    n = len(x)
    stats = np.empty(n_resamples)
    for i in range(n_resamples):
        idx = rng.integers(0, n, size=n)
        stats[i] = statistic(x[idx], y[idx])
    return float(np.percentile(stats, 2.5)), float(np.percentile(stats, 97.5))


def bootstrap_band_metrics(band_df: pd.DataFrame, n_resamples: int, seed: int) -> tuple[dict, dict, int]:
    """Chain-resampled bootstrap of pooled ROC-AUC/AUPR for both inputs at
    once (same resampled chain sets for exp and af_trimmed, so the two
    inputs' CIs are computed on directly comparable draws). Returns
    (roc_auc_draws, aupr_draws, n_undefined_draws), each draws dict keyed
    "exp"/"af_trimmed".

    n_resamples is deliberately smaller than Phase 7's chain-scalar
    bootstrap here (see config.PHASE8_BAND_BOOTSTRAP_N_RESAMPLES) -- each
    draw recomputes a pooled O(n log n) metric over real residue counts,
    not a median of ~500 scalars.
    """
    per_chain = [
        (g["is_interface_contact"].to_numpy(dtype=bool), g["exp_prob"].to_numpy(), g["af_trimmed_prob"].to_numpy())
        for _, g in band_df.groupby(["pdb_id", "chain_id"], sort=False)
    ]
    n_chains = len(per_chain)
    rng = np.random.default_rng(seed)
    roc_auc_draws = {"exp": [], "af_trimmed": []}
    aupr_draws = {"exp": [], "af_trimmed": []}
    n_undefined = 0
    for _ in range(n_resamples):
        idx = rng.integers(0, n_chains, size=n_chains)
        labels = np.concatenate([per_chain[i][0] for i in idx])
        if labels.min() == labels.max():
            n_undefined += 1
            continue
        for j, input_name in enumerate(("exp", "af_trimmed"), start=1):
            probs = np.concatenate([per_chain[i][j] for i in idx])
            roc_auc_draws[input_name].append(bm.roc_auc_score(labels, probs))
            aupr_draws[input_name].append(bm.average_precision_score(labels, probs))
    return roc_auc_draws, aupr_draws, n_undefined


# --- Primary analysis: pLDDT-band metrics -----------------------------------


def compute_band_metrics(features: pd.DataFrame) -> pd.DataFrame:
    features = features.copy()
    features["plddt_band"] = sme.assign_plddt_band(features["plddt"].to_numpy())

    rows = []
    for band in BAND_ORDER:
        band_df = features[features["plddt_band"] == band]
        n_residues = len(band_df)
        n_chains = band_df[["pdb_id", "chain_id"]].drop_duplicates().shape[0]
        undersized = n_residues < config.PLDDT_BAND_MIN_RESIDUES

        if n_residues == 0:
            logger.warning("pLDDT band %s: empty -- reported as zero rows, not silently skipped", band)
            for input_name in ("exp", "af_trimmed"):
                rows.append({"band": band, "input": input_name, "n_residues": 0, "n_chains": 0, "undersized": True})
            continue
        if undersized:
            logger.warning(
                "pLDDT band %s: %d residues (< PLDDT_BAND_MIN_RESIDUES=%d) -- reported, flagged undersized",
                band, n_residues, config.PLDDT_BAND_MIN_RESIDUES,
            )

        roc_auc_draws, aupr_draws, n_undefined = bootstrap_band_metrics(
            band_df, config.PHASE8_BAND_BOOTSTRAP_N_RESAMPLES, config.PHASE8_BOOTSTRAP_SEED
        )
        labels_arr = band_df["is_interface_contact"].to_numpy(dtype=bool)
        for input_name in ("exp", "af_trimmed"):
            point = bm.compute_binary_metrics(labels_arr, band_df[f"{input_name}_prob"].to_numpy())
            roc_vals, aupr_vals = roc_auc_draws[input_name], aupr_draws[input_name]
            rows.append({
                "band": band, "input": input_name, "n_residues": n_residues, "n_chains": n_chains,
                "undersized": undersized, "base_rate": point["base_rate"],
                "roc_auc": point["roc_auc"],
                "roc_auc_ci_lo": float(np.percentile(roc_vals, 2.5)) if roc_vals else float("nan"),
                "roc_auc_ci_hi": float(np.percentile(roc_vals, 97.5)) if roc_vals else float("nan"),
                "aupr": point["aupr"],
                "aupr_ci_lo": float(np.percentile(aupr_vals, 2.5)) if aupr_vals else float("nan"),
                "aupr_ci_hi": float(np.percentile(aupr_vals, 97.5)) if aupr_vals else float("nan"),
                "n_bootstrap_draws_undefined": n_undefined,
            })
    return pd.DataFrame(rows)


# --- Secondary analysis: shift regression -----------------------------------


def build_regression_frame(features: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    cols = ["pdb_id", "chain_id", "plddt", "local_rmsd", "rsa", "secondary_structure", "prob_shift"]
    df = features[cols].copy()
    n_total = len(df)
    missing_mask = df[["plddt", "local_rmsd", "rsa", "prob_shift"]].isna().any(axis=1)
    n_dropped = int(missing_mask.sum())
    df = df.loc[~missing_mask].reset_index(drop=True)
    return df, {"n_total": n_total, "n_dropped_missing": n_dropped}


def fit_shift_regression(df: pd.DataFrame) -> pd.DataFrame:
    """OLS of |p_af_trimmed - p_exp| on standardized pLDDT/local RMSD/RSA
    plus secondary-structure dummies (coil reference), with cluster-robust
    (CRV1) standard errors clustered by chain. Reports coefficients, SEs,
    p-values, and VIF per predictor (flagged if > config.VIF_COLLINEARITY_THRESHOLD).
    """
    y = df["prob_shift"].astype(float).to_numpy()

    predictors = pd.DataFrame(index=df.index)
    for col in REGRESSION_PREDICTORS:
        std = df[col].std(ddof=0)
        predictors[f"{col}_z"] = (df[col] - df[col].mean()) / std if std > 0 else 0.0

    sse_filled = df["secondary_structure"].where(df["secondary_structure"] != "", "unassigned")
    sse_dummies = pd.get_dummies(sse_filled, prefix="sse", drop_first=False)
    ref_col = f"sse_{SSE_REFERENCE_CATEGORY}"
    dummy_cols = [c for c in sse_dummies.columns if c != ref_col]
    predictors = pd.concat([predictors, sse_dummies[dummy_cols].astype(float)], axis=1)

    design = sm.add_constant(predictors.astype(float))
    groups = (df["pdb_id"] + "_" + df["chain_id"]).to_numpy()

    model = sm.OLS(y, design.to_numpy()).fit(cov_type="cluster", cov_kwds={"groups": groups})

    design_values = design.to_numpy()
    vifs = {}
    for i, col in enumerate(design.columns):
        if col == "const":
            continue
        vifs[col] = float(variance_inflation_factor(design_values, i))

    rows = []
    for i, term in enumerate(design.columns):
        rows.append({
            "term": term,
            "coef": float(model.params[i]),
            "cluster_robust_se": float(model.bse[i]),
            "t": float(model.tvalues[i]),
            "p_value": float(model.pvalues[i]),
            "vif": vifs.get(term, float("nan")),
            "vif_flagged": vifs.get(term, 0.0) > config.VIF_COLLINEARITY_THRESHOLD,
            "se_method": "cluster-robust (CRV1), clustered by chain",
            "n": len(df),
            "n_clusters": int(pd.unique(groups).shape[0]),
            "reference_category": f"secondary_structure={SSE_REFERENCE_CATEGORY}",
        })
    return pd.DataFrame(rows)


# --- Descriptive: chain-length vs. AUPR drop --------------------------------


def compute_length_analysis(per_chain_metrics_path: Path = bm.PER_CHAIN_METRICS_PATH) -> pd.DataFrame:
    per_chain = pd.read_csv(per_chain_metrics_path)
    wide = per_chain[per_chain["label"] == "is_interface_contact"].pivot_table(
        index=["pdb_id", "chain_id"], columns="input", values="aupr"
    ).dropna()
    wide["aupr_drop"] = wide["exp"] - wide["af_trimmed"]
    wide = wide.reset_index()

    labels_report = pd.read_csv(LABELS_REPORT_PATH)[["pdb_id", "chain_id", "n_observed"]]
    merged = wide.merge(labels_report, on=["pdb_id", "chain_id"], how="left")
    n_missing = int(merged["n_observed"].isna().sum())
    if n_missing:
        logger.warning("length analysis: %d/%d chains missing n_observed after merge, dropped", n_missing, len(merged))
    merged = merged.dropna(subset=["n_observed"])

    x = merged["n_observed"].to_numpy(dtype=float)
    y = merged["aupr_drop"].to_numpy(dtype=float)

    rho, p = spearmanr(x, y)
    rho_lo, rho_hi = bootstrap_paired_stat_ci(
        x, y, lambda a, b: spearmanr(a, b)[0], config.BOOTSTRAP_N_RESAMPLES, config.PHASE8_BOOTSTRAP_SEED
    )

    ols = sm.OLS(y, sm.add_constant(x)).fit()
    slope, slope_lo, slope_hi = ols.params[1], ols.conf_int()[1][0], ols.conf_int()[1][1]

    rows = [
        {"row_type": "correlation", "stratum": "spearman_length_vs_aupr_drop", "n": len(merged),
         "value": rho, "ci_lo": rho_lo, "ci_hi": rho_hi, "p_value": p},
        {"row_type": "regression_slope", "stratum": "aupr_drop_on_length_ols_slope", "n": len(merged),
         "value": slope, "ci_lo": slope_lo, "ci_hi": slope_hi, "p_value": ols.pvalues[1]},
        {"row_type": "summary", "stratum": "unweighted_mean_drop", "n": len(merged), "value": float(np.mean(y))},
        {"row_type": "summary", "stratum": "unweighted_median_drop", "n": len(merged), "value": float(np.median(y))},
        {"row_type": "summary", "stratum": "length_weighted_mean_drop", "n": len(merged),
         "value": float(np.average(y, weights=x))},
    ]

    quantile_labels = [f"Q{i + 1}" for i in range(config.LENGTH_ANALYSIS_N_QUANTILES)]
    merged = merged.copy()
    merged["length_quantile"] = pd.qcut(merged["n_observed"], config.LENGTH_ANALYSIS_N_QUANTILES, labels=quantile_labels)
    for q, g in merged.groupby("length_quantile", observed=True):
        rows.append({
            "row_type": "length_quantile", "stratum": str(q), "n": len(g),
            "value": float(g["aupr_drop"].mean()), "mean_length": float(g["n_observed"].mean()),
            "median_drop": float(g["aupr_drop"].median()),
        })

    return pd.DataFrame(rows)


# --- Descriptive: geometrically flagged vs. unflagged chains ----------------


def compute_flagged_chains_summary(
    features: pd.DataFrame, per_chain_metrics_path: Path = bm.PER_CHAIN_METRICS_PATH
) -> pd.DataFrame:
    meta = bm.load_chain_metadata()[["pdb_id", "chain_id", "geometrically_flagged"]].copy()
    # bm.load_chain_metadata()'s own left-join-then-fillna leaves this column
    # as object dtype (True/False Python bools, no NaN) rather than proper
    # bool -- harmless there (only ever .map()'d), but `~` on an object
    # array of bools does a bitwise Python int invert (~True == -2), not a
    # boolean negation. Cast explicitly before using `~` below.
    meta["geometrically_flagged"] = meta["geometrically_flagged"].astype(bool)

    per_chain = pd.read_csv(per_chain_metrics_path)
    wide = per_chain[per_chain["label"] == "is_interface_contact"].pivot_table(
        index=["pdb_id", "chain_id"], columns="input", values="aupr"
    ).dropna()
    wide["aupr_drop"] = wide["exp"] - wide["af_trimmed"]
    wide = wide.reset_index()

    labels_report = pd.read_csv(LABELS_REPORT_PATH)[["pdb_id", "chain_id", "n_observed", "interface_fraction_distance"]]

    per_chain_features = features.groupby(["pdb_id", "chain_id"]).agg(
        mean_plddt=("plddt", "mean"),
        mean_local_rmsd=("local_rmsd", "mean"),
        mean_global_ca_distance=("global_ca_distance", "mean"),
    ).reset_index()

    merged = wide.merge(meta, on=["pdb_id", "chain_id"], how="left")
    merged = merged.merge(labels_report, on=["pdb_id", "chain_id"], how="left")
    merged = merged.merge(per_chain_features, on=["pdb_id", "chain_id"], how="left")
    merged["geometrically_flagged"] = merged["geometrically_flagged"].fillna(False).astype(bool)

    metric_columns = {
        "aupr_drop": "aupr_drop",
        "chain_length": "n_observed",
        "interface_fraction": "interface_fraction_distance",
        "mean_plddt": "mean_plddt",
        "mean_local_rmsd": "mean_local_rmsd",
        "mean_global_ca_distance": "mean_global_ca_distance",
    }
    flagged = merged[merged["geometrically_flagged"]]
    unflagged = merged[~merged["geometrically_flagged"]]

    rows = []
    for label, col in metric_columns.items():
        f_vals = flagged[col].dropna().to_numpy()
        u_vals = unflagged[col].dropna().to_numpy()
        if len(f_vals) == 0 or len(u_vals) == 0:
            logger.warning("flagged-chain comparison: metric %s has an empty group -- no Mann-Whitney test", label)
            stat, p = float("nan"), float("nan")
        else:
            stat, p = mannwhitneyu(f_vals, u_vals, alternative="two-sided")
        rows.append({
            "metric": label, "n_flagged": len(f_vals), "n_unflagged": len(u_vals),
            "mean_flagged": float(np.mean(f_vals)) if len(f_vals) else float("nan"),
            "median_flagged": float(np.median(f_vals)) if len(f_vals) else float("nan"),
            "mean_unflagged": float(np.mean(u_vals)) if len(u_vals) else float("nan"),
            "median_unflagged": float(np.median(u_vals)) if len(u_vals) else float("nan"),
            "mannwhitney_stat": stat, "mannwhitney_p": p,
        })
    return pd.DataFrame(rows)


# --- Figures -----------------------------------------------------------------


def plot_band_metrics(band_metrics: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), facecolor="#fcfcfb")
    width = 0.35
    x = np.arange(len(BAND_ORDER))
    for ax, metric in zip(axes, ("aupr", "roc_auc")):
        ax.set_facecolor("#fcfcfb")
        for offset, input_name, color in ((-width / 2, "exp", bm.COLOR_EXP), (width / 2, "af_trimmed", bm.COLOR_AF_TRIMMED)):
            sub = band_metrics[band_metrics["input"] == input_name].set_index("band").reindex(BAND_ORDER)
            values = sub[metric].to_numpy()
            lo = values - sub[f"{metric}_ci_lo"].to_numpy()
            hi = sub[f"{metric}_ci_hi"].to_numpy() - values
            ax.bar(x + offset, values, width=width, color=color, label=input_name,
                   yerr=[lo, hi], capsize=3, error_kw={"linewidth": 1})
        ax.set_xticks(x)
        ax.set_xticklabels(BAND_ORDER)
        ax.set_xlabel("AlphaFold pLDDT band")
        ax.set_ylabel(metric.upper() if metric == "aupr" else "ROC-AUC")
        ax.set_title(f"Pooled {metric.upper() if metric == 'aupr' else 'ROC-AUC'} by pLDDT band")
        bm._style_axes(ax)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_shift_vs_local_rmsd(features: pd.DataFrame, path: Path, n_bins: int = 10) -> None:
    df = features.dropna(subset=["local_rmsd", "prob_shift"]).copy()
    df["rmsd_bin"] = pd.qcut(df["local_rmsd"], n_bins, duplicates="drop")

    rows = []
    for interval, g in df.groupby("rmsd_bin", observed=True):
        chain_means = g.groupby(["pdb_id", "chain_id"])["prob_shift"].mean().to_numpy()
        lo, hi = bootstrap_chain_scalar_ci(chain_means, np.mean, config.BOOTSTRAP_N_RESAMPLES, config.PHASE8_BOOTSTRAP_SEED)
        rows.append({"mid": interval.mid, "mean_shift": g["prob_shift"].mean(), "ci_lo": lo, "ci_hi": hi})
    binned = pd.DataFrame(rows).sort_values("mid")

    fig, ax = plt.subplots(figsize=(6, 4.5), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.plot(binned["mid"], binned["mean_shift"], color=bm.COLOR_EXP, linewidth=2, marker="o")
    ax.fill_between(binned["mid"], binned["ci_lo"], binned["ci_hi"], color=bm.COLOR_EXP, alpha=0.2, linewidth=0)
    ax.set_xlabel("Local Calpha RMSD (Å), binned")
    ax.set_ylabel("Mean |p_af_trimmed - p_exp|")
    ax.set_title("Prediction shift vs. local structural divergence")
    bm._style_axes(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_aupr_drop_vs_length(length_source_path: Path, path: Path) -> None:
    per_chain = pd.read_csv(bm.PER_CHAIN_METRICS_PATH)
    wide = per_chain[per_chain["label"] == "is_interface_contact"].pivot_table(
        index=["pdb_id", "chain_id"], columns="input", values="aupr"
    ).dropna()
    wide["aupr_drop"] = wide["exp"] - wide["af_trimmed"]
    labels_report = pd.read_csv(length_source_path)[["pdb_id", "chain_id", "n_observed"]]
    merged = wide.reset_index().merge(labels_report, on=["pdb_id", "chain_id"], how="left").dropna(subset=["n_observed"])

    x, y = merged["n_observed"].to_numpy(dtype=float), merged["aupr_drop"].to_numpy(dtype=float)
    ols = sm.OLS(y, sm.add_constant(x)).fit()
    x_line = np.linspace(x.min(), x.max(), 100)
    y_line = ols.params[0] + ols.params[1] * x_line

    fig, ax = plt.subplots(figsize=(6, 4.5), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.scatter(x, y, s=18, color=bm.COLOR_EXP, alpha=0.5, edgecolor="none")
    ax.plot(x_line, y_line, color=bm.COLOR_AF_TRIMMED, linewidth=2)
    ax.axhline(0, color="#c3c2b7", linewidth=1, linestyle="--")
    ax.set_xlabel("Chain length (observed residues)")
    ax.set_ylabel("Per-chain AUPR drop (exp - af_trimmed)")
    ax.set_title(f"AUPR drop vs. chain length (n={len(merged)})")
    bm._style_axes(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_flagged_vs_unflagged(flagged_summary: pd.DataFrame, path: Path) -> None:
    metrics_to_plot = ["mean_local_rmsd", "mean_plddt", "aupr_drop"]
    fig, axes = plt.subplots(1, len(metrics_to_plot), figsize=(11, 4), facecolor="#fcfcfb")
    for ax, metric in zip(axes, metrics_to_plot):
        row = flagged_summary[flagged_summary["metric"] == metric].iloc[0]
        ax.set_facecolor("#fcfcfb")
        ax.bar(["flagged", "unflagged"], [row["mean_flagged"], row["mean_unflagged"]],
               color=[bm.COLOR_AF_TRIMMED, bm.COLOR_EXP])
        ax.set_title(f"{metric}\n(Mann-Whitney p={row['mannwhitney_p']:.2g})")
        bm._style_axes(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --- Orchestration -----------------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "error_analysis.log"
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(name)s: %(message)s"))
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))
    root.handlers.clear()
    root.addHandler(file_handler)
    root.addHandler(console_handler)


def run() -> None:
    _setup_logging()
    ERROR_ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Secondary structure method: biotite P-SEA (mkdssp not installed, see PLAN.md Phase 8)")

    features, exclusions = load_all_structural_features()
    n_chains = features[["pdb_id", "chain_id"]].drop_duplicates().shape[0]
    logger.info("Built per-residue structural features: %d residues, %d chains, %d chain-level exclusions",
                len(features), n_chains, len(exclusions))
    write_exclusions(exclusions, STRUCTURAL_EXCLUSIONS_PATH)
    logger.info("Wrote %s", STRUCTURAL_EXCLUSIONS_PATH)

    features.to_parquet(PER_RESIDUE_FEATURES_PATH, index=False)
    logger.info("Wrote %s", PER_RESIDUE_FEATURES_PATH)

    band_metrics = compute_band_metrics(features)
    band_metrics.to_csv(BAND_METRICS_PATH, index=False)
    logger.info("Wrote %s", BAND_METRICS_PATH)

    reg_frame, reg_info = build_regression_frame(features)
    logger.info("Regression: dropped %d/%d residues missing plddt/local_rmsd/rsa/prob_shift",
                reg_info["n_dropped_missing"], reg_info["n_total"])
    regression = fit_shift_regression(reg_frame)
    regression.to_csv(REGRESSION_PATH, index=False)
    logger.info("Wrote %s", REGRESSION_PATH)

    length_analysis = compute_length_analysis()
    length_analysis.to_csv(LENGTH_ANALYSIS_PATH, index=False)
    logger.info("Wrote %s", LENGTH_ANALYSIS_PATH)

    flagged_summary = compute_flagged_chains_summary(features)
    flagged_summary.to_csv(FLAGGED_CHAINS_SUMMARY_PATH, index=False)
    logger.info("Wrote %s", FLAGGED_CHAINS_SUMMARY_PATH)

    plot_band_metrics(band_metrics, FIGURES_DIR / "error_analysis_band_metrics.png")
    plot_shift_vs_local_rmsd(features, FIGURES_DIR / "error_analysis_shift_vs_local_rmsd.png")
    plot_aupr_drop_vs_length(LABELS_REPORT_PATH, FIGURES_DIR / "error_analysis_aupr_drop_vs_length.png")
    plot_flagged_vs_unflagged(flagged_summary, FIGURES_DIR / "error_analysis_flagged_vs_unflagged.png")
    logger.info("Wrote figures to %s", FIGURES_DIR)


def main() -> None:
    run()


if __name__ == "__main__":
    main()
