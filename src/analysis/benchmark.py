"""Phase 7 (lean): benchmarking PeSTo on experimental vs. trimmed-AlphaFold
inputs, primary set.

Implements exactly the pre-registered analysis plan in PLAN.md ("Phase 7
analysis plan (pre-registered 2026-09-18...)"), committed to before this
module was run against real predictions. Do not add or change endpoints
here without updating that pre-registration first.

Runnable as: .venv/bin/python -m src.analysis.benchmark
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve

from src import config

logger = logging.getLogger(__name__)

# --- Paths -----------------------------------------------------------------

LABELS_REPORT_PATH = config.INTERIM_DATA_DIR / "labels_report.csv"
INTERFACE_LABELS_DIR = config.INTERIM_DATA_DIR / "interface_labels"
PREDICTIONS_DIR = config.PROCESSED_DATA_DIR / "predictions"
GEOMETRIC_VALIDATION_PATH = config.INTERIM_DATA_DIR / "phase4_geometric_validation.csv"

BENCHMARK_DIR = config.RESULTS_DIR / "benchmark"
FIGURES_DIR = config.RESULTS_DIR / "figures"
PER_CHAIN_METRICS_PATH = BENCHMARK_DIR / "per_chain_metrics.csv"
SUMMARY_PATH = BENCHMARK_DIR / "summary.csv"

INPUTS = ("exp", "af_trimmed")
GEOMETRIC_FLAG_MIN_RUN = 5  # same rule as scripts/validate_mapping_geometry.py's MIN_FLAGGED_RUN_LENGTH

# Validated categorical palette slots 1/2 (dataviz skill reference palette),
# used in this fixed order everywhere exp/af_trimmed are plotted together.
COLOR_EXP = "#2a78d6"
COLOR_AF_TRIMMED = "#eb6834"


# --- Data model --------------------------------------------------------


@dataclass
class ChainFrame:
    pdb_id: str
    chain_id: str
    uniprot_acc: str
    df: pd.DataFrame  # one row per shared residue: auth_seq_id, auth_ins_code, uniprot_resnum,
    # is_interface_contact, is_interface_sasa, rsa, exp_prob, af_trimmed_prob


class ResidueSetMismatch(ValueError):
    pass


def assert_same_residue_set(exp_keys: set, af_keys: set) -> None:
    """Raises ResidueSetMismatch unless the two key sets are identical.

    Pure/standalone so it can be unit-tested directly against the "both
    inputs must be evaluated on exactly the same residues" requirement.
    """
    if exp_keys != af_keys:
        raise ResidueSetMismatch(
            f"residue sets differ: {len(exp_keys - af_keys)} only in exp, "
            f"{len(af_keys - exp_keys)} only in af_trimmed"
        )


# --- Building per-chain joined frames ---------------------------------


def build_chain_frame(pdb_id: str, chain_id: str, uniprot_acc: str) -> tuple[ChainFrame | None, str]:
    """Joins residue map + labels + both predictions on exactly the shared
    residue set. Returns (ChainFrame, "") on success or (None, reason)."""
    labels_path = INTERFACE_LABELS_DIR / f"{pdb_id}_{chain_id}.parquet"
    exp_path = PREDICTIONS_DIR / "exp" / f"{pdb_id}_{chain_id}.parquet"
    af_path = PREDICTIONS_DIR / "af_trimmed" / f"{pdb_id}_{chain_id}.parquet"
    if not (labels_path.exists() and exp_path.exists() and af_path.exists()):
        return None, "missing_input_file"

    labels = pd.read_parquet(labels_path)
    labels = labels[labels["uniprot_resnum"].notna()].copy()
    labels["auth_ins_code"] = labels["auth_ins_code"].fillna("")
    labels["uniprot_resnum"] = labels["uniprot_resnum"].astype(int)

    exp = pd.read_parquet(exp_path)[["auth_seq_id", "auth_ins_code", "pesto_interface_prob"]]
    exp["auth_ins_code"] = exp["auth_ins_code"].fillna("")
    exp = exp.rename(columns={"pesto_interface_prob": "exp_prob"})

    af = pd.read_parquet(af_path)[["uniprot_resnum", "pesto_interface_prob"]]
    af["uniprot_resnum"] = af["uniprot_resnum"].astype(int)
    af = af.rename(columns={"pesto_interface_prob": "af_trimmed_prob"})

    exp_keys = set(zip(labels["auth_seq_id"], labels["auth_ins_code"])) & set(
        zip(exp["auth_seq_id"], exp["auth_ins_code"])
    )
    af_keys_uniprot = set(labels["uniprot_resnum"]) & set(af["uniprot_resnum"])
    # Translate af's uniprot-keyed set back to auth keys via the labels table,
    # so both sides of the assertion are expressed in the same key space.
    af_keys = set(
        zip(labels.loc[labels["uniprot_resnum"].isin(af_keys_uniprot), "auth_seq_id"],
            labels.loc[labels["uniprot_resnum"].isin(af_keys_uniprot), "auth_ins_code"])
    )
    full_key_set = set(zip(labels["auth_seq_id"], labels["auth_ins_code"]))

    try:
        assert_same_residue_set(exp_keys, full_key_set)
        assert_same_residue_set(af_keys, full_key_set)
    except ResidueSetMismatch as exc:
        return None, f"residue_set_mismatch:{exc}"

    merged = labels.merge(exp, on=["auth_seq_id", "auth_ins_code"], how="inner")
    merged = merged.merge(af, on="uniprot_resnum", how="inner")
    if len(merged) != len(full_key_set):
        return None, f"join_size_mismatch:expected={len(full_key_set)},got={len(merged)}"

    return ChainFrame(pdb_id, chain_id, uniprot_acc, merged), ""


# --- Metrics -------------------------------------------------------------


def compute_binary_metrics(labels: np.ndarray, probs: np.ndarray) -> dict:
    """{roc_auc, aupr, base_rate, n, n_positive, undefined_reason}.
    roc_auc/aupr are NaN (with a reason) when the label vector is single-class."""
    n = len(labels)
    n_positive = int(np.sum(labels))
    base_rate = n_positive / n if n else float("nan")
    if n_positive == 0 or n_positive == n:
        return {"roc_auc": float("nan"), "aupr": float("nan"), "base_rate": base_rate,
                "n": n, "n_positive": n_positive, "undefined_reason": "single_class"}
    return {
        "roc_auc": roc_auc_score(labels, probs), "aupr": average_precision_score(labels, probs),
        "base_rate": base_rate, "n": n, "n_positive": n_positive, "undefined_reason": "",
    }


def compute_per_chain_metrics(chain_frames: list[ChainFrame], label_col: str = "is_interface_contact") -> pd.DataFrame:
    rows = []
    for cf in chain_frames:
        for input_name in INPUTS:
            metrics = compute_binary_metrics(
                cf.df[label_col].to_numpy(dtype=bool), cf.df[f"{input_name}_prob"].to_numpy()
            )
            rows.append({"pdb_id": cf.pdb_id, "chain_id": cf.chain_id, "uniprot_acc": cf.uniprot_acc,
                         "input": input_name, "label": label_col, **metrics})
    return pd.DataFrame(rows)


def compute_pooled_metrics(chain_frames: list[ChainFrame], label_col: str = "is_interface_contact") -> pd.DataFrame:
    rows = []
    for input_name in INPUTS:
        labels = np.concatenate([cf.df[label_col].to_numpy(dtype=bool) for cf in chain_frames])
        probs = np.concatenate([cf.df[f"{input_name}_prob"].to_numpy() for cf in chain_frames])
        metrics = compute_binary_metrics(labels, probs)
        rows.append({"input": input_name, "label": label_col, "pooled": True, **metrics})
    return pd.DataFrame(rows)


def compute_spearman_and_shift(chain_frames: list[ChainFrame]) -> pd.DataFrame:
    rows = []
    for cf in chain_frames:
        exp_p, af_p = cf.df["exp_prob"].to_numpy(), cf.df["af_trimmed_prob"].to_numpy()
        if len(exp_p) >= 2 and np.std(exp_p) > 0 and np.std(af_p) > 0:
            rho, _ = spearmanr(exp_p, af_p)
        else:
            rho = float("nan")
        rows.append({
            "pdb_id": cf.pdb_id, "chain_id": cf.chain_id,
            "spearman_exp_vs_af_trimmed": rho,
            "mean_abs_prob_shift": float(np.mean(np.abs(exp_p - af_p))),
        })
    return pd.DataFrame(rows)


# --- Bootstrap / paired test -----------------------------------------------


def bootstrap_median_ci(values: np.ndarray, n_resamples: int, seed: int) -> tuple[float, float]:
    """95% percentile bootstrap CI on the median of `values`, resampling
    values (already one-per-chain) with replacement."""
    rng = np.random.default_rng(seed)
    n = len(values)
    medians = np.empty(n_resamples)
    for i in range(n_resamples):
        sample = values[rng.integers(0, n, size=n)]
        medians[i] = np.median(sample)
    return float(np.percentile(medians, 2.5)), float(np.percentile(medians, 97.5))


def paired_test(diffs: np.ndarray) -> dict:
    """Two-sided Wilcoxon signed-rank test + bootstrap CI on the median of
    `diffs` (exp minus af_trimmed, one value per chain)."""
    n = len(diffs)
    median_diff = float(np.median(diffs))
    ci_lo, ci_hi = bootstrap_median_ci(diffs, config.BOOTSTRAP_N_RESAMPLES, config.PHASE7_BOOTSTRAP_SEED)
    if np.all(diffs == 0):
        stat, p = float("nan"), 1.0
    else:
        stat, p = wilcoxon(diffs, alternative="two-sided")
    return {"n": n, "median_diff": median_diff, "ci_lo": ci_lo, "ci_hi": ci_hi,
            "wilcoxon_stat": float(stat), "p_value": float(p)}


# --- Strata ----------------------------------------------------------------


def load_chain_metadata() -> pd.DataFrame:
    labels_report = pd.read_csv(LABELS_REPORT_PATH)[
        ["pdb_id", "chain_id", "has_homomeric_partner", "interface_fraction_distance"]
    ]
    geo = pd.read_csv(GEOMETRIC_VALIDATION_PATH)
    geo["geometrically_flagged"] = geo["status"].eq("validated") & (geo["longest_far_run"] >= GEOMETRIC_FLAG_MIN_RUN)
    geo = geo[["pdb_id", "chain_id", "geometrically_flagged"]]

    meta = labels_report.merge(geo, on=["pdb_id", "chain_id"], how="left")
    meta["geometrically_flagged"] = meta["geometrically_flagged"].fillna(False)
    meta["interface_size_stratum"] = np.where(
        meta["interface_fraction_distance"] >= config.INTERFACE_FRACTION_STRATUM_THRESHOLD, "large", "small"
    )
    return meta


def compute_strata_table(per_chain: pd.DataFrame, meta: pd.DataFrame, label_col: str = "is_interface_contact") -> pd.DataFrame:
    """One row per (stratum, value): paired AUPR difference stats, n."""
    wide = per_chain[per_chain["label"] == label_col].pivot_table(
        index=["pdb_id", "chain_id"], columns="input", values="aupr"
    ).dropna()
    wide["diff"] = wide["exp"] - wide["af_trimmed"]
    wide = wide.merge(meta, on=["pdb_id", "chain_id"], how="left")

    strata = {
        "partner_type": wide["has_homomeric_partner"].map({True: "homomeric", False: "heteromeric"}),
        "interface_size": wide["interface_size_stratum"],
        "geometric_validation": wide["geometrically_flagged"].map({True: "flagged", False: "unflagged"}),
    }
    rows = []
    for stratum_name, values in strata.items():
        for value in sorted(values.dropna().unique()):
            subset = wide.loc[values == value, "diff"].to_numpy()
            if len(subset) == 0:
                continue
            result = paired_test(subset)
            rows.append({"stratum": stratum_name, "value": value, **result})
    return pd.DataFrame(rows)


def compute_surface_vs_all_table(chain_frames: list[ChainFrame], label_col: str = "is_interface_contact") -> pd.DataFrame:
    """The all-residues vs. surface-only stratum needs a different residue
    subset per chain (not a per-chain metadata lookup), so it's computed
    separately: paired AUPR difference, all residues vs. RSA-restricted."""
    rows = []
    for scope, restrict in (("all_residues", False), ("surface_only", True)):
        diffs = []
        for cf in chain_frames:
            df = cf.df
            if restrict:
                df = df[df["rsa"] >= config.SURFACE_RSA_THRESHOLD]
            exp_m = compute_binary_metrics(df[label_col].to_numpy(dtype=bool), df["exp_prob"].to_numpy())
            af_m = compute_binary_metrics(df[label_col].to_numpy(dtype=bool), df["af_trimmed_prob"].to_numpy())
            if np.isnan(exp_m["aupr"]) or np.isnan(af_m["aupr"]):
                continue
            diffs.append(exp_m["aupr"] - af_m["aupr"])
        if diffs:
            result = paired_test(np.array(diffs))
            rows.append({"stratum": "residue_scope", "value": scope, **result})
    return pd.DataFrame(rows)


# --- Figures ---------------------------------------------------------------


def _style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#c3c2b7")
    ax.spines["bottom"].set_color("#c3c2b7")
    ax.tick_params(colors="#52514e")
    ax.xaxis.label.set_color("#0b0b0b")
    ax.yaxis.label.set_color("#0b0b0b")
    ax.title.set_color("#0b0b0b")


def plot_paired_aupr_scatter(per_chain: pd.DataFrame, path: Path, label_col: str = "is_interface_contact") -> None:
    wide = per_chain[per_chain["label"] == label_col].pivot_table(
        index=["pdb_id", "chain_id"], columns="input", values="aupr"
    ).dropna()
    fig, ax = plt.subplots(figsize=(5.5, 5.5), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.plot([0, 1], [0, 1], color="#c3c2b7", linewidth=1.5, linestyle="--", zorder=1)
    ax.scatter(wide["af_trimmed"], wide["exp"], s=22, color=COLOR_EXP, alpha=0.55, edgecolor="none", zorder=2)
    ax.set_xlabel("AlphaFold (trimmed) AUPR")
    ax.set_ylabel("Experimental AUPR")
    ax.set_title(f"Per-chain AUPR, paired (n={len(wide)})")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_diff_distribution(per_chain: pd.DataFrame, path: Path, label_col: str = "is_interface_contact") -> None:
    wide = per_chain[per_chain["label"] == label_col].pivot_table(
        index=["pdb_id", "chain_id"], columns="input", values="aupr"
    ).dropna()
    diffs = wide["exp"] - wide["af_trimmed"]
    fig, ax = plt.subplots(figsize=(6, 4), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.hist(diffs, bins=40, color=COLOR_EXP, edgecolor="white", linewidth=0.5)
    ax.axvline(0, color="#52514e", linewidth=1.5, linestyle="--")
    ax.axvline(diffs.median(), color=COLOR_AF_TRIMMED, linewidth=2, label=f"median = {diffs.median():.3f}")
    ax.set_xlabel("Paired AUPR difference (exp - af_trimmed)")
    ax.set_ylabel("Number of chains")
    ax.set_title(f"Distribution of paired AUPR differences (n={len(diffs)})")
    ax.legend(frameon=False)
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_pooled_curves(chain_frames: list[ChainFrame], pr_path: Path, roc_path: Path, label_col: str = "is_interface_contact") -> None:
    labels_by_input = {}
    for input_name in INPUTS:
        labels = np.concatenate([cf.df[label_col].to_numpy(dtype=bool) for cf in chain_frames])
        probs = np.concatenate([cf.df[f"{input_name}_prob"].to_numpy() for cf in chain_frames])
        labels_by_input[input_name] = (labels, probs)

    colors = {"exp": COLOR_EXP, "af_trimmed": COLOR_AF_TRIMMED}
    display_names = {"exp": "Experimental", "af_trimmed": "AlphaFold (trimmed)"}

    fig, ax = plt.subplots(figsize=(5.5, 5), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    for input_name, (labels, probs) in labels_by_input.items():
        precision, recall, _ = precision_recall_curve(labels, probs)
        aupr = average_precision_score(labels, probs)
        ax.plot(recall, precision, color=colors[input_name], linewidth=2,
                label=f"{display_names[input_name]} (AUPR={aupr:.3f})")
    base_rate = np.mean(labels_by_input["exp"][0])
    ax.axhline(base_rate, color="#c3c2b7", linewidth=1.5, linestyle="--", label=f"base rate ({base_rate:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Pooled precision-recall (all primary-set residues)")
    ax.legend(frameon=False, loc="upper right")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(pr_path, dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 5), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.plot([0, 1], [0, 1], color="#c3c2b7", linewidth=1.5, linestyle="--")
    for input_name, (labels, probs) in labels_by_input.items():
        fpr, tpr, _ = roc_curve(labels, probs)
        auc = roc_auc_score(labels, probs)
        ax.plot(fpr, tpr, color=colors[input_name], linewidth=2,
                label=f"{display_names[input_name]} (ROC-AUC={auc:.3f})")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Pooled ROC (all primary-set residues)")
    ax.legend(frameon=False, loc="lower right")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(roc_path, dpi=150)
    plt.close(fig)


# --- Orchestration -----------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "benchmark.log"
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


def load_all_chain_frames() -> tuple[list[ChainFrame], list[dict]]:
    labels_report = pd.read_csv(LABELS_REPORT_PATH)
    primary = labels_report[(labels_report["status"] == "labeled") & (labels_report["eligible_homolog"])]

    chain_frames, exclusions = [], []
    for rep in primary.itertuples():
        cf, reason = build_chain_frame(rep.pdb_id, rep.chain_id, rep.uniprot_acc)
        if cf is None:
            exclusions.append({"pdb_id": rep.pdb_id, "chain_id": rep.chain_id, "reason": reason})
        else:
            chain_frames.append(cf)
    return chain_frames, exclusions


def run() -> None:
    _setup_logging()
    BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    chain_frames, exclusions = load_all_chain_frames()
    logger.info("Loaded %d chain frames, %d excluded before any metric computation", len(chain_frames), len(exclusions))
    for exc in exclusions:
        logger.warning("excluded %s_%s: %s", exc["pdb_id"], exc["chain_id"], exc["reason"])

    per_chain = compute_per_chain_metrics(chain_frames, "is_interface_contact")
    n_undefined = (per_chain["undefined_reason"] != "").sum()
    logger.info("Per-chain metrics: %d rows, %d undefined (single-class)", len(per_chain), n_undefined)

    spearman_shift = compute_spearman_and_shift(chain_frames)
    per_chain_full = per_chain.merge(spearman_shift, on=["pdb_id", "chain_id"], how="left")
    per_chain_full.to_csv(PER_CHAIN_METRICS_PATH, index=False)
    logger.info("Wrote %s", PER_CHAIN_METRICS_PATH)

    # --- Primary + secondary endpoints (distance labels) ---
    summary_rows = []
    for metric in ("aupr", "roc_auc"):
        wide = per_chain[per_chain["label"] == "is_interface_contact"].pivot_table(
            index=["pdb_id", "chain_id"], columns="input", values=metric
        ).dropna()
        diffs = (wide["exp"] - wide["af_trimmed"]).to_numpy()
        result = paired_test(diffs)
        endpoint = "primary_paired_aupr_diff" if metric == "aupr" else "secondary_paired_roc_auc_diff"
        summary_rows.append({"endpoint": endpoint, "label": "distance", **result})
        logger.info("%s: n=%d median_diff=%.4f CI=[%.4f, %.4f] p=%.4g",
                    endpoint, result["n"], result["median_diff"], result["ci_lo"], result["ci_hi"], result["p_value"])

    pooled = compute_pooled_metrics(chain_frames, "is_interface_contact")
    for row in pooled.to_dict("records"):
        summary_rows.append({"endpoint": f"pooled_{row['input']}", "label": "distance",
                              "n": row["n"], "median_diff": "", "ci_lo": "", "ci_hi": "",
                              "wilcoxon_stat": "", "p_value": "",
                              "roc_auc": row["roc_auc"], "aupr": row["aupr"], "base_rate": row["base_rate"]})

    # --- Robustness: primary endpoint with SASA labels ---
    per_chain_sasa = compute_per_chain_metrics(chain_frames, "is_interface_sasa")
    wide_sasa = per_chain_sasa.pivot_table(index=["pdb_id", "chain_id"], columns="input", values="aupr").dropna()
    diffs_sasa = (wide_sasa["exp"] - wide_sasa["af_trimmed"]).to_numpy()
    result_sasa = paired_test(diffs_sasa)
    summary_rows.append({"endpoint": "primary_paired_aupr_diff_ROBUSTNESS", "label": "sasa", **result_sasa})
    logger.info("robustness (SASA labels): n=%d median_diff=%.4f CI=[%.4f, %.4f] p=%.4g",
                result_sasa["n"], result_sasa["median_diff"], result_sasa["ci_lo"], result_sasa["ci_hi"], result_sasa["p_value"])

    # --- Strata ---
    meta = load_chain_metadata()
    strata_table = compute_strata_table(per_chain, meta, "is_interface_contact")
    surface_table = compute_surface_vs_all_table(chain_frames, "is_interface_contact")
    strata_full = pd.concat([strata_table, surface_table], ignore_index=True)
    strata_full.to_csv(BENCHMARK_DIR / "strata.csv", index=False)
    logger.info("Wrote %s", BENCHMARK_DIR / "strata.csv")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(SUMMARY_PATH, index=False)
    logger.info("Wrote %s", SUMMARY_PATH)

    write_csv_exclusions(exclusions, BENCHMARK_DIR / "exclusions.csv")

    # --- Figures ---
    plot_paired_aupr_scatter(per_chain, FIGURES_DIR / "paired_aupr_scatter.png")
    plot_diff_distribution(per_chain, FIGURES_DIR / "aupr_diff_distribution.png")
    plot_pooled_curves(chain_frames, FIGURES_DIR / "pooled_pr_curve.png", FIGURES_DIR / "pooled_roc_curve.png")
    logger.info("Wrote figures to %s", FIGURES_DIR)


def write_csv_exclusions(exclusions: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pdb_id", "chain_id", "reason"])
        writer.writeheader()
        writer.writerows(exclusions)


def main() -> None:
    run()


if __name__ == "__main__":
    main()
