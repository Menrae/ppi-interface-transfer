"""Phase 8b (post-hoc, exploratory -- NOT pre-registered before Phase 7/8's
results were seen; see PLAN.md's "Phase 8b" section for the post-hoc
pre-registration written before these specific numbers were computed).

Tests the claim, made informally in Phase 8's Discussion, that AlphaFold's
pLDDT is an "actionable" signal for flagging unreliable PeSTo predictions:
if a user discards low-pLDDT residues before trusting `af_trimmed`
predictions, do precision/recall/F1 on the remaining residues actually
improve, and how much of the true interface does that cost? The identical
sweep is also run on `exp` input as a reference, since pLDDT is a property
of the AlphaFold model shared by both inputs' rows -- this isolates
`af_trimmed`-specific improvement from the generic effect of discarding
hard-to-classify low-confidence residues.

Runnable as: .venv/bin/python -m src.analysis.plddt_filtering
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src import config

logger = logging.getLogger(__name__)

ERROR_ANALYSIS_DIR = config.RESULTS_DIR / "error_analysis"
FIGURES_DIR = config.RESULTS_DIR / "figures"
PER_RESIDUE_FEATURES_PATH = ERROR_ANALYSIS_DIR / "per_residue_features.parquet"
PLDDT_FILTERING_PATH = ERROR_ANALYSIS_DIR / "plddt_filtering.csv"
FIGURE_PATH = FIGURES_DIR / "error_analysis_plddt_filtering.png"

INPUTS = ("exp", "af_trimmed")


def compute_prf(labels: np.ndarray, preds: np.ndarray) -> dict:
    """Precision/recall/F1 for boolean `labels`/`preds` of equal length.
    Each metric is NaN (not 0) when its denominator is zero, so an
    undefined value is never silently treated as a real 0."""
    tp = int(np.sum(labels & preds))
    fp = int(np.sum(~labels & preds))
    fn = int(np.sum(labels & ~preds))
    precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
    recall = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    if np.isnan(precision) or np.isnan(recall) or (precision + recall) == 0:
        f1 = float("nan")
    else:
        f1 = 2 * precision * recall / (precision + recall)
    return {"precision": precision, "recall": recall, "f1": f1}


def bootstrap_prf_ci(
    per_chain: list[tuple[np.ndarray, np.ndarray]], threshold: float, n_resamples: int, seed: int
) -> tuple[dict[str, list[float]], dict[str, int]]:
    """Chain-resampled bootstrap of precision/recall/F1 on an already
    pLDDT-cutoff-filtered residue subset. `per_chain` is a list of
    (labels, probs) arrays, one per chain that has >=1 surviving residue
    at this cutoff -- chains with zero surviving residues have nothing to
    resample and are simply absent from this list.
    """
    metrics = ("precision", "recall", "f1")
    draws: dict[str, list[float]] = {m: [] for m in metrics}
    n_undefined = {m: 0 for m in metrics}
    n_chains = len(per_chain)
    if n_chains == 0:
        return draws, n_undefined

    rng = np.random.default_rng(seed)
    for _ in range(n_resamples):
        idx = rng.integers(0, n_chains, size=n_chains)
        labels = np.concatenate([per_chain[i][0] for i in idx])
        probs = np.concatenate([per_chain[i][1] for i in idx])
        preds = probs >= threshold
        result = compute_prf(labels, preds)
        for m in metrics:
            if np.isnan(result[m]):
                n_undefined[m] += 1
            else:
                draws[m].append(result[m])
    return draws, n_undefined


def compute_filtering_curve(
    features: pd.DataFrame,
    cutoffs: tuple[float, ...] = config.PLDDT_FILTER_CUTOFFS,
    threshold: float = config.PLDDT_FILTER_PREDICTION_THRESHOLD,
    n_resamples: int = config.PHASE8_BAND_BOOTSTRAP_N_RESAMPLES,
    seed: int = config.PHASE8_BOOTSTRAP_SEED,
) -> pd.DataFrame:
    total_n = len(features)
    total_positives = int(features["is_interface_contact"].sum())
    if total_n == 0:
        raise ValueError("no residues to compute a filtering curve on")

    rows = []
    for cutoff in cutoffs:
        subset = features[features["plddt"] >= cutoff]
        n_kept = len(subset)
        n_positives_kept = int(subset["is_interface_contact"].sum()) if n_kept else 0
        frac_residues_kept = n_kept / total_n
        frac_positives_kept = n_positives_kept / total_positives if total_positives else float("nan")
        n_chains_kept = subset[["pdb_id", "chain_id"]].drop_duplicates().shape[0] if n_kept else 0

        if n_kept == 0:
            logger.warning("pLDDT cutoff %s discards every residue -- reported as zero rows, not skipped", cutoff)

        for input_name in INPUTS:
            row = {
                "input": input_name, "cutoff": cutoff, "n_residues": n_kept, "n_chains": n_chains_kept,
                "frac_residues_kept": frac_residues_kept, "n_positives": n_positives_kept,
                "frac_positives_kept": frac_positives_kept,
            }
            if n_kept == 0:
                for m in ("precision", "recall", "f1"):
                    row[m] = float("nan")
                    row[f"{m}_ci_lo"] = float("nan")
                    row[f"{m}_ci_hi"] = float("nan")
                    row[f"{m}_n_undefined"] = 0
                rows.append(row)
                continue

            labels_arr = subset["is_interface_contact"].to_numpy(dtype=bool)
            probs_arr = subset[f"{input_name}_prob"].to_numpy()
            preds_arr = probs_arr >= threshold
            point = compute_prf(labels_arr, preds_arr)
            row.update(point)

            per_chain = [
                (g["is_interface_contact"].to_numpy(dtype=bool), g[f"{input_name}_prob"].to_numpy())
                for _, g in subset.groupby(["pdb_id", "chain_id"], sort=False)
            ]
            draws, n_undefined = bootstrap_prf_ci(per_chain, threshold, n_resamples, seed)
            for m in ("precision", "recall", "f1"):
                vals = draws[m]
                row[f"{m}_ci_lo"] = float(np.percentile(vals, 2.5)) if vals else float("nan")
                row[f"{m}_ci_hi"] = float(np.percentile(vals, 97.5)) if vals else float("nan")
                row[f"{m}_n_undefined"] = n_undefined[m]
            rows.append(row)

    return pd.DataFrame(rows)


# --- Figure ------------------------------------------------------------------

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREY = "#8C8C8C"
VERMILLION = "#D55E00"


def _style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_filtering_curve(curve: pd.DataFrame, path: Path, page_width_in: float = 6.6) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(page_width_in, 2.5))

    ax = axes[0]
    for input_name, color, disp in ((("exp", BLUE, "Experimental")), (("af_trimmed", ORANGE, "AlphaFold (trimmed)"))):
        sub = curve[curve["input"] == input_name].sort_values("cutoff")
        ax.plot(sub["cutoff"], sub["f1"], color=color, marker="o", markersize=4, linewidth=1.6, label=disp)
        ax.fill_between(sub["cutoff"], sub["f1_ci_lo"], sub["f1_ci_hi"], color=color, alpha=0.15, linewidth=0)
    ax.set_xlabel("pLDDT cutoff (residues kept if pLDDT ≥ cutoff)")
    ax.set_ylabel("F1 (threshold = 0.5)")
    ax.set_title("(a) F1 vs. pLDDT cutoff")
    ax.legend(frameon=False, fontsize=7)
    _style_axes(ax)

    ax = axes[1]
    ref = curve[curve["input"] == "exp"].sort_values("cutoff")
    ax.plot(ref["cutoff"], 1 - ref["frac_residues_kept"], color=GREY, marker="s", markersize=4,
            linewidth=1.6, label="All residues discarded")
    ax.plot(ref["cutoff"], 1 - ref["frac_positives_kept"], color=VERMILLION, marker="^", markersize=4,
            linewidth=1.6, label="True interface residues discarded")
    ax.set_xlabel("pLDDT cutoff")
    ax.set_ylabel("Fraction discarded")
    ax.set_title("(b) Cost of filtering")
    ax.legend(frameon=False, fontsize=7)
    _style_axes(ax)

    fig.tight_layout(pad=0.6)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# --- Orchestration -------------------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "plddt_filtering.log"
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

    features = pd.read_parquet(
        PER_RESIDUE_FEATURES_PATH,
        columns=["pdb_id", "chain_id", "plddt", "is_interface_contact", "exp_prob", "af_trimmed_prob"],
    )
    logger.info("Loaded %d residues from %s", len(features), PER_RESIDUE_FEATURES_PATH)

    curve = compute_filtering_curve(features)
    curve.to_csv(PLDDT_FILTERING_PATH, index=False)
    logger.info("Wrote %s", PLDDT_FILTERING_PATH)

    for _, r in curve.iterrows():
        logger.info(
            "%-11s cutoff=%3g n=%6d (%.1f%% of residues, %.1f%% of true interface kept) "
            "precision=%.3f recall=%.3f f1=%.3f",
            r["input"], r["cutoff"], r["n_residues"], 100 * r["frac_residues_kept"],
            100 * r["frac_positives_kept"], r["precision"], r["recall"], r["f1"],
        )

    plot_filtering_curve(curve, FIGURE_PATH)
    logger.info("Wrote %s", FIGURE_PATH)


def main() -> None:
    run()


if __name__ == "__main__":
    main()
