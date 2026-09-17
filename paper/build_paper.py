"""Builds paper/results_paper.pdf end to end: regenerates every figure and
table from results/ and data/interim/, fills the LaTeX template with the
computed numbers, and compiles the PDF.

Every number that appears in the paper is computed here from the pipeline's
own output files -- nothing in template.tex is a typed-in result. The only
hand-entered numbers anywhere in this script are literature citations (Yuan
et al.'s ESMFold figures, quoted from docs/proposal.txt) and config
thresholds, which are read from src/config.py rather than retyped.

Usage: .venv/bin/python paper/build_paper.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src import config
from src.analysis import structural_metrics as sme

PAPER_DIR = Path(__file__).resolve().parent
FIGURES_DIR = PAPER_DIR / "figures"
# If this file exists, Figure 6 uses it as-is instead of the matplotlib
# Calpha-trace fallback -- see README.md "Substituting a rendered structure
# figure" for what it should show and how to produce it.
EXAMPLE_RENDER_OVERRIDE = FIGURES_DIR / "example_chain_render.png"
TEMPLATE_PATH = PAPER_DIR / "template.tex"
TEX_OUTPUT_PATH = PAPER_DIR / "results_paper.tex"
PDF_OUTPUT_PATH = PAPER_DIR / "results_paper.pdf"

BENCHMARK_DIR = config.RESULTS_DIR / "benchmark"
ERROR_ANALYSIS_DIR = config.RESULTS_DIR / "error_analysis"
INTERIM = config.INTERIM_DATA_DIR

# --- Colorblind-safe palette (Okabe & Ito, 2008) ---------------------------
BLUE = "#0072B2"       # experimental input, throughout
ORANGE = "#E69F00"     # AlphaFold (trimmed) input, throughout
VERMILLION = "#D55E00" # flagged / secondary emphasis
GREEN = "#009E73"      # true / positive
GREY = "#8C8C8C"       # neutral / background
YELLOW = "#F0E442"     # sparing use only (band highlight)
DARK = "#1A1A1A"

PAGE_WIDTH_IN = 6.6  # matches the LaTeX text width, so figure fonts print at true size


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#4d4d4d")
    ax.spines["bottom"].set_color("#4d4d4d")
    ax.tick_params(colors="#333333", labelsize=8)
    ax.xaxis.label.set_color("#111111")
    ax.yaxis.label.set_color("#111111")
    ax.xaxis.label.set_fontsize(9)
    ax.yaxis.label.set_fontsize(9)
    ax.title.set_color("#111111")
    ax.title.set_fontsize(9.5)


plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


# --- Formatting helpers ------------------------------------------------------


def fmt_p(p: float) -> str:
    if p < 1e-3:
        exp = int(np.floor(np.log10(p)))
        mant = p / (10 ** exp)
        return rf"{mant:.1f}\times10^{{{exp}}}"
    return f"{p:.3f}"


def fmt_p_dollars(p: float) -> str:
    return f"${fmt_p(p)}$"


def fmt_ci(lo: float, hi: float, dec: int = 3) -> str:
    return f"[{lo:.{dec}f}, {hi:.{dec}f}]"


def fmt_n(n) -> str:
    return f"{int(n):,}"


def fmt_pct(x: float, dec: int = 0) -> str:
    return f"{100 * x:.{dec}f}\\%"


# =============================================================================
# Data loading
# =============================================================================


def load_attrition() -> dict:
    cand = pd.read_csv(INTERIM / "candidates.csv")
    dedup = pd.read_csv(INTERIM / "candidates_dedup.csv")
    p3 = pd.read_csv(INTERIM / "phase3_attrition.csv")
    p4 = pd.read_csv(INTERIM / "phase4_attrition.csv")
    p5 = pd.read_csv(INTERIM / "phase5_attrition.csv")
    p6 = pd.read_csv(INTERIM / "phase6_attrition.csv")

    def row(df, col):
        return int(df.loc[df["leakage_filter_mode"] == "homolog", col].iloc[0])

    stages = {
        "entries": cand["pdb_id"].nunique(),
        "chains": len(cand),
        "clusters": len(dedup),
        "leakage_free": int(dedup["eligible_homolog"].sum()),
        "structures": row(p3, "n_both_success"),
        "mapped": row(p4, "n_mapped"),
        "labeled": row(p5, "n_labeled"),
        "predicted": row(p6, "n_succeeded"),
    }
    overlap_frac = float(dedup["pesto_homolog_overlap"].mean())
    return {"stages": stages, "overlap_frac": overlap_frac}


def load_phase7() -> dict:
    summary = pd.read_csv(BENCHMARK_DIR / "summary.csv")
    strata = pd.read_csv(BENCHMARK_DIR / "strata.csv")
    exclusions = pd.read_csv(BENCHMARK_DIR / "exclusions.csv")
    per_chain = pd.read_csv(BENCHMARK_DIR / "per_chain_metrics.csv")

    def get(endpoint, label="distance"):
        r = summary[(summary["endpoint"] == endpoint) & (summary["label"] == label)]
        return r.iloc[0]

    primary = get("primary_paired_aupr_diff")
    roc = get("secondary_paired_roc_auc_diff")
    robustness = get("primary_paired_aupr_diff_ROBUSTNESS", label="sasa")
    pooled_exp = get("pooled_exp")
    pooled_af = get("pooled_af_trimmed")

    wide = per_chain[per_chain["label"] == "is_interface_contact"].pivot_table(
        index=["pdb_id", "chain_id"], columns="input", values=["aupr", "n"]
    )
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.dropna(subset=["aupr_exp", "aupr_af_trimmed"]).reset_index()
    wide["drop"] = wide["aupr_exp"] - wide["aupr_af_trimmed"]

    n_total_chains = per_chain[["pdb_id", "chain_id"]].drop_duplicates().shape[0]
    n_excluded = n_total_chains - len(wide)

    return {
        "primary": primary, "roc": roc, "robustness": robustness,
        "pooled_exp": pooled_exp, "pooled_af": pooled_af,
        "strata": strata, "exclusions": exclusions, "per_chain_wide": wide,
        "n_total_chains": n_total_chains, "n_excluded": n_excluded,
    }


def load_phase8() -> dict:
    band = pd.read_csv(ERROR_ANALYSIS_DIR / "band_metrics.csv")
    regression = pd.read_csv(ERROR_ANALYSIS_DIR / "regression.csv")
    length = pd.read_csv(ERROR_ANALYSIS_DIR / "length_analysis.csv")
    flagged = pd.read_csv(ERROR_ANALYSIS_DIR / "flagged_chains_summary.csv")
    return {"band": band, "regression": regression, "length": length, "flagged": flagged}


def load_plddt_filtering() -> pd.DataFrame:
    return pd.read_csv(ERROR_ANALYSIS_DIR / "plddt_filtering.csv")


def select_example_chain(per_chain_wide: pd.DataFrame, median_drop: float, seed: int = 0) -> dict:
    """Seeded selection of a chain with an AUPR drop close to Phase 7's
    primary-endpoint median, restricted to a length range that renders
    legibly, for the structure figure."""
    candidates = per_chain_wide[(per_chain_wide["n_exp"] >= 80) & (per_chain_wide["n_exp"] <= 220)].copy()
    candidates["dist"] = (candidates["drop"] - median_drop).abs()
    candidates = candidates.sort_values("dist")
    top_k = candidates.head(8)
    chosen = top_k.sample(n=1, random_state=seed).iloc[0]
    return {
        "pdb_id": chosen["pdb_id"], "chain_id": chosen["chain_id"],
        "n": int(chosen["n_exp"]), "aupr_exp": float(chosen["aupr_exp"]),
        "aupr_af": float(chosen["aupr_af_trimmed"]), "drop": float(chosen["drop"]),
    }


def load_chain_uniprot_acc(pdb_id: str, chain_id: str) -> str:
    labels_report = pd.read_csv(INTERIM / "labels_report.csv")
    row = labels_report[(labels_report["pdb_id"] == pdb_id) & (labels_report["chain_id"] == chain_id)].iloc[0]
    return row["uniprot_acc"]


# =============================================================================
# Figure 1: dataset attrition flow
# =============================================================================


def fig_attrition_flow(stages: dict, path: Path) -> None:
    labels = ["PDB\nentries", "Candidate\nchains", "Non-redundant\nrepresentatives",
              "Leakage-free\n(primary set)", "Both structures\ndownloaded", "Residue\nmapping OK",
              "Interface\nlabeled", "PeSTo\npredicted"]
    keys = ["entries", "chains", "clusters", "leakage_free", "structures", "mapped", "labeled", "predicted"]
    values = [stages[k] for k in keys]

    fig, ax = plt.subplots(figsize=(PAGE_WIDTH_IN, 2.15))
    n = len(labels)
    box_w, box_h = 0.86, 0.62
    xs = np.arange(n)
    for i, (x, lab, val) in enumerate(zip(xs, labels, values)):
        color = BLUE if i < 3 else (VERMILLION if i == 3 else ORANGE)
        box = FancyBboxPatch((x - box_w / 2, -box_h / 2), box_w, box_h,
                              boxstyle="round,pad=0.02,rounding_size=0.06",
                              linewidth=1.0, edgecolor="#333333", facecolor=color, alpha=0.85, zorder=2)
        ax.add_patch(box)
        ax.text(x, 0.10, f"{val:,}", ha="center", va="center", fontsize=8.3, fontweight="bold", color="white", zorder=3)
        ax.text(x, -0.20, lab, ha="center", va="center", fontsize=6.7, color="white", zorder=3, linespacing=1.25)
        if i > 0:
            arrow = FancyArrowPatch((xs[i - 1] + box_w / 2, 0), (x - box_w / 2, 0),
                                     arrowstyle="-|>", mutation_scale=9, linewidth=1.0, color="#333333", zorder=1)
            ax.add_patch(arrow)

    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(-0.55, 0.55)
    ax.axis("off")
    fig.tight_layout(pad=0.3)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 2: Phase 7 primary result (paired AUPR scatter + pooled PR curves)
# =============================================================================


def fig_phase7_primary(per_chain_wide: pd.DataFrame, phase7: dict, path: Path) -> None:
    from sklearn.metrics import average_precision_score, precision_recall_curve

    # Rebuild pooled residue-level arrays for the PR curve directly from the
    # same joined chain frames benchmark.py uses, so the curve is not a
    # copy of a cached image.
    from src.analysis import benchmark as bm
    chain_frames, _ = bm.load_all_chain_frames()

    fig, axes = plt.subplots(1, 2, figsize=(PAGE_WIDTH_IN, 2.6))

    ax = axes[0]
    ax.plot([0, 1], [0, 1], color="#BBBBBB", linewidth=1.2, linestyle="--", zorder=1)
    ax.scatter(per_chain_wide["aupr_af_trimmed"], per_chain_wide["aupr_exp"], s=10,
               color=BLUE, alpha=0.45, edgecolor="none", zorder=2)
    ax.set_xlabel("AlphaFold (trimmed) AUPR")
    ax.set_ylabel("Experimental AUPR")
    ax.set_title(f"(a) Per-chain AUPR (n={len(per_chain_wide)})")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    style_axes(ax)

    ax = axes[1]
    for input_name, color, disp in (("exp", BLUE, "Experimental"), ("af_trimmed", ORANGE, "AlphaFold (trimmed)")):
        labels_arr = np.concatenate([cf.df["is_interface_contact"].to_numpy(dtype=bool) for cf in chain_frames])
        probs_arr = np.concatenate([cf.df[f"{input_name}_prob"].to_numpy() for cf in chain_frames])
        precision, recall, _ = precision_recall_curve(labels_arr, probs_arr)
        aupr = average_precision_score(labels_arr, probs_arr)
        ax.plot(recall, precision, color=color, linewidth=1.6, label=f"{disp} ({aupr:.3f})")
    base_rate = phase7["pooled_exp"]["base_rate"]
    ax.axhline(base_rate, color="#BBBBBB", linewidth=1.0, linestyle="--", label=f"base rate ({base_rate:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("(b) Pooled precision-recall")
    ax.legend(frameon=False, loc="upper right", fontsize=6.6)
    style_axes(ax)

    fig.tight_layout(pad=0.6)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 3: Phase 7 strata forest plot
# =============================================================================

STRATUM_DISPLAY = {
    ("partner_type", "heteromeric"): "Heteromeric partner",
    ("partner_type", "homomeric"): "Homomeric partner",
    ("interface_size", "small"): "Interface fraction < 0.5",
    ("interface_size", "large"): "Interface fraction $\\geq$ 0.5",
    ("geometric_validation", "flagged"): "Geometrically flagged",
    ("geometric_validation", "unflagged"): "Geometrically unflagged",
    ("residue_scope", "all_residues"): "All residues",
    ("residue_scope", "surface_only"): "Surface residues only",
}


def fig_phase7_strata(phase7: dict, path: Path) -> None:
    primary = phase7["primary"]
    strata = phase7["strata"]

    rows = [("Primary endpoint (all chains)", primary["n"], primary["median_diff"], primary["ci_lo"], primary["ci_hi"])]
    for _, r in strata.iterrows():
        key = (r["stratum"], r["value"])
        rows.append((STRATUM_DISPLAY.get(key, f"{r['stratum']}={r['value']}"), int(r["n"]), r["median_diff"], r["ci_lo"], r["ci_hi"]))

    rows = rows[::-1]  # top-to-bottom reading order in the plot
    all_lo = min(r[3] for r in rows)
    all_hi = max(r[4] for r in rows)
    span = all_hi - all_lo
    xlim = (min(0, all_lo) - 0.06 * span, all_hi + 0.28 * span)

    fig, ax = plt.subplots(figsize=(PAGE_WIDTH_IN, 2.75))
    ys = np.arange(len(rows))
    for y, (label, n, med, lo, hi) in zip(ys, rows):
        color = VERMILLION if label.startswith("Primary") else BLUE
        ax.plot([lo, hi], [y, y], color=color, linewidth=1.6, zorder=1)
        ax.scatter([med], [y], color=color, s=22, zorder=2)
        ax.text(xlim[1] - 0.02 * span, y, f"n={n}", va="center", ha="left", fontsize=6.6, color="#333333")
    ax.axvline(0, color="#999999", linewidth=1.0, linestyle="--", zorder=0)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.6)
    ax.set_xlabel("Paired AUPR difference (experimental − AlphaFold trimmed), median with 95% CI")
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)
    fig.tight_layout(pad=0.5)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 4: Phase 8 pLDDT-band metrics
# =============================================================================


def fig_phase8_bands(band: pd.DataFrame, path: Path) -> None:
    band_order = sme.plddt_band_labels()
    fig, axes = plt.subplots(1, 2, figsize=(PAGE_WIDTH_IN, 2.5))
    width = 0.34
    x = np.arange(len(band_order))
    for ax, metric, disp in zip(axes, ("aupr", "roc_auc"), ("AUPR", "ROC-AUC")):
        for offset, input_name, color, lab in ((-width / 2, "exp", BLUE, "Experimental"), (width / 2, "af_trimmed", ORANGE, "AlphaFold (trimmed)")):
            sub = band[band["input"] == input_name].set_index("band").reindex(band_order)
            values = sub[metric].to_numpy()
            lo = values - sub[f"{metric}_ci_lo"].to_numpy()
            hi = sub[f"{metric}_ci_hi"].to_numpy() - values
            ax.bar(x + offset, values, width=width, color=color, label=lab,
                   yerr=[lo, hi], capsize=2.5, error_kw={"linewidth": 0.9, "ecolor": "#333333"})
        ax.set_xticks(x)
        ax.set_xticklabels(band_order, fontsize=7.6)
        ax.set_xlabel("AlphaFold pLDDT band")
        ax.set_ylabel(disp)
        ax.set_title(f"Pooled {disp} by pLDDT band")
        ax.set_ylim(0, 1.0)
        style_axes(ax)
    axes[0].legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.92, fontsize=7, loc="upper right")
    fig.tight_layout(pad=0.6)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 5: Phase 8 shift-vs-local-RMSD + regression coefficients
# =============================================================================

REG_TERM_DISPLAY = {
    "plddt_z": "pLDDT (z)", "local_rmsd_z": "Local RMSD (z)", "rsa_z": "RSA (z)",
    "sse_a": "Helix vs. coil", "sse_b": "Strand vs. coil",
}


def fig_phase8_shift_and_regression(features_path: Path, regression: pd.DataFrame, path: Path) -> None:
    features = pd.read_parquet(features_path, columns=["pdb_id", "chain_id", "local_rmsd", "prob_shift"])
    df = features.dropna(subset=["local_rmsd", "prob_shift"]).copy()
    n_bins = 10
    df["rmsd_bin"] = pd.qcut(df["local_rmsd"], n_bins, duplicates="drop")

    from src.analysis import error_analysis as ea
    rows = []
    for interval, g in df.groupby("rmsd_bin", observed=True):
        chain_means = g.groupby(["pdb_id", "chain_id"])["prob_shift"].mean().to_numpy()
        lo, hi = ea.bootstrap_chain_scalar_ci(chain_means, np.mean, config.BOOTSTRAP_N_RESAMPLES, config.PHASE8_BOOTSTRAP_SEED)
        rows.append({"mid": interval.mid, "mean_shift": g["prob_shift"].mean(), "ci_lo": lo, "ci_hi": hi})
    binned = pd.DataFrame(rows).sort_values("mid")

    fig, axes = plt.subplots(1, 2, figsize=(PAGE_WIDTH_IN, 2.5))

    ax = axes[0]
    ax.plot(binned["mid"], binned["mean_shift"], color=BLUE, linewidth=1.6, marker="o", markersize=3)
    ax.fill_between(binned["mid"], binned["ci_lo"], binned["ci_hi"], color=BLUE, alpha=0.18, linewidth=0)
    ax.set_xlabel("Local Cα RMSD (Å), binned")
    ax.set_ylabel(r"Mean $|p_{\mathrm{AF}} - p_{\mathrm{exp}}|$")
    ax.set_title("(a) Prediction shift vs. local RMSD")
    style_axes(ax)

    ax = axes[1]
    reg = regression[regression["term"] != "const"].copy()
    reg["display"] = reg["term"].map(REG_TERM_DISPLAY)
    reg = reg.iloc[::-1]
    ys = np.arange(len(reg))
    for y, r in zip(ys, reg.itertuples()):
        lo = r.coef - 1.96 * r.cluster_robust_se
        hi = r.coef + 1.96 * r.cluster_robust_se
        color = VERMILLION if r.p_value >= 0.05 else BLUE
        ax.plot([lo, hi], [y, y], color=color, linewidth=1.6, zorder=1)
        ax.scatter([r.coef], [y], color=color, s=22, zorder=2)
    ax.axvline(0, color="#999999", linewidth=1.0, linestyle="--", zorder=0)
    ax.set_yticks(ys)
    ax.set_yticklabels(reg["display"], fontsize=7.6)
    ax.set_xlabel(r"Standardized coefficient on $|p_{\mathrm{AF}} - p_{\mathrm{exp}}|$" + "\n(95% Wald CI, cluster-robust)")
    ax.set_title("(b) Shift regression coefficients")
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)

    fig.tight_layout(pad=0.6)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 6: example structure
# =============================================================================


def fig_example_structure(example: dict, path: Path) -> dict:
    pdb_id, chain_id = example["pdb_id"], example["chain_id"]
    uniprot_acc = load_chain_uniprot_acc(pdb_id, chain_id)

    exp_cif = config.RAW_DATA_DIR / "pdb" / f"{pdb_id}_updated.cif.gz"
    af_cif = config.RAW_DATA_DIR / "alphafold" / f"{uniprot_acc}.cif.gz"

    labels = pd.read_parquet(INTERIM / "interface_labels" / f"{pdb_id}_{chain_id}.parquet")
    exp_pred = pd.read_parquet(config.PROCESSED_DATA_DIR / "predictions" / "exp" / f"{pdb_id}_{chain_id}.parquet")
    af_pred = pd.read_parquet(config.PROCESSED_DATA_DIR / "predictions" / "af_trimmed" / f"{pdb_id}_{chain_id}.parquet")

    labels = labels.copy()
    labels["auth_ins_code"] = labels["auth_ins_code"].fillna("")
    exp_pred = exp_pred.copy()
    exp_pred["auth_ins_code"] = exp_pred["auth_ins_code"].fillna("")
    mapped = labels[labels["uniprot_resnum"].notna()].copy()
    mapped["uniprot_resnum"] = mapped["uniprot_resnum"].astype(int)

    merged = mapped.merge(exp_pred, on=["auth_seq_id", "auth_ins_code"], how="inner")
    merged = merged.merge(af_pred, on="uniprot_resnum", how="inner", suffixes=("_exp", "_af"))
    merged = merged.sort_values("auth_seq_id").reset_index(drop=True)

    if EXAMPLE_RENDER_OVERRIDE.exists():
        # A real molecular render has been dropped in externally (see
        # README.md "Substituting a rendered structure figure") -- use it
        # as-is instead of generating the matplotlib fallback. No
        # coordinates/plotting needed in this branch.
        print(f"Using externally provided structure render: {EXAMPLE_RENDER_OVERRIDE}")
        return {"uniprot_acc": uniprot_acc, "n_residues": len(merged), "rendered_externally": True}

    auth_keys = list(zip(merged["auth_seq_id"].astype(int), merged["auth_ins_code"]))
    uniprot_keys = list(merged["uniprot_resnum"].astype(int))

    exp_coords_map = sme.read_ca_coords_by_auth(exp_cif, chain_id, auth_keys)
    af_coords_map = sme.read_ca_coords_by_uniprot(af_cif, uniprot_keys)
    exp_xyz = np.array([exp_coords_map[k] for k in auth_keys])
    af_xyz = np.array([af_coords_map[k] for k in uniprot_keys])
    af_xyz_aligned = sme.kabsch_superpose(af_xyz, exp_xyz)

    is_interface = merged["is_interface_contact"].to_numpy()
    exp_prob = merged["pesto_interface_prob_exp"].to_numpy()
    af_prob = merged["pesto_interface_prob_af"].to_numpy()

    fig = plt.figure(figsize=(PAGE_WIDTH_IN, 2.4))
    elev, azim = 15, 60

    def style3d(ax):
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
            pane.fill = False
            pane.set_edgecolor((1, 1, 1, 0))
        ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.line.set_color((1, 1, 1, 0))
        try:
            ax.set_box_aspect([1, 1, 1])
        except Exception:
            pass

    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax1.plot(exp_xyz[:, 0], exp_xyz[:, 1], exp_xyz[:, 2], color="#AAAAAA", linewidth=0.7, zorder=1)
    colors1 = [GREEN if v else "#CCCCCC" for v in is_interface]
    ax1.scatter(exp_xyz[:, 0], exp_xyz[:, 1], exp_xyz[:, 2], c=colors1, s=11, depthshade=True, zorder=2, linewidths=0)
    ax1.view_init(elev=elev, azim=azim)
    style3d(ax1)
    ax1.set_title("(a) True interface", fontsize=8.5)

    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax2.plot(exp_xyz[:, 0], exp_xyz[:, 1], exp_xyz[:, 2], color="#AAAAAA", linewidth=0.7, zorder=1)
    ax2.scatter(exp_xyz[:, 0], exp_xyz[:, 1], exp_xyz[:, 2], c=exp_prob, cmap="viridis", vmin=0, vmax=1, s=11, depthshade=True, zorder=2, linewidths=0)
    ax2.view_init(elev=elev, azim=azim)
    style3d(ax2)
    ax2.set_title("(b) Experimental", fontsize=8)

    ax3 = fig.add_subplot(1, 3, 3, projection="3d")
    ax3.plot(af_xyz_aligned[:, 0], af_xyz_aligned[:, 1], af_xyz_aligned[:, 2], color="#AAAAAA", linewidth=0.7, zorder=1)
    sc3 = ax3.scatter(af_xyz_aligned[:, 0], af_xyz_aligned[:, 1], af_xyz_aligned[:, 2], c=af_prob, cmap="viridis", vmin=0, vmax=1, s=11, depthshade=True, zorder=2, linewidths=0)
    ax3.view_init(elev=elev, azim=azim)
    style3d(ax3)
    ax3.set_title("(c) AlphaFold (trimmed)", fontsize=8)

    fig.subplots_adjust(wspace=0.55)
    cbar = fig.colorbar(sc3, ax=[ax2, ax3], shrink=0.62, pad=0.03, label="Interface probability")
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label("Interface probability", fontsize=7.5)

    true_patch = mpatches.Patch(color=GREEN, label="True interface residue")
    non_patch = mpatches.Patch(color="#CCCCCC", label="Non-interface residue")
    ax1.legend(handles=[true_patch, non_patch], loc="lower center", bbox_to_anchor=(0.5, -0.16),
               fontsize=6.3, frameon=False, ncol=1)

    fig.savefig(path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return {"uniprot_acc": uniprot_acc, "n_residues": len(merged), "rendered_externally": False}


# =============================================================================
# Figure 7: pLDDT-filtering deployment check (Phase 8b, post-hoc/exploratory)
# =============================================================================


def fig_plddt_filtering(curve: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(PAGE_WIDTH_IN, 2.5))

    ax = axes[0]
    for input_name, color, disp in (("exp", BLUE, "Experimental"), ("af_trimmed", ORANGE, "AlphaFold (trimmed)")):
        sub = curve[curve["input"] == input_name].sort_values("cutoff")
        ax.plot(sub["cutoff"], sub["f1"], color=color, marker="o", markersize=4, linewidth=1.6, label=disp)
        ax.fill_between(sub["cutoff"], sub["f1_ci_lo"], sub["f1_ci_hi"], color=color, alpha=0.15, linewidth=0)
    ax.set_xlabel("pLDDT cutoff (kept if pLDDT ≥ cutoff)")
    ax.set_ylabel("F1 (fixed threshold = 0.5)")
    ax.set_title("(a) F1 vs. pLDDT cutoff")
    ax.legend(frameon=False, fontsize=7)
    style_axes(ax)

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
    style_axes(ax)

    fig.tight_layout(pad=0.6)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Tables (LaTeX, booktabs)
# =============================================================================


def strata_extremes(strata: pd.DataFrame) -> dict:
    """Min/max descriptive-stratum effect, used as prose numbers instead of
    a redundant table (the forest plot, Fig. 3, already shows all eight)."""
    lo_row = strata.loc[strata["median_diff"].idxmin()]
    hi_row = strata.loc[strata["median_diff"].idxmax()]
    return {
        "lo_label": STRATUM_DISPLAY[(lo_row["stratum"], lo_row["value"])],
        "lo_value": float(lo_row["median_diff"]), "lo_n": int(lo_row["n"]),
        "hi_label": STRATUM_DISPLAY[(hi_row["stratum"], hi_row["value"])],
        "hi_value": float(hi_row["median_diff"]), "hi_n": int(hi_row["n"]),
    }


def table_band_metrics(band: pd.DataFrame) -> str:
    band_order = sme.plddt_band_labels()
    lines = []
    lines.append(r"\begin{tabular}{lrrlrlrl}")
    lines.append(r"\toprule")
    lines.append(r"Band & $n$ res. & $n$ ch. & Base rate & Input & ROC-AUC (95\% CI) & \multicolumn{2}{l}{AUPR (95\% CI)} \\")
    lines.append(r"\midrule")
    for b in band_order:
        sub = band[band["band"] == b].set_index("input")
        n_res = int(sub.iloc[0]["n_residues"])
        n_ch = int(sub.iloc[0]["n_chains"])
        base = sub.iloc[0]["base_rate"]
        for i, (input_name, disp) in enumerate((("exp", "Experimental"), ("af_trimmed", "AlphaFold trimmed"))):
            r = sub.loc[input_name]
            b_label = b if i == 0 else ""
            n_res_label = fmt_n(n_res) if i == 0 else ""
            n_ch_label = fmt_n(n_ch) if i == 0 else ""
            base_label = f"{base:.3f}" if i == 0 else ""
            lines.append(rf"{b_label} & {n_res_label} & {n_ch_label} & {base_label} & {disp} & "
                         rf"{r['roc_auc']:.3f} {fmt_ci(r['roc_auc_ci_lo'], r['roc_auc_ci_hi'])} & "
                         rf"\multicolumn{{2}}{{l}}{{{r['aupr']:.3f} {fmt_ci(r['aupr_ci_lo'], r['aupr_ci_hi'])}}} \\")
        if b != band_order[-1]:
            lines.append(r"\addlinespace")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def table_regression(regression: pd.DataFrame) -> str:
    lines = []
    lines.append(r"\begin{tabular}{lrrrrr}")
    lines.append(r"\toprule")
    lines.append(r"Term & Coef. & Cluster-robust SE & $t$ & $p$ & VIF \\")
    lines.append(r"\midrule")
    disp = {"const": "Intercept", **REG_TERM_DISPLAY}
    for _, r in regression.iterrows():
        vif = "--" if pd.isna(r["vif"]) else f"{r['vif']:.2f}"
        lines.append(rf"{disp.get(r['term'], r['term'])} & {r['coef']:.4f} & {r['cluster_robust_se']:.4f} & "
                     rf"{r['t']:.2f} & {fmt_p_dollars(r['p_value'])} & {vif} \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def table_plddt_filtering(curve: pd.DataFrame) -> str:
    lines = []
    lines.append(r"\begin{tabular}{lrrrlll}")
    lines.append(r"\toprule")
    lines.append(r"Input & Cutoff & \% res. kept & \% interface kept & Precision (95\% CI) & Recall (95\% CI) & F1 (95\% CI) \\")
    lines.append(r"\midrule")
    disp = {"exp": "Experimental", "af_trimmed": "AlphaFold trimmed"}
    cutoffs = sorted(curve["cutoff"].unique())
    for c in cutoffs:
        for i, input_name in enumerate(("exp", "af_trimmed")):
            r = curve[(curve["cutoff"] == c) & (curve["input"] == input_name)].iloc[0]
            c_label = f"{c:g}" if i == 0 else ""
            lines.append(
                rf"{disp[input_name]} & {c_label} & {100 * r['frac_residues_kept']:.1f} & "
                rf"{100 * r['frac_positives_kept']:.1f} & "
                rf"{r['precision']:.2f} {fmt_ci(r['precision_ci_lo'], r['precision_ci_hi'], dec=2)} & "
                rf"{r['recall']:.2f} {fmt_ci(r['recall_ci_lo'], r['recall_ci_hi'], dec=2)} & "
                rf"{r['f1']:.2f} {fmt_ci(r['f1_ci_lo'], r['f1_ci_hi'], dec=2)} \\"
            )
        if c != cutoffs[-1]:
            lines.append(r"\addlinespace")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def table_flagged(flagged: pd.DataFrame) -> str:
    disp = {
        "aupr_drop": "Per-chain AUPR drop", "chain_length": "Chain length (residues)",
        "interface_fraction": "Interface fraction", "mean_plddt": "Mean pLDDT",
        "mean_local_rmsd": "Mean local RMSD (\\r{A})", "mean_global_ca_distance": "Mean global C$\\alpha$ dist. (\\r{A})",
    }
    lines = []
    lines.append(r"\begin{tabular}{lrrr}")
    lines.append(r"\toprule")
    lines.append(r"Metric & Flagged mean (median) & Unflagged mean (median) & $p$ \\")
    lines.append(r"\midrule")
    for _, r in flagged.iterrows():
        lines.append(rf"{disp.get(r['metric'], r['metric'])} & {r['mean_flagged']:.3g} ({r['median_flagged']:.3g}) & "
                     rf"{r['mean_unflagged']:.3g} ({r['median_unflagged']:.3g}) & {fmt_p_dollars(r['mannwhitney_p'])} \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


# =============================================================================
# Context assembly + templating
# =============================================================================


def build_context() -> dict:
    attrition = load_attrition()
    phase7 = load_phase7()
    phase8 = load_phase8()
    filtering = load_plddt_filtering()

    example = select_example_chain(phase7["per_chain_wide"], phase7["primary"]["median_diff"], seed=0)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig_attrition_flow(attrition["stages"], FIGURES_DIR / "fig1_attrition.pdf")
    fig_phase7_primary(phase7["per_chain_wide"], phase7, FIGURES_DIR / "fig2_phase7_primary.pdf")
    fig_phase7_strata(phase7, FIGURES_DIR / "fig3_phase7_strata.pdf")
    fig_phase8_bands(phase8["band"], FIGURES_DIR / "fig4_phase8_bands.pdf")
    fig_phase8_shift_and_regression(ERROR_ANALYSIS_DIR / "per_residue_features.parquet", phase8["regression"],
                                     FIGURES_DIR / "fig5_phase8_regression.pdf")
    example_info = fig_example_structure(example, FIGURES_DIR / "fig6_structure_example.png")
    fig_plddt_filtering(filtering, FIGURES_DIR / "fig7_plddt_filtering.pdf")

    primary, roc, robustness = phase7["primary"], phase7["roc"], phase7["robustness"]
    pooled_exp, pooled_af = phase7["pooled_exp"], phase7["pooled_af"]
    band = phase8["band"]
    length = phase8["length"]

    band_low = band[(band["band"] == "<50")]
    band_high = band[(band["band"] == ">=90")]
    gap_low = float(band_low[band_low["input"] == "exp"]["aupr"].iloc[0] - band_low[band_low["input"] == "af_trimmed"]["aupr"].iloc[0])
    gap_high = float(band_high[band_high["input"] == "exp"]["aupr"].iloc[0] - band_high[band_high["input"] == "af_trimmed"]["aupr"].iloc[0])

    rho_row = length[length["row_type"] == "correlation"].iloc[0]
    unweighted_mean = length[length["stratum"] == "unweighted_mean_drop"]["value"].iloc[0]
    unweighted_median = length[length["stratum"] == "unweighted_median_drop"]["value"].iloc[0]
    weighted_mean = length[length["stratum"] == "length_weighted_mean_drop"]["value"].iloc[0]

    regression = phase8["regression"].set_index("term")
    flagged = phase8["flagged"].set_index("metric")

    ctx = {}
    s = attrition["stages"]
    ctx.update({
        "N_ENTRIES": fmt_n(s["entries"]), "N_CHAINS": fmt_n(s["chains"]), "N_CLUSTERS": fmt_n(s["clusters"]),
        "N_LEAKAGE_FREE": fmt_n(s["leakage_free"]), "N_STRUCTURES": fmt_n(s["structures"]),
        "N_MAPPED": fmt_n(s["mapped"]), "N_LABELED": fmt_n(s["labeled"]), "N_PREDICTED": fmt_n(s["predicted"]),
        "OVERLAP_PCT": fmt_pct(attrition["overlap_frac"], 0),
        "RESOLUTION_CUTOFF": f"{config.RESOLUTION_CUTOFF_ANGSTROM:g}",
        "RELEASE_DATE_CUTOFF": config.PDB_RELEASE_DATE_CUTOFF,
        "SEQ_IDENTITY_CUTOFF": fmt_pct(config.SEQUENCE_IDENTITY_CUTOFF, 0),
        "INTERFACE_DISTANCE_CUTOFF": f"{config.INTERFACE_DISTANCE_CUTOFF_ANGSTROM:g}",
        "MIN_CHAIN_LENGTH": f"{config.MIN_CHAIN_LENGTH}",
        "MAX_PROTEIN_ENTITIES": f"{config.MAX_PROTEIN_ENTITIES}",
        "PLDDT_BANDS": ", ".join(str(b) for b in config.PLDDT_BANDS),
        "LOCAL_RMSD_RADIUS": f"{config.LOCAL_RMSD_RADIUS_ANGSTROM:g}",
        "BOOTSTRAP_N": fmt_n(config.BOOTSTRAP_N_RESAMPLES),
        "BAND_BOOTSTRAP_N": fmt_n(config.PHASE8_BAND_BOOTSTRAP_N_RESAMPLES),
        "VIF_THRESHOLD": f"{config.VIF_COLLINEARITY_THRESHOLD:g}",
        "SASA_BURIAL_CUTOFF": f"{config.SASA_BURIAL_CUTOFF_ANGSTROM2:g}",
        "MAPPING_VALIDATION_IDENTITY": fmt_pct(config.RESIDUE_MAPPING_VALIDATION_IDENTITY, 0),
    })

    ctx.update({
        "PRIMARY_N": fmt_n(primary["n"]), "PRIMARY_MEDIAN": f"{primary['median_diff']:.3f}",
        "PRIMARY_CI": fmt_ci(primary["ci_lo"], primary["ci_hi"]), "PRIMARY_P": fmt_p_dollars(primary["p_value"]),
        "ROC_MEDIAN": f"{roc['median_diff']:.3f}", "ROC_CI": fmt_ci(roc["ci_lo"], roc["ci_hi"]), "ROC_P": fmt_p_dollars(roc["p_value"]),
        "ROBUST_MEDIAN": f"{robustness['median_diff']:.3f}", "ROBUST_CI": fmt_ci(robustness["ci_lo"], robustness["ci_hi"]),
        "ROBUST_P": fmt_p_dollars(robustness["p_value"]),
        "POOLED_AUPR_EXP": f"{pooled_exp['aupr']:.3f}", "POOLED_AUPR_AF": f"{pooled_af['aupr']:.3f}",
        "POOLED_ROC_EXP": f"{pooled_exp['roc_auc']:.3f}", "POOLED_ROC_AF": f"{pooled_af['roc_auc']:.3f}",
        "POOLED_N": fmt_n(pooled_exp["n"]), "POOLED_BASE_RATE": f"{pooled_exp['base_rate']:.3f}",
        "POOLED_AUPR_DIFF": f"{pooled_exp['aupr'] - pooled_af['aupr']:.3f}",
        "N_EXCLUDED": fmt_n(phase7["n_excluded"]),
        "EXCLUDED_PCT": fmt_pct(phase7["n_excluded"] / phase7["n_total_chains"], 1),
    })
    strata_ex = strata_extremes(phase7["strata"])
    ctx.update({
        "STRATA_LO_LABEL": strata_ex["lo_label"], "STRATA_LO_VALUE": f"{strata_ex['lo_value']:.3f}", "STRATA_LO_N": fmt_n(strata_ex["lo_n"]),
        "STRATA_HI_LABEL": strata_ex["hi_label"], "STRATA_HI_VALUE": f"{strata_ex['hi_value']:.3f}", "STRATA_HI_N": fmt_n(strata_ex["hi_n"]),
    })

    ctx.update({
        "GAP_LOW": f"{gap_low:.3f}", "GAP_HIGH": f"{gap_high:.3f}", "GAP_RATIO": f"{gap_low / gap_high:.1f}",
        "BAND_TABLE": table_band_metrics(band),
        "REGRESSION_TABLE": table_regression(phase8["regression"]),
        "REG_N": fmt_n(regression.loc["plddt_z", "n"]), "REG_CLUSTERS": fmt_n(regression.loc["plddt_z", "n_clusters"]),
        "REG_N_DROPPED": fmt_n(pooled_exp["n"] - regression.loc["plddt_z", "n"]),
        "PLDDT_COEF": f"{regression.loc['plddt_z','coef']:.4f}", "PLDDT_P": fmt_p_dollars(regression.loc["plddt_z", "p_value"]),
        "RMSD_COEF": f"{regression.loc['local_rmsd_z','coef']:.4f}", "RMSD_P": fmt_p_dollars(regression.loc["local_rmsd_z", "p_value"]),
        "RSA_COEF": f"{regression.loc['rsa_z','coef']:.4f}", "RSA_P": fmt_p_dollars(regression.loc["rsa_z", "p_value"]),
        "MAX_VIF": f"{phase8['regression']['vif'].max(skipna=True):.2f}",
        "RHO": f"{rho_row['value']:.3f}", "RHO_CI": fmt_ci(rho_row["ci_lo"], rho_row["ci_hi"]), "RHO_P": f"{rho_row['p_value']:.2f}",
        "UNWEIGHTED_MEAN_DROP": f"{unweighted_mean:.3f}", "UNWEIGHTED_MEDIAN_DROP": f"{unweighted_median:.3f}",
        "WEIGHTED_MEAN_DROP": f"{weighted_mean:.3f}",
        "FLAGGED_TABLE": table_flagged(phase8["flagged"]),
        "FLAGGED_N": fmt_n(flagged.loc["aupr_drop", "n_flagged"]), "UNFLAGGED_N": fmt_n(flagged.loc["aupr_drop", "n_unflagged"]),
        "FLAGGED_PLDDT_P": fmt_p_dollars(flagged.loc["mean_plddt", "mannwhitney_p"]),
        "FLAGGED_RMSD_P": fmt_p_dollars(flagged.loc["mean_local_rmsd", "mannwhitney_p"]),
    })

    if example_info["rendered_externally"]:
        ex_figure_file = "figures/example_chain_render.png"
        ex_render_note = (
            "externally rendered (see the accompanying repository's README.md "
            "for what this image shows and how it was produced)"
        )
    else:
        ex_figure_file = "figures/fig6_structure_example.png"
        ex_render_note = (
            "no headless 3D molecular renderer could be installed in this "
            "environment, so this is a matplotlib rendering of C$\\alpha$ "
            "positions colored by value, not a rendered surface -- see Methods "
            "discussion in the accompanying repository documentation"
        )

    ctx.update({
        "EX_PDB_ID": str(example["pdb_id"]), "EX_CHAIN_ID": str(example["chain_id"]),
        "EX_UNIPROT": example_info["uniprot_acc"], "EX_N": fmt_n(example["n"]),
        "EX_AUPR_EXP": f"{example['aupr_exp']:.3f}", "EX_AUPR_AF": f"{example['aupr_af']:.3f}",
        "EX_DROP": f"{example['drop']:.3f}",
        "EX_FIGURE_FILE": ex_figure_file, "EX_RENDER_NOTE": ex_render_note,
    })

    min_cutoff, max_cutoff = min(config.PLDDT_FILTER_CUTOFFS), max(config.PLDDT_FILTER_CUTOFFS)

    def filt_row(input_name, cutoff):
        return filtering[(filtering["input"] == input_name) & (filtering["cutoff"] == cutoff)].iloc[0]

    f1_exp_0, f1_af_0 = filt_row("exp", min_cutoff), filt_row("af_trimmed", min_cutoff)
    f1_exp_max, f1_af_max = filt_row("exp", max_cutoff), filt_row("af_trimmed", max_cutoff)

    ctx.update({
        "FILTER_THRESHOLD": f"{config.PLDDT_FILTER_PREDICTION_THRESHOLD:g}",
        "FILTER_MIN_CUTOFF": f"{min_cutoff:g}", "FILTER_MAX_CUTOFF": f"{max_cutoff:g}",
        "FILTER_TABLE": table_plddt_filtering(filtering),
        "FILTER_F1_EXP_0": f"{f1_exp_0['f1']:.3f}",
        "FILTER_F1_AF_0": f"{f1_af_0['f1']:.3f}", "FILTER_F1_AF_0_CI": fmt_ci(f1_af_0["f1_ci_lo"], f1_af_0["f1_ci_hi"]),
        "FILTER_F1_EXP_MAX": f"{f1_exp_max['f1']:.3f}",
        "FILTER_F1_AF_MAX": f"{f1_af_max['f1']:.3f}", "FILTER_F1_AF_MAX_CI": fmt_ci(f1_af_max["f1_ci_lo"], f1_af_max["f1_ci_hi"]),
        "FILTER_GAP_0": f"{f1_exp_0['f1'] - f1_af_0['f1']:.3f}",
        "FILTER_GAP_MAX": f"{f1_exp_max['f1'] - f1_af_max['f1']:.3f}",
        "FILTER_AF_F1_GAIN": f"{f1_af_max['f1'] - f1_af_0['f1']:.3f}",
        "FILTER_AF_PRECISION_0": f"{f1_af_0['precision']:.3f}", "FILTER_AF_PRECISION_MAX": f"{f1_af_max['precision']:.3f}",
        "FILTER_AF_RECALL_0": f"{f1_af_0['recall']:.3f}", "FILTER_AF_RECALL_MAX": f"{f1_af_max['recall']:.3f}",
        "FILTER_FRAC_RES_DISCARDED_MAX": fmt_pct(1 - f1_af_max["frac_residues_kept"], 0),
        "FILTER_FRAC_POS_DISCARDED_MAX": fmt_pct(1 - f1_af_max["frac_positives_kept"], 0),
    })

    return ctx


def render_tex(ctx: dict) -> str:
    template = TEMPLATE_PATH.read_text()
    for key, value in ctx.items():
        token = f"@@{key}@@"
        if token not in template:
            raise ValueError(f"template.tex never uses placeholder {token}")
        template = template.replace(token, str(value))
    remaining = [ln for ln in template.splitlines() if "@@" in ln]
    if remaining:
        raise ValueError(f"unresolved placeholders left in template:\n" + "\n".join(remaining[:10]))
    return template


def compile_pdf() -> None:
    for _ in range(2):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", TEX_OUTPUT_PATH.name],
            cwd=PAPER_DIR, capture_output=True, text=True,
        )
        if result.returncode != 0:
            log_tail = "\n".join(result.stdout.splitlines()[-80:])
            raise RuntimeError(f"pdflatex failed:\n{log_tail}")
    if not PDF_OUTPUT_PATH.exists():
        raise RuntimeError("pdflatex reported success but no PDF was produced")


def main() -> None:
    print("Loading results and generating figures/tables...")
    ctx = build_context()
    print(f"Computed {len(ctx)} template values.")
    tex = render_tex(ctx)
    TEX_OUTPUT_PATH.write_text(tex)
    print(f"Wrote {TEX_OUTPUT_PATH}")
    print("Compiling PDF with pdflatex...")
    compile_pdf()
    print(f"Wrote {PDF_OUTPUT_PATH}")
    # Clean up LaTeX aux files, keep the .tex source and the PDF.
    for ext in (".aux", ".log", ".out"):
        p = PAPER_DIR / f"results_paper{ext}"
        if p.exists():
            p.unlink()


if __name__ == "__main__":
    main()
