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
    """The relation and value, for use as "$p @@X@@$" -- e.g. "= 0.012",
    "= 3.9\\times10^{-16}", or "< 10^{-10}" below
    config.P_VALUE_REPORTING_FLOOR."""
    if p < config.P_VALUE_REPORTING_FLOOR:
        return rf"< 10^{{{int(np.log10(config.P_VALUE_REPORTING_FLOOR))}}}"
    if p < 1e-3:
        exp = int(np.floor(np.log10(p)))
        mant = p / (10 ** exp)
        return rf"= {mant:.1f}\times10^{{{exp}}}"
    return f"= {p:.3f}"


def fmt_p_dollars(p: float) -> str:
    """Table-cell form: the value alone (bound kept as "<"), in math mode."""
    s = fmt_p(p)
    return f"${s[2:]}$" if s.startswith("= ") else f"${s}$"


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
              "Leakage-free\ncandidates", "Both structures\ndownloaded", "Residue\nmapping OK",
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

# The all-residues arm of the residue-scope split is the primary endpoint
# itself (same chains, same residues, same numbers), so it is shown once, as
# the primary row, rather than duplicated as a "stratum".
DUPLICATE_OF_PRIMARY = ("residue_scope", "all_residues")


def descriptive_strata(strata: pd.DataFrame) -> pd.DataFrame:
    keep = ~((strata["stratum"] == DUPLICATE_OF_PRIMARY[0]) & (strata["value"] == DUPLICATE_OF_PRIMARY[1]))
    return strata[keep]


def fig_phase7_strata(phase7: dict, path: Path) -> None:
    primary = phase7["primary"]
    strata = descriptive_strata(phase7["strata"])

    rows = [("Primary endpoint (all chains, all residues)", primary["n"], primary["median_diff"], primary["ci_lo"], primary["ci_hi"])]
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
        # CI for the same statistic that is plotted (the residue-level mean
        # within the decile), resampling chains: per-chain sums and counts
        # are resampled jointly and the pooled mean recomputed per draw.
        per_chain = g.groupby(["pdb_id", "chain_id"])["prob_shift"].agg(["sum", "count"])
        lo, hi = ea.bootstrap_paired_stat_ci(per_chain["sum"].to_numpy(), per_chain["count"].to_numpy(),
                                             lambda sums, counts: sums.sum() / counts.sum(),
                                             config.BOOTSTRAP_N_RESAMPLES, config.PHASE8_BOOTSTRAP_SEED)
        # Plot each decile at its own median local RMSD, not its interval
        # midpoint: the top decile's interval runs out to the maximum
        # (tens of Angstrom), so its midpoint sits far from where its
        # residues actually are.
        rows.append({"x": g["local_rmsd"].median(), "n": len(g), "mean_shift": g["prob_shift"].mean(),
                     "ci_lo": lo, "ci_hi": hi})
    binned = pd.DataFrame(rows).sort_values("x")

    fig, axes = plt.subplots(1, 2, figsize=(PAGE_WIDTH_IN, 2.5))

    ax = axes[0]
    ax.plot(binned["x"], binned["mean_shift"], color=BLUE, linewidth=1.6, marker="o", markersize=3)
    ax.fill_between(binned["x"], binned["ci_lo"], binned["ci_hi"], color=BLUE, alpha=0.18, linewidth=0)
    ax.set_xscale("log")
    from matplotlib.ticker import FixedLocator, NullFormatter, FormatStrFormatter
    ax.xaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2]))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Local Cα RMSD (Å), decile median (log scale)")
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
        color = GREY if r.p_value >= 0.05 else BLUE
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
    a redundant table (the forest plot, Fig. 3, already shows them all)."""
    strata = descriptive_strata(strata)
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
    # Fixed decimals per metric (not significant figures), so each row reads
    # at one consistent precision.
    decimals = {"aupr_drop": 3, "chain_length": 0, "interface_fraction": 3, "mean_plddt": 1,
                "mean_local_rmsd": 2, "mean_global_ca_distance": 2}
    lines = []
    lines.append(r"\begin{tabular}{lrrr}")
    lines.append(r"\toprule")
    lines.append(r"Metric & Flagged mean (median) & Unflagged mean (median) & $p$ \\")
    lines.append(r"\midrule")
    for _, r in flagged.iterrows():
        d = decimals.get(r["metric"], 3)
        lines.append(rf"{disp.get(r['metric'], r['metric'])} & {r['mean_flagged']:.{d}f} ({r['median_flagged']:.{d}f}) & "
                     rf"{r['mean_unflagged']:.{d}f} ({r['median_unflagged']:.{d}f}) & {fmt_p_dollars(r['mannwhitney_p'])} \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


# =============================================================================
# Provenance, attrition detail, and supplementary material
# =============================================================================

SUPPLEMENT_DIR = PAPER_DIR / "supplementary"
CHAIN_LIST_PATH = SUPPLEMENT_DIR / "table_s1_primary_set_chains.csv"
FETCH_LOG_PATH = config.LOGS_DIR / "fetch_structures.log"
SIFTS_BULK_PATH = config.RAW_DATA_DIR / "sifts" / "pdb_chain_uniprot.tsv.gz"
MMSEQS_BIN = config.EXTERNAL_DIR / "mmseqs" / "bin" / "mmseqs"
PESTO_ROOT = config.EXTERNAL_DIR / "PeSTo"
TODO = r"\todo{%s}"


def _cmd_output(cmd: list[str], cwd: Path | None = None) -> str | None:
    try:
        out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out.stdout.strip() if out.returncode == 0 and out.stdout.strip() else None


def _log_date_range(path: Path, host: str) -> tuple[str, str] | None:
    """First/last date a request to `host` appears in a run log."""
    if not path.exists():
        return None
    dates = [ln[:10] for ln in path.read_text(errors="replace").splitlines() if host in ln and ln[:4].isdigit()]
    return (dates[0], dates[-1]) if dates else None


def load_provenance() -> dict:
    import gzip
    import json
    from src.models import run_pesto

    prov = {}
    mm = _cmd_output([str(MMSEQS_BIN), "version"])
    prov["MMSEQS_VERSION"] = (r"\texttt{" + mm[:7] + "}") if mm else TODO % "MMseqs2 version"
    pesto_commit = _cmd_output(["git", "rev-parse", "HEAD"], cwd=PESTO_ROOT)
    prov["PESTO_COMMIT"] = (r"\texttt{" + pesto_commit[:7] + "}") if pesto_commit else TODO % "PeSTo commit"
    prov["PESTO_CHECKPOINT"] = r"\texttt{" + run_pesto.PESTO_CHECKPOINT_DIR.name.replace("_", r"\_") + "}"

    af_dates = _log_date_range(FETCH_LOG_PATH, "alphafold.ebi.ac.uk")
    pdbe_dates = _log_date_range(FETCH_LOG_PATH, "www.ebi.ac.uk:443 \"GET /pdbe")
    fmt_range = lambda r: (r[0] if r[0] == r[1] else f"{r[0]} to {r[1]}")
    prov["AF_DOWNLOAD_DATE"] = fmt_range(af_dates) if af_dates else TODO % "AlphaFold DB download date"
    prov["PDBE_DOWNLOAD_DATE"] = fmt_range(pdbe_dates) if pdbe_dates else TODO % "PDBe download date"
    sifts_header = None
    if SIFTS_BULK_PATH.exists():
        with gzip.open(SIFTS_BULK_PATH, "rt") as f:
            first = f.readline().strip()
        if first.startswith("#"):
            sifts_header = " ".join(first.lstrip("# ").split()).replace(" | ", "; ")
    prov["SIFTS_BULK_RELEASE"] = sifts_header if sifts_header else TODO % "SIFTS bulk-file release"

    # AlphaFold DB model metadata for the analyzed chains, from the cached
    # prediction-API responses Phase 3 used.
    labels = pd.read_csv(INTERIM / "labels_report.csv")
    primary = labels[(labels["status"] == "labeled") & labels["eligible_homolog"]]
    fetch = pd.read_csv(INTERIM / "fetch_report.csv")
    merged = primary.merge(fetch[["pdb_id", "chain_id", "release_date", "alphafold_model_version"]],
                           on=["pdb_id", "chain_id"], how="left")
    versions = sorted({int(v) for v in merged["alphafold_model_version"].dropna()})
    tools, created = set(), []
    for acc in primary["uniprot_acc"].unique():
        path = config.RAW_DATA_DIR / "alphafold" / "api" / f"{acc}.json"
        if not path.exists():
            continue
        for entry in json.loads(path.read_text()).get("body") or []:
            if entry.get("uniprotAccession") == acc:
                tools.add(entry.get("toolUsed", ""))
                if entry.get("modelCreatedDate"):
                    created.append(entry["modelCreatedDate"][:10])
    prov["AF_MODEL_VERSIONS"] = ", ".join(f"v{v}" for v in versions)
    prov["AF_TOOL_USED"] = "; ".join(sorted(t for t in tools if t)) or TODO % "AlphaFold DB toolUsed"
    prov["AF_MODEL_CREATED_RANGE"] = (f"{min(created)}" if min(created) == max(created)
                                     else f"{min(created)} to {max(created)}") if created else TODO % "model dates"
    prov["PRIMARY_RELEASE_RANGE"] = f"{merged['release_date'].min()} to {merged['release_date'].max()}"

    # Attrition detail, leakage-free ("homolog") candidates only.
    f_h = fetch[fetch["eligible_homolog"]]
    af_fail = f_h["alphafold_failure_reason"].value_counts()
    prov["N_AF_ABSENT"] = fmt_n(af_fail.get("accession_absent_from_alphafold_db", 0))
    prov["N_AF_FRAGMENTED"] = fmt_n(af_fail.get("fragmented_in_alphafold_db", 0))
    prov["N_AF_SEQ_MISMATCH"] = fmt_n(af_fail.get("sequence_mismatch", 0))
    assert af_fail.sum() == len(f_h) - f_h["alphafold_failure_reason"].isna().sum()

    mapping = pd.read_csv(INTERIM / "mapping_report.csv")
    m_h = mapping[mapping["eligible_homolog"]]
    reasons = m_h["exclusion_reason"].fillna("")
    n_nonofficial = int(reasons.str.startswith("alphafold_model_not_official_afdb").sum())
    n_coverage = int(reasons.str.startswith("coverage_below_threshold").sum())
    n_validation = int(reasons.str.startswith("fallback_failed_validation").sum())
    assert n_nonofficial + n_coverage + n_validation == int((m_h["status"] != "mapped").sum())
    prov["N_NONOFFICIAL_AF"] = fmt_n(n_nonofficial)
    nonofficial_providers = sorted({r.split("providerId=")[-1] for r in reasons if "providerId=" in r})
    prov["NONOFFICIAL_PROVIDERS"] = ", ".join(nonofficial_providers)
    prov["N_COVERAGE_EXCL"] = fmt_n(n_coverage)
    prov["N_VALIDATION_EXCL"] = fmt_n(n_validation)
    prov["MIN_MAPPED_COVERAGE"] = fmt_pct(config.MIN_MAPPED_COVERAGE, 0)

    m_final = primary[["pdb_id", "chain_id"]].merge(m_h, on=["pdb_id", "chain_id"], how="left")
    prov["N_SIFTS_MAPPED"] = fmt_n((m_final["method"] == "sifts").sum())
    prov["N_FALLBACK_MAPPED"] = fmt_n((m_final["method"] == "fallback").sum())
    prov["N_NONIDENTICAL"] = fmt_n((m_final["identity"] < 1).sum())
    prov["MIN_MAPPED_IDENTITY"] = fmt_pct(m_final["identity"].min(), 0)

    l_h = labels[labels["eligible_homolog"]]
    l_excl = l_h["exclusion_reason"].value_counts()
    prov["N_NO_PARTNER"] = fmt_n(l_excl.get("no_protein_partner_in_assembly", 0))
    prov["N_ZERO_INTERFACE"] = fmt_n(l_excl.get("zero_interface_residues", 0))
    prov["N_LABEL_EXCL"] = fmt_n((l_h["status"] != "labeled").sum())
    assert l_excl.sum() == (l_h["status"] != "labeled").sum()

    import inspect
    import biotite
    from biotite.structure import sasa as biotite_sasa
    sasa_defaults = inspect.signature(biotite_sasa).parameters  # interface_labels.py uses these defaults
    prov["BIOTITE_VERSION"] = biotite.__version__
    prov["SASA_PROBE"] = f"{sasa_defaults['probe_radius'].default:g}"
    prov["SASA_POINTS"] = fmt_n(sasa_defaults["point_number"].default)
    pymol_meta = sorted((config.EXTERNAL_DIR / "pymol_render_env" / "conda-meta").glob("pymol-open-source-*.json"))
    prov["PYMOL_VERSION"] = (pymol_meta[0].name.split("-")[3]) if pymol_meta else TODO % "PyMOL version"

    prov["GEO_FAR"] =f"{config.GEOMETRIC_FLAG_FAR_DISTANCE_ANGSTROM:g}"
    prov["GEO_RUN"] = f"{config.GEOMETRIC_FLAG_MIN_RUN_LENGTH}"
    prov["CLUSTER_MIN_COVERAGE"] = fmt_pct(config.CLUSTER_MIN_COVERAGE, 0)
    prov["SURFACE_RSA_THRESHOLD"] = f"{config.SURFACE_RSA_THRESHOLD:g}"
    prov["LOCAL_RMSD_MIN_NEIGHBORS"] = f"{config.LOCAL_RMSD_MIN_NEIGHBORS}"
    prov["P_FLOOR"] = rf"10^{{{int(np.log10(config.P_VALUE_REPORTING_FLOOR))}}}"
    prov["INTERFACE_FRACTION_THRESHOLD"] = f"{config.INTERFACE_FRACTION_STRATUM_THRESHOLD:g}"
    prov["PHASE7_SEED"] = f"{config.PHASE7_BOOTSTRAP_SEED}"
    prov["PHASE8_SEED"] = f"{config.PHASE8_BOOTSTRAP_SEED}"
    return {"ctx": prov, "primary_labels": primary, "release": merged}


POSTHOC_DIR = config.RESULTS_DIR / "posthoc"


def load_posthoc() -> dict:
    return {
        "template": pd.read_csv(POSTHOC_DIR / "template_window.csv").set_index("template_window"),
        "regression": pd.read_csv(POSTHOC_DIR / "signed_error_regression.csv"),
        "error_summary": pd.read_csv(POSTHOC_DIR / "signed_error_summary.csv").iloc[0],
        "paired_f1": pd.read_csv(POSTHOC_DIR / "paired_f1_bootstrap.csv").set_index("statistic"),
        "pooled": pd.read_csv(POSTHOC_DIR / "pooled_vs_median.csv").set_index("variant"),
        "flag_audit": pd.read_csv(POSTHOC_DIR / "flag_mapping_audit.csv").set_index("group"),
        "extreme": pd.read_csv(POSTHOC_DIR / "extreme_local_rmsd_chains.csv"),
    }


def table_regression_comparison(regression: pd.DataFrame) -> str:
    """Same design, two outcomes: the pre-specified |p_AF - p_exp| and the
    post-hoc signed error |p_AF - y| - |p_exp - y|."""
    disp = {"const": "Intercept", **REG_TERM_DISPLAY}
    shift = regression[regression["outcome"] == "prob_shift"].set_index("term")
    err = regression[regression["outcome"] == "error_delta"].set_index("term")
    lines = [r"\begin{tabular}{lrrlrrl}", r"\toprule",
             r" & \multicolumn{3}{c}{Shift $|p_{\mathrm{AF}}-p_{\mathrm{exp}}|$} & "
             r"\multicolumn{3}{c}{Signed error $|p_{\mathrm{AF}}-y|-|p_{\mathrm{exp}}-y|$} \\",
             r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
             r"Term & Coef. & $t$ & $p$ & Coef. & $t$ & $p$ \\", r"\midrule"]
    for term in shift.index:
        s, e = shift.loc[term], err.loc[term]
        lines.append(rf"{disp.get(term, term)} & {s['coef']:.4f} & {s['t']:.1f} & {fmt_p_dollars(s['p_value'])} & "
                     rf"{e['coef']:.4f} & {e['t']:.1f} & {fmt_p_dollars(e['p_value'])} \\")
    lines.append(r"\midrule")
    lines.append(rf"$R^2$ & \multicolumn{{3}}{{l}}{{{shift['r_squared'].iloc[0]:.3f}}} & "
                 rf"\multicolumn{{3}}{{l}}{{{err['r_squared'].iloc[0]:.3f}}} \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def posthoc_context(ph: dict, regression_n: int) -> dict:
    ctx = {}
    t = ph["template"]
    for key, row in (("TW_BEFORE", "could_be_template"), ("TW_AFTER", "after_model")):
        r = t.loc[row]
        ctx.update({f"{key}_N": fmt_n(r["n_chains"]), f"{key}_NS": fmt_n(r["n_scoreable"]),
                    f"{key}_MEDIAN": f"{r['median_diff']:.3f}", f"{key}_CI": fmt_ci(r["ci_lo"], r["ci_hi"]),
                    f"{key}_P": fmt_p(r["p_value"])})
    assert int(t.loc["unknown", "n_chains"]) == 0

    reg = ph["regression"]
    shift = reg[reg["outcome"] == "prob_shift"].set_index("term")
    err = reg[reg["outcome"] == "error_delta"].set_index("term")
    assert int(shift["n"].iloc[0]) == regression_n  # same residues as the pre-specified fit
    es = ph["error_summary"]
    ctx.update({
        "REG_COMPARISON_TABLE": table_regression_comparison(reg),
        "SHIFT_R2": f"{shift['r_squared'].iloc[0]:.3f}", "ERR_R2": f"{err['r_squared'].iloc[0]:.3f}",
        "ERR_PLDDT_COEF": f"{err.loc['plddt_z', 'coef']:.4f}", "ERR_PLDDT_P": fmt_p(err.loc["plddt_z", "p_value"]),
        "ERR_RMSD_COEF": f"{err.loc['local_rmsd_z', 'coef']:.4f}", "ERR_RMSD_P": fmt_p(err.loc["local_rmsd_z", "p_value"]),
        "ERR_RSA_COEF": f"{err.loc['rsa_z', 'coef']:.4f}", "ERR_RSA_P": fmt_p(err.loc["rsa_z", "p_value"]),
        "ERR_MEAN": f"{es['mean_error_delta']:.3f}", "ERR_MEAN_CI": fmt_ci(es["ci_lo"], es["ci_hi"]),
        "ERR_FRAC_WORSE": fmt_pct(es["frac_af_worse"], 0),
    })

    pf_ = ph["paired_f1"]
    for key, stat in (("PF_AF_GAIN", "af_gain"), ("PF_EXP_CHANGE", "exp_change"), ("PF_GAP_CHANGE", "gap_change")):
        if key == "PF_GAP_CHANGE":  # the other two point values are already FILTER_* keys
            ctx[key] = f"{pf_.loc[stat, 'value']:.3f}"
        ctx[f"{key}_CI"] = fmt_ci(pf_.loc[stat, "ci_lo"], pf_.loc[stat, "ci_hi"])

    pv = ph["pooled"]
    ctx.update({
        "PV_TRIM_PCT": fmt_pct(config.POSTHOC_TOP_DROP_FRACTION, 0),
        "PV_TRIM_NREMOVED": fmt_n(pv.loc["top_drop_removed_raw", "n_chains_removed"]),
        "PV_TRIM_GAP": f"{pv.loc['top_drop_removed_raw', 'gap']:.3f}",
        "PV_TRIM_CI": fmt_ci(pv.loc["top_drop_removed_raw", "gap_ci_lo"], pv.loc["top_drop_removed_raw", "gap_ci_hi"]),
        "PV_TRIM_MEAN": f"{pv.loc['per_chain_mean_drop_top_removed', 'gap']:.3f}",
        "PV_RANK_GAP": f"{pv.loc['all_chains_within_chain_rank', 'gap']:.3f}",
        "PV_RANK_CI": fmt_ci(pv.loc["all_chains_within_chain_rank", "gap_ci_lo"], pv.loc["all_chains_within_chain_rank", "gap_ci_hi"]),
        "PV_ALL_CI": fmt_ci(pv.loc["all_chains_raw", "gap_ci_lo"], pv.loc["all_chains_raw", "gap_ci_hi"]),
    })

    fa = ph["flag_audit"]
    ex = ph["extreme"]
    ctx.update({
        "FA_FLAG_FALLBACK": fmt_n(fa.loc["flagged", "n_fallback"]), "FA_FLAG_N": fmt_n(fa.loc["flagged", "n_chains"]),
        "FA_UNFLAG_FALLBACK": fmt_n(fa.loc["unflagged", "n_fallback"]), "FA_UNFLAG_N": fmt_n(fa.loc["unflagged", "n_chains"]),
        "FA_FLAG_MUT": fmt_n(fa.loc["flagged", "n_identity_below_1"]), "FA_UNFLAG_MUT": fmt_n(fa.loc["unflagged", "n_identity_below_1"]),
        "FA_FLAG_MEANID": fmt_pct(fa.loc["flagged", "mean_identity"], 1),
        "FA_UNFLAG_MEANID": fmt_pct(fa.loc["unflagged", "mean_identity"], 1),
        "FA_P": fmt_p(fa["identity_mannwhitney_p"].iloc[0]),
        "EX_THRESH": f"{config.POSTHOC_LOCAL_RMSD_AUDIT_ANGSTROM:g}",
        "EX_NRES": fmt_n(ex["n_extreme"].sum()), "EX_NCHAINS": fmt_n(len(ex)),
        "EX_NRES_FLAGGED": fmt_n(ex.loc[ex["geometrically_flagged"], "n_extreme"].sum()),
        "EX_NCHAINS_FLAGGED": fmt_n(ex["geometrically_flagged"].sum()),
        "EX_NCHAINS_SIFTS_EXACT": fmt_n(((ex["mapping_method"] == "sifts") & (ex["mapping_identity"] == 1)).sum()),
    })
    return ctx


def table_release_years(release_year_path: Path, primary_release: pd.DataFrame) -> str:
    """Leakage flag by PDB release year (all non-redundant representatives),
    plus the release-year composition of the analyzed primary set."""
    ry = pd.read_csv(release_year_path)
    analyzed = pd.to_datetime(primary_release["release_date"]).dt.year.value_counts()
    lines = [r"\begin{tabular}{lrrrr}", r"\toprule",
             r"Release year & Representatives & Leakage-free & \% flagged & Analyzed (primary set) \\", r"\midrule"]
    for _, r in ry.iterrows():
        y = int(r["release_year"])
        lines.append(rf"{y} & {fmt_n(r['n_representatives'])} & {fmt_n(r['n_survive_primary_union_30pct'])} & "
                     rf"{r['pct_flagged_primary_union_30pct']:.1f} & {fmt_n(analyzed.get(y, 0))} \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def write_chain_list(primary: pd.DataFrame, per_chain_wide: pd.DataFrame, path: Path) -> None:
    """Supplementary Table S1: every analyzed chain, with whether it enters
    the per-chain (paired) analyses or only the pooled/residue-level ones."""
    geo = pd.read_csv(INTERIM / "phase4_geometric_validation.csv")
    geo["geometrically_flagged"] = geo["status"].eq("validated") & (
        geo["longest_far_run"] >= config.GEOMETRIC_FLAG_MIN_RUN_LENGTH)
    mapping = pd.read_csv(INTERIM / "mapping_report.csv")[["pdb_id", "chain_id", "method", "identity"]]
    out = primary[["pdb_id", "chain_id", "uniprot_acc", "assembly_id", "n_observed",
                   "interface_fraction_distance", "has_homomeric_partner"]].merge(
        mapping, on=["pdb_id", "chain_id"], how="left").merge(
        geo[["pdb_id", "chain_id", "geometrically_flagged"]], on=["pdb_id", "chain_id"], how="left").merge(
        per_chain_wide[["pdb_id", "chain_id", "aupr_exp", "aupr_af_trimmed"]], on=["pdb_id", "chain_id"], how="left")
    out["in_per_chain_analyses"] = out["aupr_exp"].notna()
    out = out.rename(columns={"method": "mapping_method", "identity": "mapping_identity"})
    path.parent.mkdir(parents=True, exist_ok=True)
    out.sort_values(["pdb_id", "chain_id"]).to_csv(path, index=False)


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
        "PRIMARY_CI": fmt_ci(primary["ci_lo"], primary["ci_hi"]), "PRIMARY_P": fmt_p(primary["p_value"]),
        "ROC_MEDIAN": f"{roc['median_diff']:.3f}", "ROC_CI": fmt_ci(roc["ci_lo"], roc["ci_hi"]), "ROC_P": fmt_p(roc["p_value"]),
        "ROBUST_MEDIAN": f"{robustness['median_diff']:.3f}", "ROBUST_CI": fmt_ci(robustness["ci_lo"], robustness["ci_hi"]),
        "ROBUST_P": fmt_p(robustness["p_value"]),
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
        "PLDDT_COEF": f"{regression.loc['plddt_z','coef']:.4f}", "PLDDT_P": fmt_p(regression.loc["plddt_z", "p_value"]),
        "RMSD_COEF": f"{regression.loc['local_rmsd_z','coef']:.4f}", "RMSD_P": fmt_p(regression.loc["local_rmsd_z", "p_value"]),
        "RSA_COEF": f"{regression.loc['rsa_z','coef']:.4f}", "RSA_P": fmt_p(regression.loc["rsa_z", "p_value"]),
        "MAX_VIF": f"{phase8['regression']['vif'].max(skipna=True):.2f}",
        "PLDDT_T": f"{abs(regression.loc['plddt_z', 't']):.1f}", "RMSD_T": f"{abs(regression.loc['local_rmsd_z', 't']):.1f}",
        "RSA_T": f"{regression.loc['rsa_z', 't']:.1f}",
        "RHO": f"{rho_row['value']:.3f}", "RHO_CI": fmt_ci(rho_row["ci_lo"], rho_row["ci_hi"]), "RHO_P": fmt_p(rho_row["p_value"]),
        "UNWEIGHTED_MEAN_DROP": f"{unweighted_mean:.3f}", "UNWEIGHTED_MEDIAN_DROP": f"{unweighted_median:.3f}",
        "WEIGHTED_MEAN_DROP": f"{weighted_mean:.3f}",
        "FLAGGED_TABLE": table_flagged(phase8["flagged"]),
        "FLAGGED_N": fmt_n(flagged.loc["aupr_drop", "n_flagged"]), "UNFLAGGED_N": fmt_n(flagged.loc["aupr_drop", "n_unflagged"]),
        "FLAGGED_PLDDT_P": fmt_p(flagged.loc["mean_plddt", "mannwhitney_p"]),
        "FLAGGED_RMSD_P": fmt_p(flagged.loc["mean_local_rmsd", "mannwhitney_p"]),
    })

    if example_info["rendered_externally"]:
        ex_figure_file = "figures/example_chain_render.png"
        ex_render_note = (
            "rendered with PyMOL (cartoon representation), not a C$\\alpha$ trace; "
            "the AlphaFold model in (c) is Kabsch-superposed onto the experimental "
            "structure and all three panels share one camera orientation; "
            "produced by scripts/render\\_example\\_chain.py, see Code and Data "
            "Availability"
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

    # --- Additions for the review-1 revision --------------------------------
    prov = load_provenance()
    ctx.update(prov["ctx"])

    wide = phase7["per_chain_wide"]
    n_af_better = int((wide["drop"] < 0).sum())
    n_tie = int((wide["drop"] == 0).sum())
    pc = pd.read_csv(BENCHMARK_DIR / "per_chain_metrics.csv")
    pc = pc[(pc["label"] == "is_interface_contact") & (pc["input"] == "exp")]
    excluded = pc[pc["aupr"].isna()]
    assert len(excluded) == phase7["n_excluded"]
    ctx.update({
        "N_AF_BETTER": fmt_n(n_af_better), "PCT_AF_BETTER": fmt_pct(n_af_better / len(wide), 0),
        "N_TIE": fmt_n(n_tie), "N_EXP_BETTER": fmt_n(len(wide) - n_af_better - n_tie),
        "EXCL_RES": fmt_n(excluded["n"].sum()),
        "EXCL_RES_PCT": fmt_pct(excluded["n"].sum() / pooled_exp["n"], 1),
        "EXCL_MIN_N": fmt_n(excluded["n"].min()), "EXCL_MAX_N": fmt_n(excluded["n"].max()),
        "N_DESCRIPTIVE_SUBGROUPS": fmt_n(len(descriptive_strata(phase7["strata"]))),
    })

    def base_rate(b):
        return float(band[band["band"] == b]["base_rate"].iloc[0])

    ctx.update({"BASE_RATE_LOW": f"{base_rate('<50'):.2f}", "BASE_RATE_HIGH": f"{base_rate('>=90'):.2f}",
                "FRAC_POS_HIGH_BAND": fmt_pct(f1_af_max["frac_positives_kept"], 0)})

    gap0 = f1_exp_0["f1"] - f1_af_0["f1"]
    gapmax = f1_exp_max["f1"] - f1_af_max["f1"]
    exp_decline = f1_exp_0["f1"] - f1_exp_max["f1"]
    ctx.update({
        "FILTER_GAP_NARROW_PCT": fmt_pct((gap0 - gapmax) / gap0, 0),
        "FILTER_EXP_DECLINE": f"{exp_decline:.3f}",
        "FILTER_EXP_SHARE_PCT": fmt_pct(exp_decline / (gap0 - gapmax), 0),
    })

    ctx.update(posthoc_context(load_posthoc(), int(regression.loc["plddt_z", "n"])))
    ctx["RELEASE_YEAR_TABLE"] =table_release_years(INTERIM / "leakage_by_release_year.csv", prov["release"])
    write_chain_list(prov["primary_labels"], wide, CHAIN_LIST_PATH)
    ctx["CHAIN_LIST_FILE"] = r"\url{paper/supplementary/" + CHAIN_LIST_PATH.name + "}"
    remote = _cmd_output(["git", "remote", "get-url", "origin"], cwd=REPO_ROOT)
    ctx["REPO_URL"] = (r"\url{" + remote.removesuffix(".git") + "}") if remote else TODO % "repository URL"
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
