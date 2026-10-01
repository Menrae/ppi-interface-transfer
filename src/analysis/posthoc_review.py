"""Post-hoc analyses requested in peer review of the results paper.

NOT part of any pre-specified analysis plan: these were chosen on
2026-10-01, after every Phase 7/8/8b result was known, in response to a
reviewer's questions. Everything here is exploratory and is labeled as such
in the paper. Nothing here changes a pre-specified result; each analysis
reads already-written outputs.

1. Template window: could an experimental structure have been available as
   a template when its own AlphaFold DB model was generated? Compares each
   chain's PDB release date with the model's creation date (AlphaFold DB API
   metadata), and repeats the primary endpoint on the chains released on or
   after their model's creation date (which cannot have been templates for
   it).
2. Signed error: refits the Phase 8 regression design with
   |p_AF - y| - |p_exp - y| (positive = AlphaFold input further from the
   true label) as the outcome, alongside the original |p_AF - p_exp|, with
   R^2 for both.
3. Paired F1 bootstrap for the Phase 8b filtering check: the change in
   AlphaFold-input F1 and in the exp-AF F1 gap between the smallest and
   largest pLDDT cutoff, with both ends computed on the same resampled
   chains.
4. Pooled-vs-median: pooled AUPR gap with the chains having the largest
   per-chain drops removed, and with each chain's scores replaced by their
   within-chain percentile rank (removing cross-chain calibration
   differences).
5. Flag audit: residue-mapping method/identity of geometrically flagged vs.
   unflagged chains, and a per-chain listing of residues with extreme local
   RMSD.

Runnable as: .venv/bin/python -m src.analysis.posthoc_review
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata
from sklearn.metrics import average_precision_score

from src import config
from src.analysis import benchmark as bm
from src.analysis import error_analysis as ea
from src.analysis import plddt_filtering as pf

logger = logging.getLogger(__name__)

POSTHOC_DIR = config.RESULTS_DIR / "posthoc"
TEMPLATE_WINDOW_PATH = POSTHOC_DIR / "template_window.csv"
TEMPLATE_WINDOW_CHAINS_PATH = POSTHOC_DIR / "template_window_chains.csv"
SIGNED_ERROR_REGRESSION_PATH = POSTHOC_DIR / "signed_error_regression.csv"
SIGNED_ERROR_SUMMARY_PATH = POSTHOC_DIR / "signed_error_summary.csv"
PAIRED_F1_PATH = POSTHOC_DIR / "paired_f1_bootstrap.csv"
POOLED_VS_MEDIAN_PATH = POSTHOC_DIR / "pooled_vs_median.csv"
FLAG_AUDIT_PATH = POSTHOC_DIR / "flag_mapping_audit.csv"
LOCAL_RMSD_AUDIT_PATH = POSTHOC_DIR / "extreme_local_rmsd_chains.csv"

PER_RESIDUE_FEATURES_PATH = config.RESULTS_DIR / "error_analysis" / "per_residue_features.parquet"
PER_CHAIN_METRICS_PATH = bm.PER_CHAIN_METRICS_PATH
AF_API_CACHE_DIR = config.RAW_DATA_DIR / "alphafold" / "api"
FETCH_REPORT_PATH = config.INTERIM_DATA_DIR / "fetch_report.csv"
MAPPING_REPORT_PATH = config.INTERIM_DATA_DIR / "mapping_report.csv"
GEOMETRIC_VALIDATION_PATH = config.INTERIM_DATA_DIR / "phase4_geometric_validation.csv"

CHAIN_KEY = ["pdb_id", "chain_id"]


def _percentile_ci(draws: np.ndarray) -> tuple[float, float]:
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def load_per_chain_drops(path: Path = PER_CHAIN_METRICS_PATH) -> pd.DataFrame:
    """One row per chain with both AUPRs defined: pdb_id, chain_id,
    aupr_exp, aupr_af_trimmed, drop (exp minus af_trimmed)."""
    pc = pd.read_csv(path)
    pc = pc[pc["label"] == "is_interface_contact"]
    wide = pc.pivot_table(index=CHAIN_KEY, columns="input", values="aupr").dropna().reset_index()
    wide = wide.rename(columns={"exp": "aupr_exp", "af_trimmed": "aupr_af_trimmed"})
    wide["drop"] = wide["aupr_exp"] - wide["aupr_af_trimmed"]
    return wide[CHAIN_KEY + ["aupr_exp", "aupr_af_trimmed", "drop"]]


# --- 1. Template window -------------------------------------------------------


def model_created_date(api_body: list | None, accession: str) -> str | None:
    """Latest modelCreatedDate (YYYY-MM-DD) among the API entries for exactly
    this accession (isoform entries are ignored, as in Phase 3). The latest
    date is the conservative choice: a later model could have used more
    templates."""
    dates = [e["modelCreatedDate"][:10] for e in (api_body or [])
             if e.get("uniprotAccession") == accession and e.get("modelCreatedDate")]
    return max(dates) if dates else None


def load_model_created_dates(accessions) -> dict[str, str | None]:
    out = {}
    for acc in accessions:
        path = AF_API_CACHE_DIR / f"{acc}.json"
        body = json.loads(path.read_text()).get("body") if path.exists() else None
        out[acc] = model_created_date(body, acc)
    return out


def classify_template_window(release_date: str, created_date: str | None) -> str:
    """"could_be_template" if the PDB entry was released before its model
    was created (an upper bound: the template database predates creation),
    "after_model" if released on or after, "unknown" if no creation date."""
    if created_date is None or pd.isna(release_date):
        return "unknown"
    return "could_be_template" if str(release_date)[:10] < created_date else "after_model"


def template_window_analysis(chains: pd.DataFrame, drops: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """`chains`: pdb_id, chain_id, release_date, model_created_date (every
    analyzed chain). Returns (summary rows per window class, per-chain table)."""
    chains = chains.copy()
    chains["template_window"] = [classify_template_window(r, c) for r, c in
                                 zip(chains["release_date"], chains["model_created_date"])]
    merged = chains.merge(drops, on=CHAIN_KEY, how="left")
    rows = []
    for window in ("could_be_template", "after_model", "unknown"):
        sub = merged[merged["template_window"] == window]
        diffs = sub["drop"].dropna().to_numpy()
        row = {"template_window": window, "n_chains": len(sub), "n_scoreable": len(diffs)}
        if len(diffs) >= 2:
            row.update(bm.paired_test(diffs))
        rows.append(row)
    n_unknown = int((merged["template_window"] == "unknown").sum())
    if n_unknown:
        logger.warning("template window: %d chains have no model creation date (reported as 'unknown')", n_unknown)
    return pd.DataFrame(rows), merged


# --- 2. Signed error ----------------------------------------------------------


def add_signed_error(features: pd.DataFrame) -> pd.DataFrame:
    """error_delta = |p_AF - y| - |p_exp - y|: positive when the
    AlphaFold-input prediction is further from the true label."""
    df = features.copy()
    y = df["is_interface_contact"].astype(float)
    df["error_delta"] = (df["af_trimmed_prob"] - y).abs() - (df["exp_prob"] - y).abs()
    return df


def signed_error_analysis(features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = add_signed_error(features)
    cols = CHAIN_KEY + ["plddt", "local_rmsd", "rsa", "secondary_structure", "prob_shift", "error_delta"]
    df = df[cols]
    missing = df[["plddt", "local_rmsd", "rsa", "prob_shift", "error_delta"]].isna().any(axis=1)
    logger.info("signed-error regression: dropped %d/%d residues with a missing predictor/outcome",
                int(missing.sum()), len(df))
    df = df.loc[~missing].reset_index(drop=True)

    regression = pd.concat([ea.fit_shift_regression(df, outcome="prob_shift"),
                            ea.fit_shift_regression(df, outcome="error_delta")], ignore_index=True)

    # Residue-level mean error_delta with a chain-resampled CI (same pooled
    # statistic per draw), plus the fraction of residues where AF input is
    # further from the label.
    per_chain = df.groupby(CHAIN_KEY)["error_delta"].agg(["sum", "count"])
    lo, hi = ea.bootstrap_paired_stat_ci(per_chain["sum"].to_numpy(), per_chain["count"].to_numpy(),
                                         lambda s, c: s.sum() / c.sum(),
                                         config.BOOTSTRAP_N_RESAMPLES, config.PHASE8_BOOTSTRAP_SEED)
    summary = pd.DataFrame([{
        "n_residues": len(df), "mean_error_delta": float(df["error_delta"].mean()),
        "ci_lo": lo, "ci_hi": hi,
        "frac_af_worse": float((df["error_delta"] > 0).mean()),
        "frac_af_better": float((df["error_delta"] < 0).mean()),
    }])
    return regression, summary


# --- 3. Paired F1 bootstrap ---------------------------------------------------


def _f1(labels: np.ndarray, probs: np.ndarray, threshold: float) -> float:
    return pf.compute_prf(labels, probs >= threshold)["f1"]


def paired_f1_stats(chain_arrays: list[dict], lo_cutoff: float, hi_cutoff: float, threshold: float) -> dict:
    """F1 for both inputs at both cutoffs on one set of chains, and the
    derived differences. Each chain_arrays item: {labels, exp, af, plddt}."""
    labels = np.concatenate([c["labels"] for c in chain_arrays])
    exp = np.concatenate([c["exp"] for c in chain_arrays])
    af = np.concatenate([c["af"] for c in chain_arrays])
    plddt = np.concatenate([c["plddt"] for c in chain_arrays])
    out = {}
    for name, cutoff in (("lo", lo_cutoff), ("hi", hi_cutoff)):
        keep = plddt >= cutoff
        out[f"f1_exp_{name}"] = _f1(labels[keep], exp[keep], threshold)
        out[f"f1_af_{name}"] = _f1(labels[keep], af[keep], threshold)
    out["af_gain"] = out["f1_af_hi"] - out["f1_af_lo"]
    out["exp_change"] = out["f1_exp_hi"] - out["f1_exp_lo"]
    out["gap_lo"] = out["f1_exp_lo"] - out["f1_af_lo"]
    out["gap_hi"] = out["f1_exp_hi"] - out["f1_af_hi"]
    out["gap_change"] = out["gap_hi"] - out["gap_lo"]
    return out


def chain_arrays_from_features(features: pd.DataFrame) -> list[dict]:
    return [
        {"labels": g["is_interface_contact"].to_numpy(dtype=bool), "exp": g["exp_prob"].to_numpy(),
         "af": g["af_trimmed_prob"].to_numpy(), "plddt": g["plddt"].to_numpy()}
        for _, g in features.groupby(CHAIN_KEY, sort=True)
    ]


def paired_f1_bootstrap(features: pd.DataFrame, n_resamples: int, seed: int) -> pd.DataFrame:
    lo_cutoff, hi_cutoff = min(config.PLDDT_FILTER_CUTOFFS), max(config.PLDDT_FILTER_CUTOFFS)
    threshold = config.PLDDT_FILTER_PREDICTION_THRESHOLD
    chains = chain_arrays_from_features(features)
    point = paired_f1_stats(chains, lo_cutoff, hi_cutoff, threshold)
    rng = np.random.default_rng(seed)
    draws = {k: [] for k in point}
    n_undefined = 0
    for _ in range(n_resamples):
        idx = rng.integers(0, len(chains), size=len(chains))
        stats = paired_f1_stats([chains[i] for i in idx], lo_cutoff, hi_cutoff, threshold)
        if any(np.isnan(v) for v in stats.values()):
            n_undefined += 1
            continue
        for k, v in stats.items():
            draws[k].append(v)
    if n_undefined:
        logger.warning("paired F1 bootstrap: %d/%d draws undefined, excluded", n_undefined, n_resamples)
    rows = []
    for k, v in point.items():
        arr = np.array(draws[k])
        lo, hi = _percentile_ci(arr)
        rows.append({"statistic": k, "value": v, "ci_lo": lo, "ci_hi": hi,
                     "frac_draws_le_zero": float(np.mean(arr <= 0)), "n_draws": len(arr),
                     "n_undefined": n_undefined, "lo_cutoff": lo_cutoff, "hi_cutoff": hi_cutoff})
    return pd.DataFrame(rows)


# --- 4. Pooled vs. median -----------------------------------------------------


def within_chain_rank(features: pd.DataFrame, col: str) -> pd.Series:
    """Each residue's score replaced by its percentile rank within its own
    chain (average ranks for ties, scaled to (0, 1])."""
    return features.groupby(CHAIN_KEY)[col].transform(lambda s: pd.Series(rankdata(s) / len(s), index=s.index))


def top_drop_chains(drops: pd.DataFrame, fraction: float) -> set[tuple[str, str]]:
    n = int(round(fraction * len(drops)))
    top = drops.sort_values("drop", ascending=False).head(n)
    return set(zip(top["pdb_id"], top["chain_id"]))


def pooled_gap(df: pd.DataFrame, exp_col: str, af_col: str) -> dict:
    y = df["is_interface_contact"].to_numpy(dtype=bool)
    exp = average_precision_score(y, df[exp_col].to_numpy())
    af = average_precision_score(y, df[af_col].to_numpy())
    return {"aupr_exp": exp, "aupr_af": af, "gap": exp - af}


def bootstrap_pooled_gap(df: pd.DataFrame, exp_col: str, af_col: str, n_resamples: int, seed: int) -> tuple[float, float]:
    groups = [g for _, g in df.groupby(CHAIN_KEY, sort=True)]
    ys = [g["is_interface_contact"].to_numpy(dtype=bool) for g in groups]
    es = [g[exp_col].to_numpy() for g in groups]
    afs = [g[af_col].to_numpy() for g in groups]
    rng = np.random.default_rng(seed)
    gaps = []
    for _ in range(n_resamples):
        idx = rng.integers(0, len(groups), size=len(groups))
        y = np.concatenate([ys[i] for i in idx])
        if y.all() or not y.any():
            continue
        gaps.append(average_precision_score(y, np.concatenate([es[i] for i in idx]))
                    - average_precision_score(y, np.concatenate([afs[i] for i in idx])))
    return _percentile_ci(np.array(gaps))


def pooled_vs_median_analysis(features: pd.DataFrame, drops: pd.DataFrame, n_resamples: int, seed: int) -> pd.DataFrame:
    df = features[CHAIN_KEY + ["is_interface_contact", "exp_prob", "af_trimmed_prob"]].copy()
    df["exp_rank"] = within_chain_rank(df, "exp_prob")
    df["af_rank"] = within_chain_rank(df, "af_trimmed_prob")
    removed = top_drop_chains(drops, config.POSTHOC_TOP_DROP_FRACTION)
    keys = list(zip(df["pdb_id"], df["chain_id"]))
    trimmed = df[[k not in removed for k in keys]]
    logger.info("pooled-vs-median: removed %d chains with the largest per-chain drops (top %.0f%%)",
                len(removed), 100 * config.POSTHOC_TOP_DROP_FRACTION)

    variants = [
        ("all_chains_raw", df, "exp_prob", "af_trimmed_prob", 0),
        ("top_drop_removed_raw", trimmed, "exp_prob", "af_trimmed_prob", len(removed)),
        ("all_chains_within_chain_rank", df, "exp_rank", "af_rank", 0),
    ]
    rows = []
    for name, data, ecol, acol, n_removed in variants:
        g = pooled_gap(data, ecol, acol)
        lo, hi = bootstrap_pooled_gap(data, ecol, acol, n_resamples, seed)
        rows.append({"variant": name, "n_chains": data[CHAIN_KEY].drop_duplicates().shape[0],
                     "n_residues": len(data), "n_chains_removed": n_removed, **g, "gap_ci_lo": lo, "gap_ci_hi": hi})
    # Reference points from the per-chain distribution, for the same table.
    rows.append({"variant": "per_chain_median_drop", "n_chains": len(drops), "gap": float(drops["drop"].median())})
    rows.append({"variant": "per_chain_mean_drop", "n_chains": len(drops), "gap": float(drops["drop"].mean())})
    remaining = drops[[k not in removed for k in zip(drops["pdb_id"], drops["chain_id"])]]
    rows.append({"variant": "per_chain_mean_drop_top_removed", "n_chains": len(remaining),
                 "n_chains_removed": len(removed), "gap": float(remaining["drop"].mean())})
    return pd.DataFrame(rows)


# --- 5. Flag / extreme-local-RMSD audit ---------------------------------------


def flag_mapping_audit(chains: pd.DataFrame) -> pd.DataFrame:
    """`chains`: one row per scoreable chain with geometrically_flagged,
    mapping_method, mapping_identity."""
    rows = []
    for flagged, g in chains.groupby("geometrically_flagged"):
        rows.append({
            "group": "flagged" if flagged else "unflagged", "n_chains": len(g),
            "n_fallback": int((g["mapping_method"] == "fallback").sum()),
            "n_identity_below_1": int((g["mapping_identity"] < 1).sum()),
            "mean_identity": float(g["mapping_identity"].mean()),
            "min_identity": float(g["mapping_identity"].min()),
        })
    f = chains.loc[chains["geometrically_flagged"], "mapping_identity"]
    u = chains.loc[~chains["geometrically_flagged"], "mapping_identity"]
    p = float(mannwhitneyu(f, u, alternative="two-sided").pvalue) if len(f) and len(u) else float("nan")
    out = pd.DataFrame(rows)
    out["identity_mannwhitney_p"] = p
    return out


def extreme_local_rmsd_audit(features: pd.DataFrame, chain_info: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Per chain with >=1 residue above `threshold` local RMSD: how many
    such residues, plus the chain's flag/mapping/global-fit context."""
    far = features[features["local_rmsd"] > threshold]
    per_chain = far.groupby(CHAIN_KEY).agg(n_extreme=("local_rmsd", "size"),
                                           max_local_rmsd=("local_rmsd", "max")).reset_index()
    sizes = features.groupby(CHAIN_KEY).size().rename("n_residues").reset_index()
    per_chain = per_chain.merge(sizes, on=CHAIN_KEY).merge(chain_info, on=CHAIN_KEY, how="left")
    return per_chain.sort_values("n_extreme", ascending=False)


# --- Orchestration ------------------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    fh = logging.FileHandler(config.LOGS_DIR / "posthoc_review.log")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(name)s: %(message)s"))
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))
    root.handlers.clear()
    root.addHandler(fh)
    root.addHandler(ch)


def load_chain_info() -> pd.DataFrame:
    """Every analyzed (primary-set, labeled) chain: release date, AlphaFold
    model creation date, mapping method/identity, geometric flag."""
    labels = pd.read_csv(bm.LABELS_REPORT_PATH)
    primary = labels[(labels["status"] == "labeled") & labels["eligible_homolog"]][CHAIN_KEY + ["uniprot_acc"]]
    fetch = pd.read_csv(FETCH_REPORT_PATH)[CHAIN_KEY + ["release_date"]]
    mapping = pd.read_csv(MAPPING_REPORT_PATH)[CHAIN_KEY + ["method", "identity"]].rename(
        columns={"method": "mapping_method", "identity": "mapping_identity"})
    geo = pd.read_csv(GEOMETRIC_VALIDATION_PATH)
    geo["geometrically_flagged"] = geo["status"].eq("validated") & (
        geo["longest_far_run"] >= config.GEOMETRIC_FLAG_MIN_RUN_LENGTH)
    geo = geo[CHAIN_KEY + ["geometrically_flagged", "median_distance"]].rename(
        columns={"median_distance": "global_median_ca_distance"})
    info = primary.merge(fetch, on=CHAIN_KEY, how="left").merge(mapping, on=CHAIN_KEY, how="left").merge(
        geo, on=CHAIN_KEY, how="left")
    info["geometrically_flagged"] = info["geometrically_flagged"].fillna(False).astype(bool)
    dates = load_model_created_dates(info["uniprot_acc"].unique())
    info["model_created_date"] = info["uniprot_acc"].map(dates)
    return info


def run() -> None:
    _setup_logging()
    POSTHOC_DIR.mkdir(parents=True, exist_ok=True)
    logger.info("Post-hoc review analyses (exploratory; not pre-specified)")

    features = pd.read_parquet(PER_RESIDUE_FEATURES_PATH)
    drops = load_per_chain_drops()
    info = load_chain_info()
    logger.info("Loaded %d residues, %d chains (%d scoreable)", len(features), len(info), len(drops))

    summary, per_chain = template_window_analysis(info[CHAIN_KEY + ["release_date", "model_created_date"]], drops)
    summary.to_csv(TEMPLATE_WINDOW_PATH, index=False)
    per_chain.to_csv(TEMPLATE_WINDOW_CHAINS_PATH, index=False)
    logger.info("template window:\n%s", summary.to_string(index=False))

    regression, err_summary = signed_error_analysis(features)
    regression.to_csv(SIGNED_ERROR_REGRESSION_PATH, index=False)
    err_summary.to_csv(SIGNED_ERROR_SUMMARY_PATH, index=False)
    logger.info("signed-error regression:\n%s", regression[["outcome", "term", "coef", "t", "p_value", "r_squared"]].to_string(index=False))

    paired = paired_f1_bootstrap(features, config.PHASE8_BAND_BOOTSTRAP_N_RESAMPLES, config.PHASE8_BOOTSTRAP_SEED)
    paired.to_csv(PAIRED_F1_PATH, index=False)
    logger.info("paired F1 bootstrap:\n%s", paired.to_string(index=False))

    pvm = pooled_vs_median_analysis(features, drops, config.PHASE8_BAND_BOOTSTRAP_N_RESAMPLES, config.PHASE8_BOOTSTRAP_SEED)
    pvm.to_csv(POOLED_VS_MEDIAN_PATH, index=False)
    logger.info("pooled vs median:\n%s", pvm.to_string(index=False))

    scoreable = info.merge(drops[CHAIN_KEY], on=CHAIN_KEY, how="inner")
    audit = flag_mapping_audit(scoreable)
    audit.to_csv(FLAG_AUDIT_PATH, index=False)
    logger.info("flag mapping audit:\n%s", audit.to_string(index=False))

    extreme = extreme_local_rmsd_audit(
        features, info[CHAIN_KEY + ["geometrically_flagged", "mapping_method", "mapping_identity",
                                    "global_median_ca_distance"]], config.POSTHOC_LOCAL_RMSD_AUDIT_ANGSTROM)
    extreme.to_csv(LOCAL_RMSD_AUDIT_PATH, index=False)
    logger.info("extreme local RMSD (> %g A): %d residues in %d chains",
                config.POSTHOC_LOCAL_RMSD_AUDIT_ANGSTROM, int(extreme["n_extreme"].sum()), len(extreme))


def main() -> None:
    run()


if __name__ == "__main__":
    main()
