"""Part C sanity check (not the benchmark): pooled ROC-AUC of PeSTo's
experimental-input predictions against Phase 5's distance-based interface
labels, on the primary set.

This exists purely to confirm the interface channel and the residue joins
are right before any real benchmarking happens -- PeSTo's paper reports
high protein-interface accuracy on experimental structures, so a pooled
AUC anywhere near 0.5 here means something upstream (channel index, key
join, alt-loc handling) is wrong and must be fixed before Phase 7.

Usage: .venv/bin/python scripts/phase6_sanity_check.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import config
from src.models.run_pesto import output_path_for, is_valid_output

LABELS_REPORT_PATH = config.INTERIM_DATA_DIR / "labels_report.csv"
INTERFACE_LABELS_DIR = config.INTERIM_DATA_DIR / "interface_labels"


def main() -> None:
    report = pd.read_csv(LABELS_REPORT_PATH)
    primary = report[(report["status"] == "labeled") & (report["eligible_homolog"])]

    all_labels, all_probs = [], []
    per_chain_auc = []
    n_no_exp_prediction = 0
    n_single_class = 0

    for rep in primary.itertuples():
        exp_path = output_path_for("exp", rep.pdb_id, rep.chain_id)
        if not is_valid_output(exp_path):
            n_no_exp_prediction += 1
            continue

        labels_df = pd.read_parquet(INTERFACE_LABELS_DIR / f"{rep.pdb_id}_{rep.chain_id}.parquet")
        preds_df = pd.read_parquet(exp_path)
        merged = labels_df.merge(preds_df, on=["auth_seq_id", "auth_ins_code"], how="inner")

        all_labels.extend(merged["is_interface_contact"].tolist())
        all_probs.extend(merged["pesto_interface_prob"].tolist())

        if merged["is_interface_contact"].nunique() < 2:
            n_single_class += 1
            continue
        per_chain_auc.append(roc_auc_score(merged["is_interface_contact"], merged["pesto_interface_prob"]))

    print(f"Primary set: {len(primary)} chains; {n_no_exp_prediction} without a valid exp prediction yet; "
          f"{n_single_class} skipped from per-chain AUC (single-class labels)")
    print(f"Pooled residues: {len(all_labels)}, {sum(all_labels)} interface ({sum(all_labels) / len(all_labels):.1%})")

    pooled_auc = roc_auc_score(all_labels, all_probs)
    print(f"\nPOOLED ROC-AUC (exp predictions vs. Phase 5 distance labels): {pooled_auc:.4f}")
    print(f"Per-chain AUC: n={len(per_chain_auc)}, mean={pd.Series(per_chain_auc).mean():.4f}, "
          f"median={pd.Series(per_chain_auc).median():.4f}")

    if pooled_auc < 0.6:
        print("\n*** STOP: pooled AUC is near chance. Investigate the channel index, key join, or "
              "input preparation before proceeding to Phase 7. ***")
    else:
        print("\nSanity check passed: pooled AUC is well above chance.")


if __name__ == "__main__":
    main()
