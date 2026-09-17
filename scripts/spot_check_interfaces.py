"""Hand-verification helper for Phase 5's interface labels.

Picks N random *labeled* chains from the primary (eligible_homolog) set,
seeded for reproducibility, and for each prints a ChimeraX command block
that opens the chosen biological assembly directly from RCSB (independent
of our own local assembly-generation code, for a genuinely independent
check), colors the distance-labeled interface residues on the
representative chain, and shows the partner chains -- so a human can see
whether the labeled patch actually faces the partner.

This is a manual spot-check tool, not part of the Phase 5 pipeline itself
-- it reads Phase 5's already-written outputs (labels_report.csv, the
per-chain parquet files) and does not recompute anything.

Usage: .venv/bin/python scripts/spot_check_interfaces.py --n 3 --seed 0
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import config
from src.data.interface_labels import INTERFACE_LABELS_DIR, LABELS_REPORT_PATH


def pick_chains(n: int, seed: int) -> pd.DataFrame:
    report = pd.read_csv(LABELS_REPORT_PATH)
    primary_labeled = report[(report["status"] == "labeled") & (report["eligible_homolog"])]
    if len(primary_labeled) < n:
        raise SystemExit(f"Only {len(primary_labeled)} labeled primary-set chains, need {n}")
    return primary_labeled.sample(n=n, random_state=seed).reset_index(drop=True)


def interface_residues(pdb_id: str, chain_id: str) -> pd.DataFrame:
    path = INTERFACE_LABELS_DIR / f"{pdb_id}_{chain_id}.parquet"
    df = pd.read_parquet(path)
    return df[df["is_interface_contact"]].sort_values("auth_seq_id")


def print_table(chains: pd.DataFrame) -> None:
    cols = [
        "pdb_id", "chain_id", "uniprot_acc", "assembly_id", "n_partner_chains",
        "has_homomeric_partner", "n_observed", "interface_fraction_distance", "jaccard",
    ]
    print(chains[cols].to_string(index=False))


def print_chimerax_block(pdb_id: str, chain_id: str, assembly_id: str, residues: pd.DataFrame) -> None:
    res_spec = ",".join(str(int(r)) for r in residues["auth_seq_id"])
    open_cmd = (
        f"open {pdb_id} assembly {assembly_id}" if assembly_id != "asymmetric_unit" else f"open {pdb_id}"
    )
    print(f"# --- {pdb_id}_{chain_id} (assembly {assembly_id}, {len(residues)} interface residues) ---")
    print(open_cmd)
    print("hide atoms")
    print("show cartoons")
    print("color #1 gray")
    print(f"color #1/{chain_id} cyan")
    if res_spec:
        print(f"color #1/{chain_id}:{res_spec} orange")
        print(f"show #1/{chain_id}:{res_spec} atoms")
    print(f"view #1/{chain_id}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=3, help="number of chains to spot-check")
    parser.add_argument("--seed", type=int, default=0, help="random seed for reproducible picks")
    args = parser.parse_args()

    chains = pick_chains(args.n, args.seed)
    print(f"Spot-check: {args.n} primary-set labeled chains, seed={args.seed}\n")
    print_table(chains)
    print()
    print("# ChimeraX command block -- paste into the ChimeraX command line.")
    print("# Each block fetches the chosen biological assembly directly from RCSB")
    print("# (independent of this project's own local assembly-generation code) and")
    print("# highlights the distance-labeled interface residues on the representative chain.")
    print()
    for _, chain in chains.iterrows():
        residues = interface_residues(chain["pdb_id"], chain["chain_id"])
        print_chimerax_block(chain["pdb_id"], chain["chain_id"], str(chain["assembly_id"]), residues)


if __name__ == "__main__":
    main()
