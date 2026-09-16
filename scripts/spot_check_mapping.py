"""Hand-verification helper for Phase 4's residue mapping.

Picks N random chains from the primary (eligible_homolog) mapped set,
seeded for reproducibility, and for each picks 3 random *kept* mapped
residues. Prints a plain table plus a ChimeraX command block that opens
both structures and selects those residues, so a human can visually
confirm the experimental residue and the AlphaFold residue at the
"same" (author vs. UniProt-mapped) position really do correspond.

This is a manual spot-check tool, not part of the Phase 4 pipeline itself
-- it reads Phase 4's already-written outputs (mapping_report.csv, the
per-chain parquet files) and does not recompute anything.

Usage: .venv/bin/python scripts/spot_check_mapping.py --n 3 --seed 0
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import config
from src.data.fetch_structures import ALPHAFOLD_DIR, PDB_UPDATED_DIR

MAPPING_REPORT_PATH = config.INTERIM_DATA_DIR / "mapping_report.csv"
RESIDUE_MAPS_DIR = config.INTERIM_DATA_DIR / "residue_mappings"


def pick_chains(n: int, seed: int) -> pd.DataFrame:
    report = pd.read_csv(MAPPING_REPORT_PATH)
    primary_mapped = report[(report["status"] == "mapped") & (report["eligible_homolog"])]
    if len(primary_mapped) < n:
        raise SystemExit(
            f"Only {len(primary_mapped)} mapped primary-set chains available, need {n}"
        )
    return primary_mapped.sample(n=n, random_state=seed).reset_index(drop=True)


def pick_residues(pdb_id: str, chain_id: str, uniprot_acc: str, seed: int, k: int = 3) -> pd.DataFrame:
    path = RESIDUE_MAPS_DIR / f"{pdb_id}_{chain_id}.parquet"
    df = pd.read_parquet(path)
    kept = df[df["unmapped_reason"] == ""]
    if len(kept) < k:
        raise SystemExit(f"{pdb_id}_{chain_id} has only {len(kept)} kept residues, need {k}")
    return kept.sample(n=k, random_state=seed).sort_values("auth_seq_id").reset_index(drop=True)


def format_auth_residue(auth_seq_id: int, ins_code: str) -> str:
    return f"{auth_seq_id}{ins_code}" if ins_code else str(auth_seq_id)


def chimerax_residue_spec(chain_id: str, auth_seq_id: int, ins_code: str) -> str:
    number = format_auth_residue(auth_seq_id, ins_code)
    return f"/{chain_id}:{number}"


def print_table(rows: list[dict]) -> None:
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))


def print_chimerax_block(pdb_id: str, chain_id: str, uniprot_acc: str, residues: pd.DataFrame) -> None:
    pdb_path = (PDB_UPDATED_DIR / f"{pdb_id.upper()}_updated.cif.gz").resolve()
    af_path = (ALPHAFOLD_DIR / f"{uniprot_acc}.cif.gz").resolve()

    exp_specs = ",".join(
        format_auth_residue(row["auth_seq_id"], row["auth_ins_code"]) for _, row in residues.iterrows()
    )
    af_specs = ",".join(str(int(row["alphafold_resnum"])) for _, row in residues.iterrows())

    print(f"# --- {pdb_id}_{chain_id} vs AlphaFold {uniprot_acc} ---")
    print(f"open {pdb_path}")
    print(f"open {af_path}")
    print(f"select #1/{chain_id}:{exp_specs}")
    print(f"select add #2/A:{af_specs}")
    print("color sel orange")
    print("show sel atoms")
    print("view sel")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=3, help="number of chains to spot-check")
    parser.add_argument("--seed", type=int, default=0, help="random seed for reproducible picks")
    parser.add_argument("--k", type=int, default=3, help="number of residues per chain")
    args = parser.parse_args()

    chains = pick_chains(args.n, args.seed)

    all_rows = []
    chimerax_blocks = []
    for _, chain in chains.iterrows():
        pdb_id, chain_id, uniprot_acc = chain["pdb_id"], chain["chain_id"], chain["uniprot_acc"]
        residues = pick_residues(pdb_id, chain_id, uniprot_acc, seed=args.seed, k=args.k)
        for _, r in residues.iterrows():
            all_rows.append(
                {
                    "pdb_id": pdb_id,
                    "chain": chain_id,
                    "auth_residue": format_auth_residue(r["auth_seq_id"], r["auth_ins_code"]),
                    "exp_residue_name": r["residue_name"],
                    "uniprot_position": int(r["uniprot_resnum"]),
                    "af_residue_name": r["alphafold_residue_name"],
                    "match": r["match_flag"],
                    "method": r["method"],
                }
            )
        chimerax_blocks.append((pdb_id, chain_id, uniprot_acc, residues))

    print(f"Spot-check: {args.n} chains x {args.k} residues, seed={args.seed}\n")
    print_table(all_rows)
    print()
    print("# ChimeraX command block -- paste into the ChimeraX command line.")
    print("# #1 = experimental structure, #2 = AlphaFold model, per chain checked below.")
    print()
    for pdb_id, chain_id, uniprot_acc, residues in chimerax_blocks:
        print_chimerax_block(pdb_id, chain_id, uniprot_acc, residues)


if __name__ == "__main__":
    main()
