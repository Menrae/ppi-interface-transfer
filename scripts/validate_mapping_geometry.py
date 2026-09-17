"""Geometric validation of Phase 4's residue mapping, for the primary set.

Phase 4 validated its residue mapping by *sequence* identity only (SIFTS or
the fallback alignment, gated at >=90% identity -- PLAN.md Phase 4). That
catches numbering/frame errors, but not a case where the mapping is
numerically self-consistent yet the two structures don't actually
correspond in 3D (e.g. a genuinely different but similar-sequence protein,
or a domain that's been swapped/misassigned) -- high sequence identity
does not guarantee low geometric distance.

This script closes that gap: for every primary-set Phase-4-mapped chain,
superpose the experimental chain's Calpha atoms onto the AlphaFold model's
Calpha atoms at the same UniProt positions (Kabsch/SVD, one rigid-body fit
per chain), then look at the post-superposition per-residue Calpha
distance. A chain passes if the fit is good globally; it's flagged if a
long contiguous run of residues stays far apart despite passing Phase 4's
identity gate and (usually) fitting well everywhere else -- exactly the
"high identity, wrong geometry" case sequence-only validation can't see.

Not part of the Phase 4 pipeline itself -- reads Phase 4's already-written
outputs and the cached structures, and does not change any Phase 4 file.

Usage: .venv/bin/python scripts/validate_mapping_geometry.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import config

MAPPING_REPORT_PATH = config.INTERIM_DATA_DIR / "mapping_report.csv"
RESIDUE_MAPS_DIR = config.INTERIM_DATA_DIR / "residue_mappings"
PDB_UPDATED_DIR = config.RAW_DATA_DIR / "pdb"
ALPHAFOLD_DIR = config.RAW_DATA_DIR / "alphafold"
OUTPUT_PATH = config.INTERIM_DATA_DIR / "phase4_geometric_validation.csv"

# A run of at least this many *consecutive* mapped residues (consecutive in
# the sorted-by-UniProt-position list, not necessarily consecutive
# auth_seq_id) all beyond FAR_DISTANCE_ANGSTROM after superposition is
# flagged as a likely mapping/register problem rather than ordinary local
# flexibility (which is usually isolated single residues or short loops).
FAR_DISTANCE_ANGSTROM = 10.0
MIN_FLAGGED_RUN_LENGTH = 5


def get_ca_coords(cif_path: Path, chain_id: str | None, key_type: str, keys: list) -> dict:
    """key_type: "auth" -> keys are (auth_seq_id, ins_code); "uniprot" -> keys are uniprot_resnum ints."""
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    if chain_id is None:
        chain_id = st[0][0].name if len(list(st[0])) == 1 else None
        if chain_id is None:
            raise ValueError(f"expected a single chain in {cif_path}")
    chain = st[0][chain_id]

    coords = {}
    for res in chain:
        if res.label_seq is None:
            continue
        key = (res.seqid.num, res.seqid.icode.strip()) if key_type == "auth" else res.seqid.num
        for atom in res:
            if atom.name == "CA" and not atom.element.is_hydrogen:
                coords[key] = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                break
    return {k: coords[k] for k in keys if k in coords}


def kabsch_superpose(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rigid-body superpose `mobile` onto `target` (both (n,3)); returns
    mobile's coordinates after the best-fit rotation+translation."""
    mobile_c = mobile - mobile.mean(axis=0)
    target_c = target - target.mean(axis=0)
    H = mobile_c.T @ target_c
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return (R @ mobile_c.T).T + target.mean(axis=0)


def longest_far_run(distances: np.ndarray, threshold: float) -> int:
    longest = current = 0
    for d in distances:
        if d > threshold:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def validate_chain(pdb_id: str, chain_id: str, uniprot_acc: str, af_model_path: Path) -> dict:
    mapping = pd.read_parquet(RESIDUE_MAPS_DIR / f"{pdb_id}_{chain_id}.parquet")
    mapped = mapping[mapping["uniprot_resnum"].notna()].sort_values("uniprot_resnum")
    if len(mapped) < 5:
        return {"status": "skipped", "reason": "too_few_mapped_residues", "n_mapped": len(mapped)}

    auth_keys = list(zip(mapped["auth_seq_id"].astype(int), mapped["auth_ins_code"].fillna("")))
    uniprot_keys = list(mapped["uniprot_resnum"].astype(int))

    exp_coords = get_ca_coords(PDB_UPDATED_DIR / f"{pdb_id}_updated.cif.gz", chain_id, "auth", auth_keys)
    af_coords = get_ca_coords(af_model_path, None, "uniprot", uniprot_keys)

    common = [
        (a_key, u_key) for a_key, u_key in zip(auth_keys, uniprot_keys)
        if a_key in exp_coords and u_key in af_coords
    ]
    if len(common) < 5:
        return {"status": "skipped", "reason": "too_few_ca_atoms_found", "n_mapped": len(mapped)}

    exp_xyz = np.array([exp_coords[a] for a, _ in common])
    af_xyz = np.array([af_coords[u] for _, u in common])

    exp_superposed = kabsch_superpose(exp_xyz, af_xyz)
    distances = np.linalg.norm(exp_superposed - af_xyz, axis=1)

    return {
        "status": "validated",
        "n_mapped": len(common),
        "median_distance": float(np.median(distances)),
        "mean_distance": float(np.mean(distances)),
        "max_distance": float(np.max(distances)),
        "longest_far_run": longest_far_run(distances, FAR_DISTANCE_ANGSTROM),
        "frac_far": float(np.mean(distances > FAR_DISTANCE_ANGSTROM)),
    }


def main() -> None:
    mapping_report = pd.read_csv(MAPPING_REPORT_PATH)
    primary = mapping_report[(mapping_report["status"] == "mapped") & (mapping_report["eligible_homolog"])]
    print(f"Validating {len(primary)} primary-set Phase-4-mapped chains...")

    rows = []
    for i, rep in enumerate(primary.itertuples()):
        af_path = ALPHAFOLD_DIR / f"{rep.uniprot_acc}.cif.gz"
        try:
            result = validate_chain(rep.pdb_id, rep.chain_id, rep.uniprot_acc, af_path)
        except Exception as exc:  # noqa: BLE001
            result = {"status": "error", "reason": str(exc)}
        result.update({"pdb_id": rep.pdb_id, "chain_id": rep.chain_id, "uniprot_acc": rep.uniprot_acc})
        rows.append(result)
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(primary)}")

    fields = ["pdb_id", "chain_id", "uniprot_acc", "status", "reason", "n_mapped",
              "median_distance", "mean_distance", "max_distance", "longest_far_run", "frac_far"]
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)

    df = pd.DataFrame(rows)
    validated = df[df["status"] == "validated"]
    flagged = validated[validated["longest_far_run"] >= MIN_FLAGGED_RUN_LENGTH]

    print(f"\nValidated: {len(validated)}/{len(df)} chains ({df['status'].value_counts().to_dict()})")
    print(f"Median distance distribution: {validated['median_distance'].describe().to_dict()}")
    print(f"\nFlagged ({MIN_FLAGGED_RUN_LENGTH}+ consecutive residues >{FAR_DISTANCE_ANGSTROM}A after superposition): "
          f"{len(flagged)}/{len(validated)}")
    if len(flagged):
        print(flagged[["pdb_id", "chain_id", "uniprot_acc", "median_distance", "longest_far_run", "frac_far"]]
              .to_string(index=False))
    print(f"\nWrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
