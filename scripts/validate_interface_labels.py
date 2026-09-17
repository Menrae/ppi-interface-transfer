"""Programmatic (no-ChimeraX) verification of Phase 5's interface labels.

Attempts PDBe's PISA/interface-annotation API first, for an independent
external check against a seeded sample of primary-set chains. As of this
session, PDBe's public API didn't serve any of the interface-annotation
endpoint patterns tried (see PROGRESS.md for what was attempted) -- this
script therefore falls back to the two purely-geometric checks specified
as the fallback: (1) most labeled interface residues have another
interface residue nearby (spatial contiguity -- real interfaces are
patches, not scattered single residues), and (2) the interface patch's
centroid sits closer to the partner chain than the centroid of
comparably-exposed non-interface surface residues does (the patch should
actually face the partner, not just be "some exposed residues").

Not part of the Phase 5 pipeline -- reads Phase 5's already-written
outputs and the cached assemblies, and does not change any Phase 5 file.

Usage: .venv/bin/python scripts/validate_interface_labels.py --n 20 --seed 0
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import config
from src.data.interface_labels import (
    choose_assembly, build_assembly_model, find_self_copy, representative_entity_id,
    collect_partners, ChainExcluded,
)

LABELS_REPORT_PATH = config.INTERIM_DATA_DIR / "labels_report.csv"
INTERFACE_LABELS_DIR = config.INTERIM_DATA_DIR / "interface_labels"
PDB_UPDATED_DIR = config.RAW_DATA_DIR / "pdb"
OUTPUT_PATH = config.INTERIM_DATA_DIR / "phase5_spatial_verification.csv"

CONTACT_NEIGHBOR_DISTANCE_ANGSTROM = 8.0


def try_pdbe_pisa(pdb_id: str) -> dict | None:
    """Best-effort PDBe interface-annotation fetch; returns None if
    unavailable (never raises, never blocks the geometric fallback)."""
    import requests

    candidates = [
        f"https://www.ebi.ac.uk/pdbe/api/pisa/interfacelist/{pdb_id.lower()}/1",
        f"https://www.ebi.ac.uk/pdbe/api/pisa/interfacelist/{pdb_id.lower()}",
        f"https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pisa/interface_residues/{pdb_id.lower()}",
    ]
    for url in candidates:
        try:
            response = requests.get(url, timeout=10)
        except requests.RequestException:
            continue
        if response.status_code == 200:
            return response.json()
    return None


def get_ca_by_key(chain: gemmi.Chain) -> dict:
    coords = {}
    for res in chain:
        if res.label_seq is None:
            continue
        key = (res.seqid.num, res.seqid.icode.strip())
        for atom in res:
            if atom.name == "CA" and not atom.element.is_hydrogen:
                coords[key] = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                break
    return coords


def validate_chain(pdb_id: str, chain_id: str) -> dict:
    labels = pd.read_parquet(INTERFACE_LABELS_DIR / f"{pdb_id}_{chain_id}.parquet")

    st = gemmi.read_structure(str(PDB_UPDATED_DIR / f"{pdb_id}_updated.cif.gz"))
    st.setup_entities()
    asu_chain = st[0][chain_id]

    assembly = choose_assembly(st, chain_id)
    model = build_assembly_model(st, assembly)
    self_chain = asu_chain if assembly is None else find_self_copy(model, asu_chain, chain_id)
    rep_entity_id = representative_entity_id(asu_chain)
    partners = collect_partners(model, self_chain.name, rep_entity_id)
    if not partners["partner_atoms"]:
        return {"status": "skipped", "reason": "no_partner_atoms_found"}
    partner_centroid = np.mean([a.pos for a in partners["partner_atoms"]], axis=0)

    ca_by_key = get_ca_by_key(asu_chain)

    interface_keys = [(int(r.auth_seq_id), r.auth_ins_code or "") for r in labels[labels["is_interface_contact"]].itertuples()]
    surface_keys = [
        (int(r.auth_seq_id), r.auth_ins_code or "")
        for r in labels[(~labels["is_interface_contact"]) & (labels["rsa"] >= config.SURFACE_RSA_THRESHOLD)].itertuples()
    ]

    interface_ca = np.array([ca_by_key[k] for k in interface_keys if k in ca_by_key])
    surface_ca = np.array([ca_by_key[k] for k in surface_keys if k in ca_by_key])

    if len(interface_ca) == 0:
        return {"status": "skipped", "reason": "no_interface_residues_with_ca"}

    if len(interface_ca) > 1:
        tree = cKDTree(interface_ca)
        dists, _ = tree.query(interface_ca, k=2)
        frac_with_neighbor = float(np.mean(dists[:, 1] <= CONTACT_NEIGHBOR_DISTANCE_ANGSTROM))
    else:
        frac_with_neighbor = None  # single-residue interface; contiguity undefined

    interface_centroid = interface_ca.mean(axis=0)
    dist_interface_to_partner = float(np.linalg.norm(interface_centroid - partner_centroid))

    if len(surface_ca) > 0:
        surface_centroid = surface_ca.mean(axis=0)
        dist_surface_to_partner = float(np.linalg.norm(surface_centroid - partner_centroid))
        centroid_check_passed = dist_interface_to_partner < dist_surface_to_partner
    else:
        dist_surface_to_partner = None
        centroid_check_passed = None

    return {
        "status": "validated",
        "n_interface_residues": len(interface_ca),
        "n_surface_noninterface_residues": len(surface_ca),
        "frac_interface_with_neighbor_within_8A": frac_with_neighbor,
        "dist_interface_centroid_to_partner": dist_interface_to_partner,
        "dist_surface_centroid_to_partner": dist_surface_to_partner,
        "centroid_check_passed": centroid_check_passed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    report = pd.read_csv(LABELS_REPORT_PATH)
    primary_labeled = report[(report["status"] == "labeled") & (report["eligible_homolog"])]
    sample = primary_labeled.sample(n=args.n, random_state=args.seed)

    print(f"Attempting PDBe PISA/interface API for {args.n} sample entries...")
    pisa_hits = 0
    for pdb_id in sample["pdb_id"].unique():
        if try_pdbe_pisa(pdb_id) is not None:
            pisa_hits += 1
    print(f"PDBe interface-annotation API reachable for {pisa_hits}/{sample['pdb_id'].nunique()} entries tried "
          f"(0 means none of the tried endpoint patterns returned data this session -- falling back to the "
          f"geometric checks below for all {args.n} chains).")

    rows = []
    for rep in sample.itertuples():
        try:
            result = validate_chain(rep.pdb_id, rep.chain_id)
        except ChainExcluded as exc:
            result = {"status": "skipped", "reason": exc.reason}
        except Exception as exc:  # noqa: BLE001
            result = {"status": "error", "reason": str(exc)}
        result.update({"pdb_id": rep.pdb_id, "chain_id": rep.chain_id, "uniprot_acc": rep.uniprot_acc})
        rows.append(result)

    fields = ["pdb_id", "chain_id", "uniprot_acc", "status", "reason", "n_interface_residues",
              "n_surface_noninterface_residues", "frac_interface_with_neighbor_within_8A",
              "dist_interface_centroid_to_partner", "dist_surface_centroid_to_partner", "centroid_check_passed"]
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)

    df = pd.DataFrame(rows)
    validated = df[df["status"] == "validated"]
    print(f"\nValidated {len(validated)}/{len(df)} sampled chains")
    print(validated[["pdb_id", "chain_id", "n_interface_residues", "frac_interface_with_neighbor_within_8A",
                      "dist_interface_centroid_to_partner", "dist_surface_centroid_to_partner", "centroid_check_passed"]]
          .to_string(index=False))
    contiguity_vals = validated["frac_interface_with_neighbor_within_8A"].dropna()
    print(f"\nContiguity: mean fraction of interface residues with another interface residue "
          f"within {CONTACT_NEIGHBOR_DISTANCE_ANGSTROM}A = {contiguity_vals.mean():.3f} "
          f"(n={len(contiguity_vals)} chains with >1 interface residue)")
    centroid_vals = validated["centroid_check_passed"].dropna()
    print(f"Centroid check passed: {centroid_vals.sum()}/{len(centroid_vals)} chains "
          f"(interface patch centroid closer to partner than non-interface surface centroid)")
    print(f"\nWrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
