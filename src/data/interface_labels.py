"""Phase 5: ground-truth interface labels from the experimental complexes.

For every chain with a validated Phase 4 residue mapping, computes two
independent interface labels against the chain's protein partners in its
*biological assembly* (not the raw asymmetric unit -- see
``choose_assembly``/``build_assembly_model`` below): a primary heavy-atom
distance label and a robustness ΔSASA label, plus the relative solvent
accessibility (RSA) of the isolated chain. Labels are restricted to
exactly the residue set in the chain's Phase 4 residue map (same
auth_seq_id/ins_code keys), asserted to match exactly -- never a silently
assumed correspondence.

Runnable as: .venv/bin/python -m src.data.interface_labels
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from biotite.structure import AtomArray, sasa
from scipy.spatial import cKDTree

from src import config

logger = logging.getLogger(__name__)

# --- Paths -----------------------------------------------------------------

MAPPING_REPORT_PATH = config.INTERIM_DATA_DIR / "mapping_report.csv"
RESIDUE_MAPS_DIR = config.INTERIM_DATA_DIR / "residue_mappings"
INTERFACE_LABELS_DIR = config.INTERIM_DATA_DIR / "interface_labels"
LABELS_REPORT_PATH = config.INTERIM_DATA_DIR / "labels_report.csv"
ATTRITION_PATH = config.INTERIM_DATA_DIR / "phase5_attrition.csv"
AGREEMENT_PATH = config.INTERIM_DATA_DIR / "phase5_label_agreement.csv"

# --- Reference data ---------------------------------------------------------

# Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) "Maximum
# Allowed Solvent Accessibilities of Residues in Proteins", PLOS ONE --
# "Theoretical" column, one-letter code. Used as the RSA denominator.
# Modified residues are looked up via their parent one-letter code
# (gemmi.find_tabulated_residue), same convention as Phase 4.
MAX_SASA_ANGSTROM2 = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

# --- Exclusion reasons (logged distinctly, never silently dropped) ---------

CHAIN_EXCL_NOT_MAPPED = "phase4_mapping_not_available"
CHAIN_EXCL_READ_ERROR = "structure_read_error"
CHAIN_EXCL_CHAIN_NOT_IN_STRUCTURE = "chain_not_in_structure"
CHAIN_EXCL_NO_ASSEMBLY_FOR_CHAIN = "chain_not_in_any_assembly"
CHAIN_EXCL_SELF_COPY_NOT_FOUND = "assembly_self_copy_not_found"
CHAIN_EXCL_RESIDUE_KEY_MISMATCH = "residue_key_mismatch_with_phase4_map"
CHAIN_EXCL_NO_PROTEIN_PARTNER = "no_protein_partner_in_assembly"
CHAIN_EXCL_ZERO_INTERFACE_RESIDUES = "zero_interface_residues"

PARQUET_FIELDS = [
    "auth_seq_id",
    "auth_ins_code",
    "uniprot_resnum",
    "is_interface_contact",
    "min_partner_distance",
    "delta_sasa",
    "is_interface_sasa",
    "rsa",
]

LABELS_REPORT_FIELDS = [
    "cluster_id",
    "pdb_id",
    "chain_id",
    "uniprot_acc",
    "eligible_homolog",
    "eligible_exact_train",
    "eligible_none",
    "assembly_id",
    "n_partner_chains",
    "has_homomeric_partner",
    "n_observed",
    "interface_fraction_distance",
    "interface_fraction_sasa",
    "interface_fraction_distance_surface",
    "n_surface_residues",
    "jaccard",
    "n_partner_chains_asu_only",
    "would_exclude_if_asu",
    "status",
    "exclusion_reason",
]

ATTRITION_FIELDS = [
    "leakage_filter_mode",
    "n_eligible_representatives",
    "n_labeled",
    "n_excluded",
    "note",
]


class ChainExcluded(Exception):
    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


@dataclass(frozen=True)
class HeavyAtom:
    name: str
    element: str
    comp_id: str  # parent residue's 3-letter name, needed for ProtOr VdW radii
    pos: np.ndarray  # shape (3,)


@dataclass
class ProteinResidue:
    auth_seq_id: int
    ins_code: str
    comp_id: str
    entity_id: str
    atoms: list  # list[HeavyAtom]

    @property
    def coords(self) -> np.ndarray:
        return np.array([a.pos for a in self.atoms]) if self.atoms else np.empty((0, 3))


# --- Heavy-atom / alt-loc resolution (gemmi structured API) ---------------


def resolve_heavy_atoms(residue: gemmi.Residue) -> list[HeavyAtom]:
    """Heavy (non-hydrogen) atoms for one residue, resolving alternate
    locations deterministically: for atoms sharing the same atom *name*,
    the highest-occupancy conformer wins (ties broken by ascending
    altloc) -- same rule as Phase 4's residue-level altloc resolution,
    generalized here to the full atom set for geometry. An atom with no
    alternate (a single row, whatever its own occupancy) is always kept.
    Keeps each atom's own name/element -- required for biotite's ProtOr
    van-der-Waals radii, which are looked up per (residue name, atom
    name), not a single element-wide constant.
    """
    by_name: dict[str, list[gemmi.Atom]] = {}
    for atom in residue:
        if atom.element.is_hydrogen:
            continue
        by_name.setdefault(atom.name, []).append(atom)

    chosen = []
    for atoms in by_name.values():
        best = atoms[0] if len(atoms) == 1 else sorted(atoms, key=lambda a: (-a.occ, a.altloc))[0]
        chosen.append(best)
    return [
        HeavyAtom(
            name=a.name, element=a.element.name, comp_id=residue.name,
            pos=np.array([a.pos.x, a.pos.y, a.pos.z]),
        )
        for a in chosen
    ]


def one_letter_code(comp_id: str) -> str | None:
    info = gemmi.find_tabulated_residue(comp_id)
    if info is None or not info.is_amino_acid():
        return None
    return info.one_letter_code.upper()


# --- Assembly resolution ----------------------------------------------------


def _chain_subchains(chain: gemmi.Chain) -> set[str]:
    return {r.subchain for r in chain}


def _generator_covers_chain(gen: gemmi.Assembly.Gen, chain_id: str, subchains: set[str]) -> bool:
    return chain_id in gen.chains or bool(subchains & set(gen.subchains))


def choose_assembly(st: gemmi.Structure, chain_id: str) -> gemmi.Assembly | None:
    """Pick the assembly to use for chain_id: the (author-preferred, else
    first-listed) assembly whose generators cover this chain. Returns None
    if the file defines no assemblies at all -- the documented fallback is
    then to treat the asymmetric unit itself as "the assembly" (a common,
    legitimate case: many simple depositions have no separate
    pdbx_struct_assembly category because the asymmetric unit already *is*
    the biological unit). Raises ChainExcluded if assemblies exist but none
    of them actually cover this chain.
    """
    if not st.assemblies:
        return None
    subchains = _chain_subchains(st[0][chain_id])
    candidates = [
        asm for asm in st.assemblies
        if any(_generator_covers_chain(gen, chain_id, subchains) for gen in asm.generators)
    ]
    if not candidates:
        raise ChainExcluded(CHAIN_EXCL_NO_ASSEMBLY_FOR_CHAIN)
    order = {id(asm): i for i, asm in enumerate(st.assemblies)}
    candidates.sort(key=lambda a: (0 if a.author_determined else 1, order[id(a)]))
    return candidates[0]


def build_assembly_model(st: gemmi.Structure, assembly: gemmi.Assembly | None) -> gemmi.Model:
    if assembly is None:
        return st[0]
    return gemmi.make_assembly(assembly, st[0], gemmi.HowToNameCopiedChain.AddNumber, None)


def _matching_chains(model: gemmi.Model, chain_id: str) -> list[gemmi.Chain]:
    """Chains in model generated from chain_id by gemmi's AddNumber naming
    (original name + a 1-based copy index, e.g. "A" -> "A1", "A2", ...)."""
    pattern = re.compile(r"^" + re.escape(chain_id) + r"\d+$")
    return [c for c in model if c.name == chain_id or pattern.fullmatch(c.name)]


def _first_heavy_atom(chain: gemmi.Chain):
    """(seqid.num, icode, atom_name, pos) for the first heavy atom in chain, or None."""
    for res in chain:
        for atom in res:
            if not atom.element.is_hydrogen:
                return res.seqid.num, res.seqid.icode, atom.name, atom.pos
    return None


def _is_identity_copy(candidate: gemmi.Chain, asu_chain: gemmi.Chain) -> bool:
    """True if candidate's coordinates exactly match the deposited
    (untransformed) asymmetric-unit chain -- i.e. it was generated by the
    identity operator, not a rotation/translation. One matching heavy atom
    is decisive evidence either way (rigid-body transforms move every atom
    identically, or none at all)."""
    ref = _first_heavy_atom(asu_chain)
    if ref is None:
        return False
    seqnum, icode, atom_name, pos = ref
    for res in candidate:
        if res.seqid.num != seqnum or res.seqid.icode != icode:
            continue
        for atom in res:
            if atom.name == atom_name and not atom.element.is_hydrogen:
                return atom.pos.dist(pos) < 1e-2
    return False


def find_self_copy(model: gemmi.Model, asu_chain: gemmi.Chain, chain_id: str) -> gemmi.Chain:
    """Pick which chain in the (possibly symmetry-expanded) assembly model
    represents chain_id for labeling purposes.

    A chain's auth_seq_id/residue identity is unaffected by which rigid-
    body operator generated it (make_assembly preserves residue metadata,
    only atom coordinates move) -- so when chain_id appears via exactly one
    copy in the chosen assembly, that copy is unambiguously "the
    representative," identity-transformed or not (real depositions
    sometimes reach a chain only via a non-identity operator, e.g. when a
    biological trimer is built from one chain via the identity operator
    plus a second chain via a crystal-symmetry operator -- confirmed on
    real data, see PROGRESS.md).

    Only when chain_id appears via *multiple* copies in one assembly (a
    true multi-copy homomer within that assembly) is a tie-break needed --
    documented here: prefer the identity-transformed copy if one exists,
    else the lowest AddNumber copy index, both fully deterministic.
    """
    candidates = _matching_chains(model, chain_id)
    if not candidates:
        raise ChainExcluded(CHAIN_EXCL_SELF_COPY_NOT_FOUND)
    if len(candidates) == 1:
        return candidates[0]

    for candidate in candidates:
        if _is_identity_copy(candidate, asu_chain):
            return candidate
    return sorted(candidates, key=lambda c: c.name)[0]


def representative_entity_id(asu_chain: gemmi.Chain) -> str | None:
    for res in asu_chain:
        if res.label_seq is not None:
            return res.entity_id
    return None


# --- Residue/atom pools -----------------------------------------------------


def representative_residues(self_chain: gemmi.Chain, phase4_keys: set[tuple[int, str]]) -> dict:
    """Residue-keyed heavy-atom coordinates for the representative chain,
    restricted to exactly phase4_keys (PLAN.md Phase 5 item 6). Raises
    ChainExcluded if the geometry-derived residue set doesn't match
    phase4_keys exactly -- the assert this session's brief asked for."""
    found: dict[tuple[int, str], ProteinResidue] = {}
    for res in self_chain:
        if res.label_seq is None:
            continue  # not part of the polymer sequence (water/ligand sharing this auth chain)
        key = (res.seqid.num, res.seqid.icode.strip())
        atoms = resolve_heavy_atoms(res)
        found[key] = ProteinResidue(key[0], key[1], res.name, res.entity_id, atoms)

    missing_from_geometry = phase4_keys - set(found)
    extra_in_geometry = set(found) - phase4_keys
    if missing_from_geometry or extra_in_geometry:
        raise ChainExcluded(
            f"{CHAIN_EXCL_RESIDUE_KEY_MISMATCH}:missing={len(missing_from_geometry)},"
            f"extra={len(extra_in_geometry)}"
        )
    return found


def collect_partners(model: gemmi.Model, self_chain_name: str, rep_entity_id: str | None) -> dict:
    """Protein-polymer partner residues/atoms in the assembly, excluding
    the representative's own (self) copy. Excludes ligands, ions, water,
    and nucleic acids by checking each residue's own chemical-component
    data (gemmi.find_tabulated_residue(...).is_amino_acid()) -- not chain-
    or entity-level metadata, which can be absent/unreliable -- so partner
    chains need not have passed any candidate-selection filter themselves.
    """
    partner_chain_names: set[str] = set()
    homomeric = False
    partner_atoms: list[HeavyAtom] = []
    for chain in model:
        if chain.name == self_chain_name:
            continue
        chain_has_protein_residue = False
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            atoms = resolve_heavy_atoms(res)
            if not atoms:
                continue
            chain_has_protein_residue = True
            partner_atoms.extend(atoms)
            if rep_entity_id is not None and res.entity_id == rep_entity_id:
                homomeric = True
        if chain_has_protein_residue:
            partner_chain_names.add(chain.name)
    return {
        "partner_chain_names": partner_chain_names,
        "homomeric": homomeric,
        "partner_atoms": partner_atoms,
    }


def count_asu_only_partner_chains(st: gemmi.Structure, chain_id: str) -> int:
    """Sanity comparison (PLAN.md Phase 5 item 7): how many protein partner
    chains chain_id would have if contacts were computed in the raw
    asymmetric unit instead of the chosen biological assembly."""
    count = 0
    for chain in st[0]:
        if chain.name == chain_id:
            continue
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and info.is_amino_acid():
                count += 1
                break
    return count


# --- Distance label (KD-tree) -----------------------------------------------


def compute_distance_labels(residues: dict, partner_atoms: list) -> dict:
    """{key: (min_partner_distance, is_interface_contact)} via a KD-tree
    over pooled partner heavy atoms, so this scales to large assemblies."""
    if not partner_atoms:
        return {key: (float("inf"), False) for key in residues}
    partner_coords = np.array([a.pos for a in partner_atoms])
    tree = cKDTree(partner_coords)
    result = {}
    for key, residue in residues.items():
        if not residue.atoms:
            result[key] = (float("inf"), False)
            continue
        distances, _ = tree.query(residue.coords, k=1)
        min_distance = float(distances.min())
        result[key] = (min_distance, min_distance <= config.INTERFACE_DISTANCE_CUTOFF_ANGSTROM)
    return result


# --- SASA / RSA (biotite) ---------------------------------------------------


def _build_atom_array(atoms: list) -> AtomArray:
    """Build a biotite AtomArray with the (element, res_name, atom_name)
    metadata biotite's ProtOr van-der-Waals radii need (looked up per
    (residue name, atom name), not a bare element symbol) -- see PLAN.md
    Phase 5 for why a coordinates-only array silently/loudly breaks this.
    """
    arr = AtomArray(len(atoms))
    arr.coord = np.array([a.pos for a in atoms])
    arr.element = np.array([a.element for a in atoms])
    arr.atom_name = np.array([a.name for a in atoms])
    arr.res_name = np.array([a.comp_id for a in atoms])
    arr.res_id = np.zeros(len(atoms), dtype=int)
    arr.chain_id = np.array(["A"] * len(atoms))
    return arr


def compute_sasa_labels(residues: dict, ordered_keys: list, partner_atoms: list) -> dict:
    """{key: (delta_sasa, is_interface_sasa, rsa)} using biotite's
    Shrake-Rupley implementation (freesasa isn't installable in this
    container -- recorded risk, PLAN.md Sec 6), isolated vs. in the
    assembly context (representative + protein partners only, no
    ligands/waters -- consistent with the distance label's partner set).
    """
    rep_atom_lists = [residues[k].atoms for k in ordered_keys]
    rep_atom_counts = [len(a) for a in rep_atom_lists]
    rep_atoms_flat = [a for atoms in rep_atom_lists for a in atoms]

    isolated_array = _build_atom_array(rep_atoms_flat)
    isolated_sasa = sasa(isolated_array, vdw_radii="ProtOr")

    if partner_atoms:
        complex_array = _build_atom_array(rep_atoms_flat + list(partner_atoms))
        mask = np.zeros(len(complex_array), dtype=bool)
        mask[: len(rep_atoms_flat)] = True
        complex_sasa_full = sasa(complex_array, atom_filter=mask, vdw_radii="ProtOr")
        complex_sasa = complex_sasa_full[: len(rep_atoms_flat)]
    else:
        complex_sasa = isolated_sasa

    result = {}
    offset = 0
    for key, n_atoms in zip(ordered_keys, rep_atom_counts):
        residue_isolated = float(np.nansum(isolated_sasa[offset : offset + n_atoms]))
        residue_complex = float(np.nansum(complex_sasa[offset : offset + n_atoms]))
        delta = residue_isolated - residue_complex
        letter = one_letter_code(residues[key].comp_id)
        max_sasa = MAX_SASA_ANGSTROM2.get(letter) if letter else None
        rsa = (residue_isolated / max_sasa) if max_sasa else float("nan")
        result[key] = (delta, delta > config.SASA_BURIAL_CUTOFF_ANGSTROM2, rsa)
        offset += n_atoms
    return result


# --- Agreement (Jaccard + confusion matrix) --------------------------------


def jaccard_and_confusion(distance_flags: list[bool], sasa_flags: list[bool]) -> dict:
    tp = sum(1 for d, s in zip(distance_flags, sasa_flags) if d and s)
    fp = sum(1 for d, s in zip(distance_flags, sasa_flags) if d and not s)
    fn = sum(1 for d, s in zip(distance_flags, sasa_flags) if not d and s)
    tn = sum(1 for d, s in zip(distance_flags, sasa_flags) if not d and not s)
    union = tp + fp + fn
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "jaccard": (tp / union) if union else 0.0}


# --- One representative, start to finish ------------------------------------


def process_representative(pdb_cif_path: Path, pdb_id: str, chain_id: str, phase4_df: pd.DataFrame) -> tuple[list[dict] | None, dict]:
    base_report = {
        "assembly_id": "",
        "n_partner_chains": 0,
        "has_homomeric_partner": False,
        "n_observed": 0,
        "interface_fraction_distance": 0.0,
        "interface_fraction_sasa": 0.0,
        "interface_fraction_distance_surface": 0.0,
        "n_surface_residues": 0,
        "jaccard": 0.0,
        "n_partner_chains_asu_only": 0,
        "would_exclude_if_asu": False,
        "status": "excluded",
        "exclusion_reason": "",
    }
    phase4_keys = set(zip(phase4_df["auth_seq_id"].astype(int), phase4_df["auth_ins_code"].fillna("")))
    uniprot_by_key = {
        (int(row.auth_seq_id), row.auth_ins_code or ""): (
            None if pd.isna(row.uniprot_resnum) else int(row.uniprot_resnum)
        )
        for row in phase4_df.itertuples()
    }

    try:
        st = gemmi.read_structure(str(pdb_cif_path))
        st.setup_entities()
        if chain_id not in [c.name for c in st[0]]:
            raise ChainExcluded(CHAIN_EXCL_CHAIN_NOT_IN_STRUCTURE)
        asu_chain = st[0][chain_id]

        base_report["n_partner_chains_asu_only"] = count_asu_only_partner_chains(st, chain_id)
        base_report["would_exclude_if_asu"] = base_report["n_partner_chains_asu_only"] == 0

        assembly = choose_assembly(st, chain_id)
        base_report["assembly_id"] = assembly.name if assembly is not None else "asymmetric_unit"
        model = build_assembly_model(st, assembly)
        self_chain = asu_chain if assembly is None else find_self_copy(model, asu_chain, chain_id)

        residues = representative_residues(self_chain, phase4_keys)
        base_report["n_observed"] = len(residues)

        rep_entity_id = representative_entity_id(asu_chain)
        partners = collect_partners(model, self_chain.name, rep_entity_id)
        base_report["n_partner_chains"] = len(partners["partner_chain_names"])
        base_report["has_homomeric_partner"] = partners["homomeric"]
        if not partners["partner_chain_names"]:
            raise ChainExcluded(CHAIN_EXCL_NO_PROTEIN_PARTNER)

        distance_labels = compute_distance_labels(residues, partners["partner_atoms"])
        ordered_keys = sorted(residues)
        sasa_labels = compute_sasa_labels(residues, ordered_keys, partners["partner_atoms"])

        n_interface_distance = sum(1 for k in ordered_keys if distance_labels[k][1])
        if n_interface_distance == 0:
            raise ChainExcluded(CHAIN_EXCL_ZERO_INTERFACE_RESIDUES)

    except ChainExcluded as exc:
        base_report["exclusion_reason"] = exc.reason
        logger.warning("excluded %s_%s: %s", pdb_id, chain_id, exc.reason)
        return None, base_report
    except Exception as exc:  # noqa: BLE001 -- never crash the batch on one bad file
        base_report["exclusion_reason"] = f"{CHAIN_EXCL_READ_ERROR}:{exc}"
        logger.warning("excluded %s_%s: %s: %s", pdb_id, chain_id, CHAIN_EXCL_READ_ERROR, exc)
        return None, base_report

    rows = []
    distance_flags, sasa_flags = [], []
    n_surface = 0
    n_interface_distance_surface = 0
    for key in ordered_keys:
        min_dist, is_interface_dist = distance_labels[key]
        delta_sasa, is_interface_sasa, rsa = sasa_labels[key]
        distance_flags.append(is_interface_dist)
        sasa_flags.append(is_interface_sasa)
        is_surface = (not np.isnan(rsa)) and rsa >= config.SURFACE_RSA_THRESHOLD
        if is_surface:
            n_surface += 1
            if is_interface_dist:
                n_interface_distance_surface += 1
        rows.append(
            {
                "auth_seq_id": key[0],
                "auth_ins_code": key[1],
                "uniprot_resnum": uniprot_by_key.get(key),
                "is_interface_contact": is_interface_dist,
                "min_partner_distance": min_dist,
                "delta_sasa": delta_sasa,
                "is_interface_sasa": is_interface_sasa,
                "rsa": rsa,
            }
        )

    agreement = jaccard_and_confusion(distance_flags, sasa_flags)
    n = len(ordered_keys)
    base_report.update(
        {
            "interface_fraction_distance": sum(distance_flags) / n,
            "interface_fraction_sasa": sum(sasa_flags) / n,
            "interface_fraction_distance_surface": (
                n_interface_distance_surface / n_surface if n_surface else 0.0
            ),
            "n_surface_residues": n_surface,
            "jaccard": agreement["jaccard"],
            "status": "labeled",
        }
    )
    return rows, base_report


# --- Loading Phase 4 outputs / priority ordering ---------------------------


def _csv_bool(value: str) -> bool:
    return value == "True"


def load_mapped_representatives(path: Path = MAPPING_REPORT_PATH) -> list[dict]:
    """Load mapping_report.csv rows with status == 'mapped' (Phase 5's
    scope: all chains with a validated Phase 4 mapping)."""
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["status"] != "mapped":
                continue
            rows.append(
                {
                    "cluster_id": row["cluster_id"],
                    "pdb_id": row["pdb_id"],
                    "chain_id": row["chain_id"],
                    "uniprot_acc": row["uniprot_acc"],
                    "eligible_homolog": _csv_bool(row["eligible_homolog"]),
                    "eligible_exact_train": _csv_bool(row["eligible_exact_train"]),
                    "eligible_none": _csv_bool(row["eligible_none"]),
                }
            )
    return rows


def sort_by_priority(rows: list[dict]) -> list[dict]:
    def priority_rank(row: dict) -> tuple[int, str]:
        if row["eligible_homolog"]:
            rank = 0
        elif row["eligible_exact_train"]:
            rank = 1
        else:
            rank = 2
        return (rank, row["cluster_id"])

    return sorted(rows, key=priority_rank)


# --- Output writers / attrition --------------------------------------------


def write_parquet(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=PARQUET_FIELDS).to_parquet(path, index=False)


def write_csv(rows: list[dict], fields: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def compute_attrition(report_rows: list[dict]) -> list[dict]:
    attrition = []
    for mode in config.LEAKAGE_FILTER_MODES:
        col = f"eligible_{mode}"
        subset = [r for r in report_rows if r[col]]
        n_labeled = sum(1 for r in subset if r["status"] == "labeled")
        attrition.append(
            {
                "leakage_filter_mode": mode,
                "n_eligible_representatives": len(subset),
                "n_labeled": n_labeled,
                "n_excluded": len(subset) - n_labeled,
                "note": f"meets_target({config.MIN_TEST_CHAINS_TARGET})="
                f"{n_labeled >= config.MIN_TEST_CHAINS_TARGET}",
            }
        )
    return attrition


# --- Orchestration -----------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "interface_labels.log"

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
    rows = sort_by_priority(load_mapped_representatives())
    logger.info("Loaded %d Phase 4-mapped representatives", len(rows))

    report_rows = []
    all_distance_flags: list[bool] = []
    all_sasa_flags: list[bool] = []
    for i, rep in enumerate(rows):
        parquet_path = RESIDUE_MAPS_DIR / f"{rep['pdb_id']}_{rep['chain_id']}.parquet"
        pdb_cif_path = config.RAW_DATA_DIR / "pdb" / f"{rep['pdb_id']}_updated.cif.gz"
        base = {
            "cluster_id": rep["cluster_id"], "pdb_id": rep["pdb_id"], "chain_id": rep["chain_id"],
            "uniprot_acc": rep["uniprot_acc"], "eligible_homolog": rep["eligible_homolog"],
            "eligible_exact_train": rep["eligible_exact_train"], "eligible_none": rep["eligible_none"],
        }
        if not parquet_path.exists():
            base.update(
                {"status": "excluded", "exclusion_reason": CHAIN_EXCL_NOT_MAPPED, "assembly_id": "",
                 "n_partner_chains": 0, "has_homomeric_partner": False, "n_observed": 0,
                 "interface_fraction_distance": 0.0, "interface_fraction_sasa": 0.0,
                 "interface_fraction_distance_surface": 0.0, "n_surface_residues": 0, "jaccard": 0.0,
                 "n_partner_chains_asu_only": 0, "would_exclude_if_asu": False}
            )
            report_rows.append(base)
            continue

        phase4_df = pd.read_parquet(parquet_path)
        parquet_rows, result = process_representative(pdb_cif_path, rep["pdb_id"], rep["chain_id"], phase4_df)
        base.update(result)
        report_rows.append(base)

        if parquet_rows is not None:
            write_parquet(parquet_rows, INTERFACE_LABELS_DIR / f"{rep['pdb_id']}_{rep['chain_id']}.parquet")
            all_distance_flags.extend(r["is_interface_contact"] for r in parquet_rows)
            all_sasa_flags.extend(r["is_interface_sasa"] for r in parquet_rows)

        if (i + 1) % 250 == 0 or (i + 1) == len(rows):
            n_labeled_so_far = sum(1 for r in report_rows if r["status"] == "labeled")
            logger.info("Processed %d/%d representatives (%d labeled)", i + 1, len(rows), n_labeled_so_far)

    write_csv(report_rows, LABELS_REPORT_FIELDS, LABELS_REPORT_PATH)
    logger.info("Wrote labels report (%d rows) to %s", len(report_rows), LABELS_REPORT_PATH)

    attrition = compute_attrition(report_rows)
    write_csv(attrition, ATTRITION_FIELDS, ATTRITION_PATH)
    logger.info("Wrote per-leakage-mode attrition to %s", ATTRITION_PATH)
    for row in attrition:
        logger.info(
            "mode=%-14s n_eligible=%-6d n_labeled=%-6d %s",
            row["leakage_filter_mode"], row["n_eligible_representatives"], row["n_labeled"], row["note"],
        )

    pooled = jaccard_and_confusion(all_distance_flags, all_sasa_flags)
    write_csv([{"leakage_filter_mode": "none (pooled, all labeled chains)", **pooled}],
              ["leakage_filter_mode", "tp", "fp", "fn", "tn", "jaccard"], AGREEMENT_PATH)
    logger.info("Pooled distance-vs-SASA agreement: %s", pooled)

    exclusion_counts: dict[str, int] = {}
    for r in report_rows:
        if r["exclusion_reason"]:
            key = r["exclusion_reason"].split(":", 1)[0]
            exclusion_counts[key] = exclusion_counts.get(key, 0) + 1
    logger.info("Exclusion-reason breakdown: %s", exclusion_counts)

    n_would_exclude_asu = sum(1 for r in report_rows if r.get("would_exclude_if_asu"))
    logger.info(
        "Asymmetric-unit sanity comparison: %d/%d representatives would have zero protein "
        "partners if the assembly were the ASU instead of the chosen biological assembly",
        n_would_exclude_asu, len(report_rows),
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Phase 5: ground-truth interface labels")
    parser.parse_args(argv)
    run()


if __name__ == "__main__":
    main()
