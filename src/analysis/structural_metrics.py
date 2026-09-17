"""Phase 8: per-residue structural/confidence features for error analysis.

Four features per mapped residue, feeding src.analysis.error_analysis:

1. pLDDT -- read from the AlphaFold mmCIF's ``_atom_site.B_iso_or_equiv``
   column (gemmi's ``atom.b_iso``), confirmed (not assumed) by direct
   inspection: values fall in AlphaFold's documented 0-100 pLDDT range and
   are identical across all atoms of one residue.
2. Local RMSD -- Calpha RMSD over the mapped residues within
   ``config.LOCAL_RMSD_RADIUS_ANGSTROM`` of a given residue *in the
   experimental structure*, with a separate Kabsch superposition fit per
   neighborhood (not one whole-chain global fit). See PLAN.md Phase 8
   analysis plan for why local, not global, superposition is used.
   ``global_ca_distance`` (one global fit per chain, reusing Phase 4's
   `scripts/validate_mapping_geometry.py` method) is kept as a separate,
   complementary measure.
3. RSA -- reused directly from Phase 5's per-residue output, not
   recomputed here (see src.analysis.error_analysis, which reads it off
   the same interface_labels parquet Phase 7 already joins).
4. Secondary structure -- ``mkdssp``/DSSP is confirmed not installed in
   this container (no root access); falls back to
   ``biotite.structure.annotate_sse`` (P-SEA, 3-state: helix/strand/coil)
   on the experimental chain's Calpha trace.

The numeric core (Kabsch superposition, local-RMSD-per-residue, pLDDT-band
assignment) takes plain numpy arrays so it can be unit-tested against
synthetic coordinates without any real structure files. Thin gemmi/biotite
I/O helpers sit alongside it.
"""

from __future__ import annotations

import logging
from pathlib import Path

import gemmi
import numpy as np

try:
    import biotite.structure as bstruc
except ImportError:  # pragma: no cover - biotite is a hard requirement elsewhere
    bstruc = None

from src import config

logger = logging.getLogger(__name__)

SSE_LABELS = ("a", "b", "c")  # helix, strand, coil (biotite's P-SEA convention)


# --- Pure numeric core (unit-tested with synthetic coordinates) -----------


def kabsch_superpose(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rigid-body (rotation + translation) superposition of `mobile` (n,3)
    onto `target` (n,3); returns mobile's coordinates after the best fit."""
    mobile_c = mobile - mobile.mean(axis=0)
    target_c = target - target.mean(axis=0)
    h = mobile_c.T @ target_c
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    diag = np.diag([1.0, 1.0, d])
    r = vt.T @ diag @ u.T
    return (r @ mobile_c.T).T + target.mean(axis=0)


def global_ca_distance_per_residue(exp_coords: np.ndarray, af_coords: np.ndarray) -> np.ndarray:
    """One whole-chain Kabsch fit (exp onto af); returns per-residue
    post-superposition distance, same method as
    scripts/validate_mapping_geometry.py (reused, not re-derived)."""
    exp_superposed = kabsch_superpose(exp_coords, af_coords)
    return np.linalg.norm(exp_superposed - af_coords, axis=1)


def local_rmsd_per_residue(
    exp_coords: np.ndarray,
    af_coords: np.ndarray,
    radius: float = config.LOCAL_RMSD_RADIUS_ANGSTROM,
    min_neighbors: int = config.LOCAL_RMSD_MIN_NEIGHBORS,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """For each residue i (row of exp_coords/af_coords, same order/keys in
    both), finds its neighborhood (residues within `radius` of i's Calpha
    in the EXPERIMENTAL structure, always including i itself), superposes
    that neighborhood's exp coords onto its af coords with its own Kabsch
    fit, and returns the RMSD of that local fit.

    Returns (local_rmsd, n_neighbors, reasons) -- local_rmsd/n_neighbors are
    float/int arrays of length n (NaN/0 where undefined); reasons[i] is ""
    when defined, else a short string (e.g. "too_few_neighbors").
    """
    n = len(exp_coords)
    local_rmsd = np.full(n, np.nan)
    n_neighbors = np.zeros(n, dtype=int)
    reasons = [""] * n
    if n == 0:
        return local_rmsd, n_neighbors, reasons

    # Pairwise distances in the experimental structure define neighborhoods.
    diffs = exp_coords[:, None, :] - exp_coords[None, :, :]
    dist_matrix = np.linalg.norm(diffs, axis=-1)

    for i in range(n):
        neighbor_idx = np.where(dist_matrix[i] <= radius)[0]
        n_neighbors[i] = len(neighbor_idx)
        if len(neighbor_idx) < min_neighbors:
            reasons[i] = "too_few_neighbors"
            continue
        local_exp = exp_coords[neighbor_idx]
        local_af = af_coords[neighbor_idx]
        superposed = kabsch_superpose(local_exp, local_af)
        sq_dev = np.sum((superposed - local_af) ** 2, axis=1)
        local_rmsd[i] = float(np.sqrt(np.mean(sq_dev)))

    return local_rmsd, n_neighbors, reasons


def plddt_band_labels(bands: tuple[float, ...] = config.PLDDT_BANDS) -> list[str]:
    """Ordered band labels for `bands` (ascending edges), e.g. (50,70,90) ->
    ["<50", "50-70", "70-90", ">=90"]."""
    edges = list(bands)
    labels = [f"<{edges[0]:g}"]
    labels += [f"{lo:g}-{hi:g}" for lo, hi in zip(edges[:-1], edges[1:])]
    labels += [f">={edges[-1]:g}"]
    return labels


def assign_plddt_band(plddt: np.ndarray, bands: tuple[float, ...] = config.PLDDT_BANDS) -> np.ndarray:
    """Assigns each pLDDT value to a band string, following AlphaFold's own
    convention that band edges belong to the upper band (e.g. exactly 70.0
    is "70-90", not "50-70"). `bands` must be sorted ascending."""
    plddt = np.asarray(plddt, dtype=float)
    labels = plddt_band_labels(bands)
    band_idx = np.searchsorted(np.asarray(bands, dtype=float), plddt, side="right")
    return np.array(labels, dtype=object)[band_idx]


# --- Thin I/O helpers (gemmi/biotite) --------------------------------------


def read_ca_coords_by_auth(cif_path: Path, chain_id: str, keys: list[tuple[int, str]]) -> dict:
    """keys: list of (auth_seq_id, auth_ins_code). Returns {key: (3,) ndarray}."""
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    chain = st[0][chain_id]
    coords = {}
    for res in chain:
        key = (res.seqid.num, res.seqid.icode.strip())
        for atom in res:
            if atom.name == "CA" and not atom.element.is_hydrogen:
                coords[key] = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                break
    return {k: coords[k] for k in keys if k in coords}


def read_ca_coords_by_uniprot(cif_path: Path, keys: list[int]) -> dict:
    """AlphaFold models: single chain, residue numbering is UniProt-native.
    Returns {uniprot_resnum: (3,) ndarray}."""
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    chain = st[0][0]
    coords = {}
    for res in chain:
        for atom in res:
            if atom.name == "CA" and not atom.element.is_hydrogen:
                coords[res.seqid.num] = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                break
    return {k: coords[k] for k in keys if k in coords}


def read_plddt_by_uniprot(cif_path: Path, keys: list[int]) -> dict:
    """Returns {uniprot_resnum: plddt}. Logs (does not raise) if a
    residue's atoms disagree on b_iso beyond floating tolerance -- expected
    to never happen per this session's inspection, but checked rather than
    assumed."""
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    chain = st[0][0]
    plddt = {}
    for res in chain:
        if res.seqid.num not in keys:
            continue
        b_values = [atom.b_iso for atom in res if not atom.element.is_hydrogen]
        if not b_values:
            continue
        if max(b_values) - min(b_values) > 1e-6:
            logger.warning(
                "residue %s in %s: per-atom B_iso_or_equiv (pLDDT) values disagree (%s)",
                res.seqid.num, cif_path, b_values,
            )
        plddt[res.seqid.num] = float(b_values[0])
    return plddt


def read_secondary_structure_by_auth(cif_path: Path, chain_id: str) -> tuple[dict, str]:
    """Returns ({(auth_seq_id, auth_ins_code): sse_label}, method_used).
    method_used is always "biotite_psea" in this container (mkdssp is
    confirmed unavailable, see PLAN.md Phase 8 analysis plan) -- kept as a
    return value, not a hardcoded assumption, so a future DSSP install is
    a one-line change here rather than a silent switch.
    """
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    chain = st[0][chain_id]

    keys, res_ids, res_names = [], [], []
    coords = []
    for res in chain:
        for atom in res:
            if atom.name == "CA" and not atom.element.is_hydrogen:
                keys.append((res.seqid.num, res.seqid.icode.strip()))
                res_ids.append(res.seqid.num)
                res_names.append(res.name)
                coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                break

    if not keys:
        return {}, "biotite_psea"

    n = len(keys)
    atom_array = bstruc.AtomArray(n)
    atom_array.coord = np.array(coords)
    atom_array.chain_id = np.array([chain_id] * n)
    atom_array.res_id = np.array(res_ids)
    atom_array.res_name = np.array(res_names)
    atom_array.atom_name = np.array(["CA"] * n)
    atom_array.element = np.array(["C"] * n)
    atom_array.hetero = np.zeros(n, dtype=bool)

    sse = bstruc.annotate_sse(atom_array)
    return dict(zip(keys, sse)), "biotite_psea"
