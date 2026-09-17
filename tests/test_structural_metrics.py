"""Tests for src.analysis.structural_metrics. No network access."""

from pathlib import Path

import numpy as np
import pytest

from src.analysis import structural_metrics as sm

AA3 = "ALA"


def write_cif(path: Path, chain: str, res_ids, coords, b_isos=None, res_names=None) -> None:
    """Minimal single-chain, single-model mmCIF with real coordinates and
    B_iso_or_equiv values, parseable by gemmi.read_structure -- mirrors the
    tag set gemmi needs (confirmed against a real AlphaFold DB mmCIF this
    session; align_residues.py's own fixtures omit coordinates entirely
    since Phase 4 never needs geometry)."""
    if b_isos is None:
        b_isos = [50.0] * len(res_ids)
    if res_names is None:
        res_names = [AA3] * len(res_ids)
    tags = [
        "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
        "_atom_site.label_atom_id", "_atom_site.label_alt_id", "_atom_site.label_comp_id",
        "_atom_site.label_asym_id", "_atom_site.label_entity_id", "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
        "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_formal_charge", "_atom_site.auth_seq_id", "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id", "_atom_site.auth_atom_id", "_atom_site.pdbx_PDB_model_num",
    ]
    lines = ["data_test", "loop_"] + tags
    for i, (res_id, (x, y, z)) in enumerate(zip(res_ids, coords)):
        lines.append(
            f"ATOM {i + 1} C CA . {res_names[i]} {chain} 1 {res_id} ? "
            f"{x:.3f} {y:.3f} {z:.3f} 1.00 {b_isos[i]:.2f} ? {res_id} {res_names[i]} {chain} CA 1"
        )
    lines.append("#")
    path.write_text("\n".join(lines) + "\n")


def helix_like_coords(n: int, spacing: float = 3.8) -> np.ndarray:
    """A gently curved backbone (not collinear) so Kabsch fits are
    well-defined, spaced ~one Calpha-Calpha bond length apart."""
    t = np.arange(n) * 0.3
    x = spacing * np.arange(n)
    y = 3.0 * np.sin(t)
    z = 3.0 * np.cos(t)
    return np.column_stack([x, y, z])


# --- kabsch_superpose -------------------------------------------------------


def test_kabsch_superpose_recovers_known_rotation():
    rng = np.random.default_rng(0)
    points = rng.normal(size=(20, 3))
    theta = 0.7
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ])
    translation = np.array([5.0, -2.0, 1.0])
    target = (rotation @ points.T).T + translation

    superposed = sm.kabsch_superpose(points, target)
    assert superposed == pytest.approx(target, abs=1e-8)


# --- local_rmsd_per_residue --------------------------------------------------


def test_local_rmsd_near_zero_for_rigidly_moved_copy():
    exp_coords = helix_like_coords(40)
    theta = 0.4
    rotation = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)],
    ])
    af_coords = (rotation @ exp_coords.T).T + np.array([10.0, -5.0, 2.0])

    local_rmsd, n_neighbors, reasons = sm.local_rmsd_per_residue(exp_coords, af_coords, radius=10.0, min_neighbors=3)

    assert np.all(n_neighbors >= 3)
    assert all(r == "" for r in reasons)
    assert np.nanmax(local_rmsd) < 1e-6


def test_local_rmsd_large_only_near_locally_perturbed_region():
    exp_coords = helix_like_coords(60)
    af_coords = exp_coords.copy()  # identity "prediction" -- isolates the perturbation's effect

    perturbed_idx = slice(28, 33)
    rng = np.random.default_rng(1)
    af_coords[perturbed_idx] += rng.normal(scale=8.0, size=(5, 3))

    local_rmsd, n_neighbors, reasons = sm.local_rmsd_per_residue(exp_coords, af_coords, radius=10.0, min_neighbors=3)

    # Far from the perturbation (neighborhoods, defined on EXPERIMENTAL
    # distances, don't reach the perturbed residues at all): near zero.
    assert local_rmsd[0] < 0.5
    assert local_rmsd[-1] < 0.5
    # At the perturbed region itself: clearly elevated.
    assert local_rmsd[30] > 2.0
    assert local_rmsd[30] > local_rmsd[0]
    assert local_rmsd[30] > local_rmsd[-1]


def test_local_rmsd_too_few_neighbors_flagged_not_dropped():
    # Two isolated residues, far apart -- radius doesn't reach any neighbor
    # beyond self for either, well under min_neighbors=3.
    exp_coords = np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]])
    af_coords = exp_coords.copy()
    local_rmsd, n_neighbors, reasons = sm.local_rmsd_per_residue(exp_coords, af_coords, radius=10.0, min_neighbors=3)
    assert np.all(np.isnan(local_rmsd))
    assert list(n_neighbors) == [1, 1]
    assert all(r == "too_few_neighbors" for r in reasons)


# --- global_ca_distance_per_residue -----------------------------------------


def test_global_ca_distance_zero_for_translated_copy():
    exp_coords = helix_like_coords(15)
    af_coords = exp_coords + np.array([3.0, -1.0, 0.5])
    distances = sm.global_ca_distance_per_residue(exp_coords, af_coords)
    assert distances == pytest.approx(np.zeros(15), abs=1e-8)


def test_global_ca_distance_large_for_locally_perturbed_residue():
    exp_coords = helix_like_coords(15)
    af_coords = exp_coords.copy()
    af_coords[7] += np.array([15.0, 0.0, 0.0])
    distances = sm.global_ca_distance_per_residue(exp_coords, af_coords)
    assert distances[7] > distances[0]
    assert distances[7] > 5.0


# --- assign_plddt_band -------------------------------------------------------


def test_assign_plddt_band_boundary_values():
    values = np.array([0.0, 49.999, 50.0, 69.999, 70.0, 89.999, 90.0, 100.0])
    bands = sm.assign_plddt_band(values, bands=(50, 70, 90))
    assert list(bands) == ["<50", "<50", "50-70", "50-70", "70-90", "70-90", ">=90", ">=90"]


# --- I/O helpers -------------------------------------------------------------


def test_read_ca_coords_by_auth(tmp_path):
    path = tmp_path / "exp.cif"
    coords = helix_like_coords(3)
    write_cif(path, "A", [10, 11, 12], coords)
    result = sm.read_ca_coords_by_auth(path, "A", [(10, ""), (11, ""), (12, ""), (99, "")])
    assert set(result.keys()) == {(10, ""), (11, ""), (12, "")}
    assert result[(10, "")] == pytest.approx(coords[0])


def test_read_ca_coords_by_uniprot(tmp_path):
    path = tmp_path / "af.cif"
    coords = helix_like_coords(3)
    write_cif(path, "A", [1, 2, 3], coords)
    result = sm.read_ca_coords_by_uniprot(path, [1, 2, 3, 50])
    assert set(result.keys()) == {1, 2, 3}
    assert result[2] == pytest.approx(coords[1], abs=1e-3)


def test_read_plddt_by_uniprot(tmp_path):
    path = tmp_path / "af.cif"
    coords = helix_like_coords(3)
    write_cif(path, "A", [1, 2, 3], coords, b_isos=[42.5, 91.0, 55.5])
    result = sm.read_plddt_by_uniprot(path, [1, 2, 3])
    assert result == {1: pytest.approx(42.5), 2: pytest.approx(91.0), 3: pytest.approx(55.5)}


def test_read_secondary_structure_by_auth_returns_method_and_labels(tmp_path):
    path = tmp_path / "exp.cif"
    coords = helix_like_coords(10)
    write_cif(path, "A", list(range(1, 11)), coords)
    labels, method = sm.read_secondary_structure_by_auth(path, "A")
    assert method == "biotite_psea"
    assert set(labels.keys()) == {(i, "") for i in range(1, 11)}
    assert all(v in ("a", "b", "c", "") for v in labels.values())
