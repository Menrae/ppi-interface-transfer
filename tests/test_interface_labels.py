from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import config
from src.data import interface_labels as il

# --- mmCIF fixture builder ---------------------------------------------

ATOM_SITE_TAGS = [
    "_atom_site.group_PDB", "_atom_site.id", "_atom_site.auth_asym_id", "_atom_site.auth_seq_id",
    "_atom_site.pdbx_PDB_ins_code", "_atom_site.label_seq_id", "_atom_site.label_alt_id",
    "_atom_site.label_comp_id", "_atom_site.auth_comp_id", "_atom_site.label_asym_id",
    "_atom_site.label_entity_id", "_atom_site.type_symbol", "_atom_site.label_atom_id",
    "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.occupancy",
    "_atom_site.pdbx_PDB_model_num",
]

OPER_LIST_TAGS = [
    "_pdbx_struct_oper_list.id", "_pdbx_struct_oper_list.type", "_pdbx_struct_oper_list.name",
    "_pdbx_struct_oper_list.symmetry_operation",
    "_pdbx_struct_oper_list.matrix[1][1]", "_pdbx_struct_oper_list.matrix[1][2]",
    "_pdbx_struct_oper_list.matrix[1][3]", "_pdbx_struct_oper_list.vector[1]",
    "_pdbx_struct_oper_list.matrix[2][1]", "_pdbx_struct_oper_list.matrix[2][2]",
    "_pdbx_struct_oper_list.matrix[2][3]", "_pdbx_struct_oper_list.vector[2]",
    "_pdbx_struct_oper_list.matrix[3][1]", "_pdbx_struct_oper_list.matrix[3][2]",
    "_pdbx_struct_oper_list.matrix[3][3]", "_pdbx_struct_oper_list.vector[3]",
]

_ATOM_ID = [0]


def atom_row(chain, auth_seq_id, comp_id, atom_name, element, x, y, z, occ=1.0, alt=".", ins="?", model=1, label_asym=None):
    _ATOM_ID[0] += 1
    label_asym = label_asym or chain
    group = "HETATM" if comp_id in ("HOH", "SO4") else "ATOM"
    return (
        group, _ATOM_ID[0], chain, auth_seq_id, ins, auth_seq_id, alt, comp_id, comp_id,
        label_asym, "1", element, atom_name, f"{x:.3f}", f"{y:.3f}", f"{z:.3f}", occ, model,
    )


IDENTITY_OP = (1, "'identity operation'", "1_555", "x,y,z", 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)


def translation_op(op_id, dx, dy, dz):
    return (op_id, "'crystal symmetry operation'", f"{op_id}_555", "x,y,z", 1.0, 0.0, 0.0, dx, 0.0, 1.0, 0.0, dy, 0.0, 0.0, 1.0, dz)


def write_structure_cif(path: Path, atom_rows: list[tuple], assembly: dict | None = None) -> None:
    lines = ["data_test", "#", "loop_"] + ATOM_SITE_TAGS
    for row in atom_rows:
        lines.append(" ".join(str(v) for v in row))
    lines.append("#")
    if assembly is not None:
        lines += [
            "_pdbx_struct_assembly.id 1",
            "_pdbx_struct_assembly.details author_defined_assembly",
            "_pdbx_struct_assembly.method_details ?",
            "_pdbx_struct_assembly.oligomeric_details ?",
            "_pdbx_struct_assembly.oligomeric_count ?",
            "#",
            "_pdbx_struct_assembly_gen.assembly_id 1",
            f"_pdbx_struct_assembly_gen.oper_expression {assembly['oper_expression']}",
            f"_pdbx_struct_assembly_gen.asym_id_list {assembly['asym_id_list']}",
            "#",
            "loop_",
        ] + OPER_LIST_TAGS
        for op in assembly["operators"]:
            lines.append(" ".join(str(v) for v in op))
        lines.append("#")
    path.write_text("\n".join(lines) + "\n")


def make_phase4_df(keys: list[tuple[int, str]], uniprot_start: int = 100) -> pd.DataFrame:
    rows = []
    for i, (auth_seq_id, ins_code) in enumerate(keys):
        rows.append(
            {
                "auth_seq_id": auth_seq_id, "auth_ins_code": ins_code,
                "label_seq_id": i + 1, "residue_name": "ALA", "is_observed": True,
                "uniprot_acc": "P00000", "uniprot_resnum": uniprot_start + i,
                "alphafold_resnum": uniprot_start + i, "alphafold_residue_name": "ALA",
                "match_flag": True, "method": "sifts", "unmapped_reason": "",
            }
        )
    return pd.DataFrame(rows)


@pytest.fixture(autouse=True)
def _reset_atom_id_counter():
    _ATOM_ID[0] = 0
    yield


# --- one_letter_code ---------------------------------------------------


def test_one_letter_code_standard_and_modified():
    assert il.one_letter_code("ALA") == "A"
    assert il.one_letter_code("MSE") == "M"


def test_one_letter_code_non_amino_acid_is_none():
    assert il.one_letter_code("HOH") is None


# --- jaccard_and_confusion ------------------------------------------------


def test_jaccard_and_confusion_matches_hand_count():
    distance = [True, True, False, False, True]
    sasa_flags = [True, False, False, True, True]
    result = il.jaccard_and_confusion(distance, sasa_flags)
    assert result == {"tp": 2, "fp": 1, "fn": 1, "tn": 1, "jaccard": pytest.approx(2 / 4)}


def test_jaccard_and_confusion_empty_union_is_zero():
    assert il.jaccard_and_confusion([False, False], [False, False])["jaccard"] == 0.0


# --- compute_attrition / sort_by_priority ---------------------------------


def test_compute_attrition_per_leakage_mode():
    rows = [
        {"eligible_homolog": True, "eligible_exact_train": True, "eligible_none": True, "status": "labeled"},
        {"eligible_homolog": False, "eligible_exact_train": True, "eligible_none": True, "status": "excluded"},
    ]
    attrition = il.compute_attrition(rows)
    by_mode = {r["leakage_filter_mode"]: r for r in attrition}
    assert by_mode["homolog"]["n_labeled"] == 1
    assert by_mode["exact_train"]["n_labeled"] == 1
    assert by_mode["exact_train"]["n_excluded"] == 1


def test_sort_by_priority_primary_first():
    rows = [
        {"eligible_homolog": False, "eligible_exact_train": False, "cluster_id": "3CCC_A"},
        {"eligible_homolog": True, "eligible_exact_train": True, "cluster_id": "1AAA_A"},
    ]
    ordered = il.sort_by_priority(rows)
    assert [r["cluster_id"] for r in ordered] == ["1AAA_A", "3CCC_A"]


# --- choose_assembly -------------------------------------------------------


def test_choose_assembly_none_when_no_assemblies_defined(tmp_path):
    path = tmp_path / "test.cif"
    rows = [atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0)]
    write_structure_cif(path, rows, assembly=None)
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    assert il.choose_assembly(st, "A") is None


def test_choose_assembly_raises_when_chain_not_covered(tmp_path):
    path = tmp_path / "test.cif"
    rows = [
        atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0),
        atom_row("B", 1, "GLY", "CA", "C", 50, 50, 50),
    ]
    write_structure_cif(
        path, rows,
        assembly={"oper_expression": "1", "asym_id_list": "B", "operators": [IDENTITY_OP]},
    )
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    with pytest.raises(il.ChainExcluded) as exc:
        il.choose_assembly(st, "A")
    assert exc.value.reason == il.CHAIN_EXCL_NO_ASSEMBLY_FOR_CHAIN


# --- end-to-end: known contacts inside/outside cutoff, H/water/ligand ----


def test_contacts_just_inside_and_outside_cutoff_ignoring_h_water_ligand(tmp_path):
    cutoff = config.INTERFACE_DISTANCE_CUTOFF_ANGSTROM
    rows = [
        # Representative chain A: residue 1 at origin, residue 2 far away (isolated).
        atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0),
        atom_row("A", 2, "ALA", "CA", "C", 30, 0, 0),
        # Partner chain B: one atom just inside the cutoff of A's residue 1.
        atom_row("B", 1, "GLY", "CA", "C", cutoff - 0.1, 0, 0),
        # A hydrogen on B even closer -- must be ignored (would otherwise be
        # the true nearest neighbor to A's residue 1).
        atom_row("B", 1, "GLY", "H", "H", 0.01, 0, 0),
        # Water and a ligand extremely close to residue 2 -- must be ignored
        # entirely (never counted as a partner), so residue 2 stays non-interface.
        atom_row("B", 101, "HOH", "O", "O", 30.01, 0, 0),
        atom_row("B", 201, "SO4", "S", "S", 30.02, 0, 0),
    ]
    write_structure_cif(tmp_path / "test.cif", rows, assembly=None)
    phase4_df = make_phase4_df([(1, ""), (2, "")])

    parquet_rows, report = il.process_representative(tmp_path / "test.cif", "TEST", "A", phase4_df)

    assert report["status"] == "labeled"
    by_seq = {r["auth_seq_id"]: r for r in parquet_rows}
    assert by_seq[1]["is_interface_contact"] is True
    assert by_seq[1]["min_partner_distance"] == pytest.approx(cutoff - 0.1, abs=1e-3)
    assert by_seq[2]["is_interface_contact"] is False  # water/ligand proximity must not count


def test_contact_exactly_at_cutoff_boundary_and_just_outside(tmp_path):
    cutoff = config.INTERFACE_DISTANCE_CUTOFF_ANGSTROM
    rows = [
        atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0),
        atom_row("B", 1, "GLY", "CA", "C", cutoff + 0.1, 0, 0),
    ]
    write_structure_cif(tmp_path / "test.cif", rows, assembly=None)
    phase4_df = make_phase4_df([(1, "")])

    parquet_rows, report = il.process_representative(tmp_path / "test.cif", "TEST", "A", phase4_df)

    assert report["status"] == "excluded"
    assert report["exclusion_reason"] == il.CHAIN_EXCL_ZERO_INTERFACE_RESIDUES


# --- end-to-end: partner only via assembly symmetry operator -------------


def test_partner_only_exists_after_assembly_symmetry_operator(tmp_path):
    rows = [
        atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0),
        atom_row("A", 2, "ALA", "CA", "C", 0, 0, 20),
    ]
    assembly = {
        "oper_expression": "1,2",
        "asym_id_list": "A",
        "operators": [IDENTITY_OP, translation_op(2, 3.0, 0.0, 0.0)],
    }
    write_structure_cif(tmp_path / "test.cif", rows, assembly=assembly)
    phase4_df = make_phase4_df([(1, ""), (2, "")])

    # Sanity: the raw asymmetric unit alone (no assembly expansion) has no
    # protein partner at all for chain A -- only the assembly provides one.
    import gemmi
    st = gemmi.read_structure(str(tmp_path / "test.cif"))
    st.setup_entities()
    assert il.count_asu_only_partner_chains(st, "A") == 0

    parquet_rows, report = il.process_representative(tmp_path / "test.cif", "TEST", "A", phase4_df)

    assert report["status"] == "labeled"
    assert report["n_partner_chains_asu_only"] == 0
    assert report["would_exclude_if_asu"] is True
    assert report["has_homomeric_partner"] is True
    assert report["n_partner_chains"] == 1
    # Every representative residue finds its own symmetry-translated copy
    # 3.0 A away -- well inside the cutoff.
    assert all(r["is_interface_contact"] for r in parquet_rows)


# --- end-to-end: ASU crystal contact not in the chosen assembly ----------


def test_asu_contact_excluded_when_not_in_assembly(tmp_path):
    rows = [
        atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0),
        # Chain C sits right next to A in the raw file but is NOT part of
        # the assembly's asym_id_list below -- a crystal-packing contact.
        atom_row("C", 1, "GLY", "CA", "C", 2.0, 0, 0),
        # Chain D IS part of the assembly, a bit further but still in range,
        # so the chain is still labeled (not excluded for lack of any partner).
        atom_row("D", 1, "SER", "CA", "C", 3.5, 0, 0),
    ]
    assembly = {"oper_expression": "1", "asym_id_list": "A,D", "operators": [IDENTITY_OP]}
    write_structure_cif(tmp_path / "test.cif", rows, assembly=assembly)
    phase4_df = make_phase4_df([(1, "")])

    parquet_rows, report = il.process_representative(tmp_path / "test.cif", "TEST", "A", phase4_df)

    assert report["status"] == "labeled"
    assert report["n_partner_chains"] == 1  # D only, not C
    # Chain C would have been closer than D -- if it had leaked in, the
    # nearest-partner distance would be 2.0, not 3.5.
    assert parquet_rows[0]["min_partner_distance"] == pytest.approx(3.5, abs=1e-3)
    assert parquet_rows[0]["is_interface_contact"] is True


# --- exact join between labels and the residue map ------------------------


def test_representative_residues_matches_phase4_keys_exactly(tmp_path):
    rows = [atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0), atom_row("A", 2, "ALA", "CA", "C", 5, 0, 0)]
    write_structure_cif(tmp_path / "test.cif", rows, assembly=None)
    import gemmi
    st = gemmi.read_structure(str(tmp_path / "test.cif"))
    st.setup_entities()
    residues = il.representative_residues(st[0]["A"], {(1, ""), (2, "")})
    assert set(residues) == {(1, ""), (2, "")}


def test_representative_residues_raises_on_missing_key(tmp_path):
    rows = [atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0)]
    write_structure_cif(tmp_path / "test.cif", rows, assembly=None)
    import gemmi
    st = gemmi.read_structure(str(tmp_path / "test.cif"))
    st.setup_entities()
    # phase4_keys claims a residue (auth_seq_id=2) that isn't actually there.
    with pytest.raises(il.ChainExcluded) as exc:
        il.representative_residues(st[0]["A"], {(1, ""), (2, "")})
    assert exc.value.reason.startswith(il.CHAIN_EXCL_RESIDUE_KEY_MISMATCH)


def test_representative_residues_raises_on_extra_geometry_residue(tmp_path):
    rows = [atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0), atom_row("A", 2, "ALA", "CA", "C", 5, 0, 0)]
    write_structure_cif(tmp_path / "test.cif", rows, assembly=None)
    import gemmi
    st = gemmi.read_structure(str(tmp_path / "test.cif"))
    st.setup_entities()
    # phase4_keys is missing residue 2, which IS present in the structure.
    with pytest.raises(il.ChainExcluded) as exc:
        il.representative_residues(st[0]["A"], {(1, "")})
    assert exc.value.reason.startswith(il.CHAIN_EXCL_RESIDUE_KEY_MISMATCH)


def test_process_representative_writes_exact_join_end_to_end(tmp_path):
    rows = [
        atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0),
        atom_row("A", 2, "ALA", "CA", "C", 5, 0, 0),
        atom_row("B", 1, "GLY", "CA", "C", 1.0, 0, 0),
    ]
    write_structure_cif(tmp_path / "test.cif", rows, assembly=None)
    phase4_df = make_phase4_df([(1, ""), (2, "")])

    parquet_rows, report = il.process_representative(tmp_path / "test.cif", "TEST", "A", phase4_df)
    assert report["status"] == "labeled"
    assert {(r["auth_seq_id"], r["auth_ins_code"]) for r in parquet_rows} == {(1, ""), (2, "")}
    assert {r["uniprot_resnum"] for r in parquet_rows} == {100, 101}


# --- alt-loc resolution (occupancy tie-break) via resolve_heavy_atoms -----


def test_resolve_heavy_atoms_picks_highest_occupancy_altloc(tmp_path):
    rows = [
        atom_row("A", 1, "HIS", "CA", "C", 0, 0, 0, occ=0.4, alt="A"),
        atom_row("A", 1, "HIS", "CA", "C", 1, 1, 1, occ=0.6, alt="B"),
    ]
    write_structure_cif(tmp_path / "test.cif", rows, assembly=None)
    import gemmi
    st = gemmi.read_structure(str(tmp_path / "test.cif"))
    st.setup_entities()
    res = list(st[0]["A"])[0]
    atoms = il.resolve_heavy_atoms(res)
    assert len(atoms) == 1
    assert atoms[0].pos.tolist() == [1.0, 1.0, 1.0]  # occupancy 0.6 wins


def test_resolve_heavy_atoms_ignores_hydrogen(tmp_path):
    rows = [
        atom_row("A", 1, "ALA", "CA", "C", 0, 0, 0),
        atom_row("A", 1, "ALA", "HA", "H", 0, 0, 1),
    ]
    write_structure_cif(tmp_path / "test.cif", rows, assembly=None)
    import gemmi
    st = gemmi.read_structure(str(tmp_path / "test.cif"))
    st.setup_entities()
    res = list(st[0]["A"])[0]
    atoms = il.resolve_heavy_atoms(res)
    assert len(atoms) == 1
    assert atoms[0].element == "C"
