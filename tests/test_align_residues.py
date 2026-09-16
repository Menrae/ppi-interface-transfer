from pathlib import Path

import pandas as pd
import pytest

from src import config
from src.data import align_residues as ar

# --- mmCIF fixture builders -------------------------------------------------

_TAGS_WITH_SIFTS = ar._ATOM_SITE_BASE_TAGS + ar._ATOM_SITE_SIFTS_TAGS
_TAGS_NO_SIFTS = ar._ATOM_SITE_BASE_TAGS

AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "E": "GLU",
    "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def write_cif(path: Path, tags: list[str], rows: list[tuple]) -> None:
    lines = ["data_test", "loop_"] + list(tags)
    for row in rows:
        lines.append(" ".join(str(v) for v in row))
    lines.append("#")
    path.write_text("\n".join(lines) + "\n")


def atom_row(
    auth_seq_id,
    comp_id="ALA",
    ins_code=".",
    label_seq_id=None,
    alt_id=".",
    occupancy="1.000",
    model_num="1",
    sifts_name="?",
    sifts_acc="?",
    sifts_num="?",
    group="ATOM",
    chain="A",
    with_sifts=True,
):
    if label_seq_id is None:
        label_seq_id = auth_seq_id
    row = (group, chain, auth_seq_id, ins_code, label_seq_id, alt_id, comp_id, comp_id, occupancy, model_num)
    if with_sifts:
        row = row + (sifts_name, sifts_acc, sifts_num)
    return row


def write_experimental_cif(path: Path, rows: list[tuple], with_sifts: bool = True) -> None:
    tags = _TAGS_WITH_SIFTS if with_sifts else _TAGS_NO_SIFTS
    write_cif(path, tags, rows)


def write_alphafold_cif(
    path: Path, sequence: str, start: int, accession: str, chain: str = "A", model_num: str = "1"
) -> None:
    rows = []
    for offset, letter in enumerate(sequence):
        resnum = start + offset
        rows.append(
            atom_row(
                resnum, comp_id=AA3[letter], chain=chain, model_num=model_num,
                sifts_name="UNP", sifts_acc=accession, sifts_num=resnum,
            )
        )
    write_experimental_cif(path, rows)


def sifts_row(auth_seq_id, uniprot_pos, comp_id="ALA", accession="P00000", **kwargs):
    return atom_row(auth_seq_id, comp_id=comp_id, sifts_name="UNP", sifts_acc=accession, sifts_num=uniprot_pos, **kwargs)


def unmapped_row(auth_seq_id, comp_id="ALA", **kwargs):
    return atom_row(auth_seq_id, comp_id=comp_id, sifts_name="?", sifts_acc="?", sifts_num="?", **kwargs)


@pytest.fixture(autouse=True)
def _default_af_provider_cache(tmp_path, monkeypatch):
    """process_representative's first step checks AlphaFold DB's cached
    providerId (Phase 3 JSON) before touching any mmCIF -- pre-populate a
    "GDM" (official) entry for every accession used across this file's
    fixtures so tests that aren't specifically about provider checking
    aren't incidentally blocked by it."""
    cache_dir = tmp_path / "af_api_cache"
    monkeypatch.setattr(ar, "ALPHAFOLD_API_CACHE_DIR", cache_dir)
    for acc in ("P00000", "Q99999", "WRONGACC"):
        _write_af_api_cache(acc, "GDM")
    return cache_dir


# --- one_letter_code -----------------------------------------------------


def test_one_letter_code_standard_residue():
    assert ar.one_letter_code("ALA") == "A"


def test_one_letter_code_modified_residue_resolves_to_parent():
    assert ar.one_letter_code("MSE") == "M"  # selenomethionine -> MET


def test_one_letter_code_non_amino_acid_is_none():
    assert ar.one_letter_code("HOH") is None
    assert ar.one_letter_code("NAG") is None


# --- reduce_to_residues: alt-locs, ins codes, non-polymer exclusion -------


def test_reduce_to_residues_picks_highest_occupancy_altloc():
    rows = [
        {"pdbx_PDB_model_num": "1", "label_seq_id": "5", "auth_seq_id": "5", "pdbx_PDB_ins_code": "?",
         "label_alt_id": "A", "auth_comp_id": "HIS", "occupancy": "0.4",
         "pdbx_sifts_xref_db_name": "UNP", "pdbx_sifts_xref_db_acc": "P00000", "pdbx_sifts_xref_db_num": "5"},
        {"pdbx_PDB_model_num": "1", "label_seq_id": "5", "auth_seq_id": "5", "pdbx_PDB_ins_code": "?",
         "label_alt_id": "B", "auth_comp_id": "GLU", "occupancy": "0.6",
         "pdbx_sifts_xref_db_name": "UNP", "pdbx_sifts_xref_db_acc": "P00000", "pdbx_sifts_xref_db_num": "5"},
    ]
    residues, model_nums = ar.reduce_to_residues(rows)
    assert model_nums == {"1"}
    assert len(residues) == 1
    assert residues[0].comp_id == "GLU"  # 0.6 occupancy beats 0.4


def test_reduce_to_residues_tie_occupancy_breaks_by_alt_id():
    rows = [
        {"pdbx_PDB_model_num": "1", "label_seq_id": "5", "auth_seq_id": "5", "pdbx_PDB_ins_code": "?",
         "label_alt_id": "B", "auth_comp_id": "GLU", "occupancy": "0.5",
         "pdbx_sifts_xref_db_name": "?", "pdbx_sifts_xref_db_acc": "?", "pdbx_sifts_xref_db_num": "?"},
        {"pdbx_PDB_model_num": "1", "label_seq_id": "5", "auth_seq_id": "5", "pdbx_PDB_ins_code": "?",
         "label_alt_id": "A", "auth_comp_id": "HIS", "occupancy": "0.5",
         "pdbx_sifts_xref_db_name": "?", "pdbx_sifts_xref_db_acc": "?", "pdbx_sifts_xref_db_num": "?"},
    ]
    residues, _ = ar.reduce_to_residues(rows)
    assert residues[0].comp_id == "HIS"  # alt_id "A" sorts before "B"


def test_reduce_to_residues_excludes_non_polymer_rows():
    rows = [
        {"pdbx_PDB_model_num": "1", "label_seq_id": ".", "auth_seq_id": "501", "pdbx_PDB_ins_code": "?",
         "label_alt_id": ".", "auth_comp_id": "HOH", "occupancy": "1.0",
         "pdbx_sifts_xref_db_name": "?", "pdbx_sifts_xref_db_acc": "?", "pdbx_sifts_xref_db_num": "?"},
        {"pdbx_PDB_model_num": "1", "label_seq_id": "1", "auth_seq_id": "1", "pdbx_PDB_ins_code": "?",
         "label_alt_id": ".", "auth_comp_id": "ALA", "occupancy": "1.0",
         "pdbx_sifts_xref_db_name": "UNP", "pdbx_sifts_xref_db_acc": "P00000", "pdbx_sifts_xref_db_num": "1"},
    ]
    residues, _ = ar.reduce_to_residues(rows)
    assert len(residues) == 1
    assert residues[0].comp_id == "ALA"


def test_reduce_to_residues_insertion_code_disambiguates_same_auth_seq_id():
    rows = [
        {"pdbx_PDB_model_num": "1", "label_seq_id": "10", "auth_seq_id": "10", "pdbx_PDB_ins_code": "?",
         "label_alt_id": ".", "auth_comp_id": "ALA", "occupancy": "1.0",
         "pdbx_sifts_xref_db_name": "UNP", "pdbx_sifts_xref_db_acc": "P00000", "pdbx_sifts_xref_db_num": "10"},
        {"pdbx_PDB_model_num": "1", "label_seq_id": "11", "auth_seq_id": "10", "pdbx_PDB_ins_code": "A",
         "label_alt_id": ".", "auth_comp_id": "GLY", "occupancy": "1.0",
         "pdbx_sifts_xref_db_name": "UNP", "pdbx_sifts_xref_db_acc": "P00000", "pdbx_sifts_xref_db_num": "11"},
    ]
    residues, _ = ar.reduce_to_residues(rows)
    assert len(residues) == 2
    by_ins = {r.ins_code: r for r in residues}
    assert by_ins[""].comp_id == "ALA"
    assert by_ins["A"].comp_id == "GLY"


# --- load_chain_residues / load_alphafold_residues (I/O, fixture files) --


def test_load_chain_residues_reads_fixture_cif(tmp_path):
    path = tmp_path / "test.cif"
    rows = [sifts_row(1, 1), sifts_row(2, 2), sifts_row(3, 3)]
    write_experimental_cif(path, rows)
    residues = ar.load_chain_residues(path, "A")
    assert [r.auth_seq_id for r in residues] == [1, 2, 3]


def test_load_chain_residues_missing_chain_raises_chain_not_found(tmp_path):
    path = tmp_path / "test.cif"
    write_experimental_cif(path, [sifts_row(1, 1)])
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.load_chain_residues(path, "Z")
    assert exc.value.reason == ar.CHAIN_EXCL_CHAIN_NOT_FOUND


def test_load_chain_residues_multiple_models_raises(tmp_path):
    path = tmp_path / "test.cif"
    rows = [sifts_row(1, 1, model_num="1"), sifts_row(1, 1, model_num="2")]
    write_experimental_cif(path, rows)
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.load_chain_residues(path, "A")
    assert exc.value.reason.startswith(ar.CHAIN_EXCL_MULTIPLE_MODELS)


def test_load_alphafold_residues_multichain_raises(tmp_path):
    path = tmp_path / "af.cif"
    rows = [
        atom_row(1, chain="A", sifts_name="UNP", sifts_acc="P00000", sifts_num=1),
        atom_row(1, chain="B", sifts_name="UNP", sifts_acc="P00000", sifts_num=1),
    ]
    write_experimental_cif(path, rows)
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.load_alphafold_residues(path)
    assert exc.value.reason.startswith(ar.CHAIN_EXCL_MULTICHAIN_ALPHAFOLD)


def test_load_raw_atom_site_rows_handles_missing_sifts_columns(tmp_path):
    path = tmp_path / "test.cif"
    rows = [atom_row(1, with_sifts=False)]
    write_experimental_cif(path, rows, with_sifts=False)
    residues, _ = ar.reduce_to_residues(ar.load_raw_atom_site_rows(path, "A"))
    assert residues[0].sifts_db_name is None


# --- verify_official_alphafold_provider ------------------------------------


def _write_af_api_cache(acc: str, provider_id: str | None, matches_accession: bool = True):
    ar.ALPHAFOLD_API_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = ar.ALPHAFOLD_API_CACHE_DIR / f"{acc}.json"
    entry_acc = acc if matches_accession else f"{acc}-2"
    body = [{"uniprotAccession": entry_acc, "providerId": provider_id}] if provider_id is not None else []
    path.write_text(__import__("json").dumps({"status_code": 200, "body": body}))
    return path


def test_verify_official_alphafold_provider_accepts_gdm(tmp_path, monkeypatch):
    monkeypatch.setattr(ar, "ALPHAFOLD_API_CACHE_DIR", tmp_path)
    _write_af_api_cache("P00000", "GDM")
    ar.verify_official_alphafold_provider("P00000")  # does not raise


def test_verify_official_alphafold_provider_rejects_community_submission(tmp_path, monkeypatch):
    monkeypatch.setattr(ar, "ALPHAFOLD_API_CACHE_DIR", tmp_path)
    _write_af_api_cache("P00000", "VR3D")
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.verify_official_alphafold_provider("P00000")
    assert exc.value.reason.startswith(ar.CHAIN_EXCL_NOT_OFFICIAL_ALPHAFOLD)
    assert "VR3D" in exc.value.reason


def test_verify_official_alphafold_provider_rejects_when_no_exact_match(tmp_path, monkeypatch):
    monkeypatch.setattr(ar, "ALPHAFOLD_API_CACHE_DIR", tmp_path)
    _write_af_api_cache("P00000", "GDM", matches_accession=False)
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.verify_official_alphafold_provider("P00000")
    assert exc.value.reason.startswith(ar.CHAIN_EXCL_NOT_OFFICIAL_ALPHAFOLD)


# --- validate_alphafold_numbering (AF side, never assumed) ----------------


def test_validate_alphafold_numbering_success(tmp_path):
    path = tmp_path / "af.cif"
    write_alphafold_cif(path, "MSSTLHS", start=1, accession="P00000")
    residues = ar.load_alphafold_residues(path)
    by_pos = ar.validate_alphafold_numbering(residues, "P00000", 1, 7)
    assert by_pos[1].comp_id == "MET"
    assert by_pos[7].comp_id == "SER"


def test_validate_alphafold_numbering_range_mismatch_raises(tmp_path):
    path = tmp_path / "af.cif"
    write_alphafold_cif(path, "MSSTLHS", start=1, accession="P00000")
    residues = ar.load_alphafold_residues(path)
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.validate_alphafold_numbering(residues, "P00000", 1, 100)  # claims a longer range than the file has
    assert exc.value.reason.startswith(ar.CHAIN_EXCL_ALPHAFOLD_RANGE_MISMATCH)


def test_validate_alphafold_numbering_wrong_accession_raises(tmp_path):
    path = tmp_path / "af.cif"
    write_alphafold_cif(path, "MSSTLHS", start=1, accession="P00000")
    residues = ar.load_alphafold_residues(path)
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.validate_alphafold_numbering(residues, "Q11111", 1, 7)  # different accession than embedded
    assert exc.value.reason == ar.CHAIN_EXCL_ALPHAFOLD_NOT_UNIPROT_NATIVE


def test_validate_alphafold_numbering_missing_sifts_columns_raises(tmp_path):
    path = tmp_path / "af.cif"
    rows = [atom_row(i, comp_id="ALA", with_sifts=False) for i in range(1, 6)]
    write_experimental_cif(path, rows, with_sifts=False)
    residues = ar.load_alphafold_residues(path)
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.validate_alphafold_numbering(residues, "P00000", 1, 5)
    assert exc.value.reason == ar.CHAIN_EXCL_ALPHAFOLD_NOT_UNIPROT_NATIVE


# --- map_via_sifts ----------------------------------------------------------


def test_map_via_sifts_maps_and_matches(tmp_path):
    af_path = tmp_path / "af.cif"
    write_alphafold_cif(af_path, "MSSTL", start=1, accession="P00000")
    af_by_pos = ar.validate_alphafold_numbering(ar.load_alphafold_residues(af_path), "P00000", 1, 5)

    exp_residues = [
        ar.Residue(1, "", 1, "MET", "UNP", "P00000", 1),
        ar.Residue(2, "", 2, "SER", "UNP", "P00000", 2),
    ]
    mapped = ar.map_via_sifts(exp_residues, af_by_pos, "P00000")
    assert [m.uniprot_resnum for m in mapped] == [1, 2]
    assert all(m.match_flag for m in mapped)
    assert all(m.unmapped_reason == "" for m in mapped)


def test_map_via_sifts_no_mapping_reason():
    exp_residues = [ar.Residue(1, "", 1, "ALA", None, None, None)]
    mapped = ar.map_via_sifts(exp_residues, {}, "P00000")
    assert mapped[0].unmapped_reason == ar.REASON_NO_SIFTS_MAPPING


def test_map_via_sifts_wrong_accession_reason():
    exp_residues = [ar.Residue(1, "", 1, "ALA", "UNP", "Q99999", 1)]
    mapped = ar.map_via_sifts(exp_residues, {}, "P00000")
    assert mapped[0].unmapped_reason == ar.REASON_SIFTS_WRONG_ACCESSION


def test_map_via_sifts_outside_alphafold_range_reason():
    exp_residues = [ar.Residue(1, "", 1, "ALA", "UNP", "P00000", 500)]
    mapped = ar.map_via_sifts(exp_residues, {1: ar.Residue(1, "", 1, "ALA", "UNP", "P00000", 1)}, "P00000")
    assert mapped[0].unmapped_reason == ar.REASON_OUTSIDE_AF_RANGE


def test_map_via_sifts_non_amino_acid_reason():
    exp_residues = [ar.Residue(1, "", 1, "NAG", "UNP", "P00000", 1)]
    af_by_pos = {1: ar.Residue(1, "", 1, "ALA", "UNP", "P00000", 1)}
    mapped = ar.map_via_sifts(exp_residues, af_by_pos, "P00000")
    assert mapped[0].unmapped_reason == ar.REASON_NON_AMINO_ACID


def test_map_via_sifts_modified_residue_matches_parent():
    # MSE (selenomethionine) at a position where AlphaFold has plain MET.
    exp_residues = [ar.Residue(1, "", 1, "MSE", "UNP", "P00000", 1)]
    af_by_pos = {1: ar.Residue(1, "", 1, "MET", "UNP", "P00000", 1)}
    mapped = ar.map_via_sifts(exp_residues, af_by_pos, "P00000")
    assert mapped[0].match_flag is True
    assert mapped[0].residue.comp_id == "MSE"


def test_map_via_sifts_records_engineered_mutation_without_failing():
    # SIFTS correctly maps this residue to position 3, but the observed
    # residue (ALA) differs from AlphaFold's (GLY) there -- an engineered
    # point mutation, not a mapping error: kept, flagged, not excluded.
    exp_residues = [ar.Residue(3, "", 3, "ALA", "UNP", "P00000", 3)]
    af_by_pos = {3: ar.Residue(3, "", 3, "GLY", "UNP", "P00000", 3)}
    mapped = ar.map_via_sifts(exp_residues, af_by_pos, "P00000")
    assert mapped[0].match_flag is False
    assert mapped[0].uniprot_resnum == 3
    assert mapped[0].unmapped_reason == ""  # kept, not excluded


# --- naive offset mapping would be wrong; ours is not ---------------------


def test_sifts_mapping_correct_where_naive_offset_would_be_wrong(tmp_path):
    """auth_seq_id jumps non-sequentially (-2..1, then 100..102) while the
    true UniProt numbering is a plain contiguous 1..6. A naive
    "uniprot = auth_seq_id + constant_offset" derived from the first
    residue would compute a wrong position for the later, jumped block;
    the SIFTS-driven mapping (no offset assumption anywhere) gets it right.
    """
    af_path = tmp_path / "af.cif"
    write_alphafold_cif(af_path, "MSSTLK", start=1, accession="P00000")
    af_by_pos = ar.validate_alphafold_numbering(ar.load_alphafold_residues(af_path), "P00000", 1, 6)

    exp_residues = [
        ar.Residue(-2, "", -2, "MET", "UNP", "P00000", 1),
        ar.Residue(-1, "", -1, "SER", "UNP", "P00000", 2),
        ar.Residue(0, "", 0, "SER", "UNP", "P00000", 3),
        ar.Residue(1, "", 1, "THR", "UNP", "P00000", 4),
        ar.Residue(100, "", 100, "LEU", "UNP", "P00000", 5),
        ar.Residue(101, "", 101, "LYS", "UNP", "P00000", 6),
    ]
    mapped = ar.map_via_sifts(exp_residues, af_by_pos, "P00000")

    naive_offset = exp_residues[0].sifts_db_num - exp_residues[0].auth_seq_id  # = 3
    naive_prediction_for_last = exp_residues[-1].auth_seq_id + naive_offset  # 101 + 3 = 104 (wrong)
    assert naive_prediction_for_last != mapped[-1].uniprot_resnum
    assert mapped[-1].uniprot_resnum == 6
    assert all(m.match_flag for m in mapped)


# --- map_via_fallback --------------------------------------------------


def test_map_via_fallback_aligns_identical_sequence():
    exp_residues = [ar.Residue(i, "", i, AA3[c], None, None, None) for i, c in enumerate("MSSTLHS", start=1)]
    af_by_pos = {i: ar.Residue(i, "", i, AA3[c], "UNP", "P00000", i) for i, c in enumerate("MSSTLHS", start=1)}
    mapped = ar.map_via_fallback(exp_residues, af_by_pos)
    assert [m.uniprot_resnum for m in mapped] == [1, 2, 3, 4, 5, 6, 7]
    assert all(m.match_flag for m in mapped)


def test_map_via_fallback_trims_expression_tag():
    # 4-residue His-ish tag (no SIFTS mapping in real life either) followed
    # by the real sequence, which exactly matches AlphaFold's.
    tag = "HHHH"
    real = "MSSTLHS"
    exp_residues = [
        ar.Residue(i, "", i, AA3[c], None, None, None)
        for i, c in enumerate(tag + real, start=1)
    ]
    af_by_pos = {i: ar.Residue(i, "", i, AA3[c], "UNP", "P00000", i) for i, c in enumerate(real, start=1)}
    mapped = ar.map_via_fallback(exp_residues, af_by_pos)

    tag_mapped = mapped[: len(tag)]
    real_mapped = mapped[len(tag):]
    assert all(m.unmapped_reason == ar.REASON_FALLBACK_NO_PARTNER for m in tag_mapped)
    assert [m.uniprot_resnum for m in real_mapped] == [1, 2, 3, 4, 5, 6, 7]
    assert all(m.match_flag for m in real_mapped)


# --- choose_mapping: SIFTS vs fallback selection ---------------------------


def test_choose_mapping_uses_sifts_when_it_validates():
    exp_residues = [ar.Residue(i, "", i, "ALA", "UNP", "P00000", i) for i in range(1, 21)]
    af_by_pos = {i: ar.Residue(i, "", i, "ALA", "UNP", "P00000", i) for i in range(1, 21)}
    mapped, method = ar.choose_mapping(exp_residues, af_by_pos, "P00000")
    assert method == "sifts"
    assert all(m.match_flag for m in mapped)


def test_choose_mapping_falls_back_when_sifts_absent():
    exp_residues = [ar.Residue(i, "", i, AA3[c], None, None, None) for i, c in enumerate("MSSTLHSVFFT", start=1)]
    af_by_pos = {i: ar.Residue(i, "", i, AA3[c], "UNP", "P00000", i) for i, c in enumerate("MSSTLHSVFFT", start=1)}
    mapped, method = ar.choose_mapping(exp_residues, af_by_pos, "P00000")
    assert method == "fallback"
    assert all(m.match_flag for m in mapped)


def test_choose_mapping_falls_back_when_sifts_wrong_accession_for_whole_chain():
    real = "MSSTLHSVFFT"
    exp_residues = [
        ar.Residue(i, "", i, AA3[c], "UNP", "WRONGACC", i) for i, c in enumerate(real, start=1)
    ]
    af_by_pos = {i: ar.Residue(i, "", i, AA3[c], "UNP", "P00000", i) for i, c in enumerate(real, start=1)}
    mapped, method = ar.choose_mapping(exp_residues, af_by_pos, "P00000")
    assert method == "fallback"
    assert all(m.match_flag for m in mapped)


def test_choose_mapping_raises_when_both_methods_fail():
    # Completely unrelated sequence with no correct SIFTS mapping either --
    # neither method can validate.
    exp_residues = [ar.Residue(i, "", i, "TRP", None, None, None) for i in range(1, 21)]
    af_by_pos = {i: ar.Residue(i, "", i, "ASP", "UNP", "P00000", i) for i in range(1, 21)}
    with pytest.raises(ar.ChainExcluded) as exc:
        ar.choose_mapping(exp_residues, af_by_pos, "P00000")
    assert exc.value.reason.startswith(ar.CHAIN_EXCL_FALLBACK_FAILED)


# --- score_mapping / coverage ---------------------------------------------


def test_score_mapping_computes_coverage_and_identity():
    mapped = [
        ar.MappedResidue(ar.Residue(1, "", 1, "ALA", "UNP", "P00000", 1), 1, 1, "ALA", True, ""),
        ar.MappedResidue(ar.Residue(2, "", 2, "GLY", "UNP", "P00000", 2), 2, 2, "SER", False, ""),
        ar.MappedResidue(ar.Residue(3, "", 3, "HIS", None, None, None), None, None, None, None, ar.REASON_NO_SIFTS_MAPPING),
    ]
    scores = ar.score_mapping(mapped)
    assert scores["n_observed"] == 3
    assert scores["n_kept"] == 2
    assert scores["coverage"] == pytest.approx(2 / 3)
    assert scores["n_compared"] == 2
    assert scores["n_mutations"] == 1
    assert scores["identity"] == pytest.approx(0.5)


# --- process_representative: end-to-end on fixture files ------------------


def _write_full_fixture(tmp_path, exp_rows, af_sequence, af_start, accession="P00000", with_sifts=True):
    exp_path = tmp_path / "exp.cif"
    af_path = tmp_path / "af.cif"
    write_experimental_cif(exp_path, exp_rows, with_sifts=with_sifts)
    write_alphafold_cif(af_path, af_sequence, start=af_start, accession=accession)
    return exp_path, af_path


def _rep(pdb_cif_path, alphafold_cif_path, af_start, af_end, accession="P00000"):
    return {
        "cluster_id": "TEST_A", "pdb_id": "TEST", "chain_id": "A", "uniprot_acc": accession,
        "pdb_cif_path": str(pdb_cif_path), "alphafold_cif_path": str(alphafold_cif_path),
        "alphafold_uniprot_start": af_start, "alphafold_uniprot_end": af_end,
        "eligible_homolog": True, "eligible_exact_train": True, "eligible_none": True,
    }


def test_process_representative_success_via_sifts(tmp_path):
    real = "MSSTLHSVFFTLKVS"
    exp_rows = [sifts_row(i, i, comp_id=AA3[c]) for i, c in enumerate(real, start=1)]
    exp_path, af_path = _write_full_fixture(tmp_path, exp_rows, real, af_start=1)

    rows, report = ar.process_representative(_rep(exp_path, af_path, 1, len(real)))
    assert report["status"] == "mapped"
    assert report["method"] == "sifts"
    assert report["coverage"] == pytest.approx(1.0)
    assert report["identity"] == pytest.approx(1.0)
    assert len(rows) == len(real)


def test_process_representative_success_via_fallback(tmp_path):
    real = "MSSTLHSVFFTLKVS"
    exp_rows = [unmapped_row(i, comp_id=AA3[c]) for i, c in enumerate(real, start=1)]
    exp_path, af_path = _write_full_fixture(tmp_path, exp_rows, real, af_start=1, with_sifts=True)

    rows, report = ar.process_representative(_rep(exp_path, af_path, 1, len(real)))
    assert report["status"] == "mapped"
    assert report["method"] == "fallback"
    assert report["identity"] == pytest.approx(1.0)


def test_process_representative_excluded_when_alphafold_numbering_invalid(tmp_path):
    real = "MSSTLHSVFFTLKVS"
    exp_rows = [sifts_row(i, i, comp_id=AA3[c]) for i, c in enumerate(real, start=1)]
    exp_path, af_path = _write_full_fixture(tmp_path, exp_rows, real, af_start=1, accession="P00000")

    # Claim a different accession than what's embedded in the AF fixture.
    rows, report = ar.process_representative(_rep(exp_path, af_path, 1, len(real), accession="Q99999"))
    assert report["status"] == "excluded"
    assert report["exclusion_reason"] == ar.CHAIN_EXCL_ALPHAFOLD_NOT_UNIPROT_NATIVE
    assert rows is None


def test_process_representative_excluded_below_coverage_threshold(tmp_path):
    real = "MSSTLHSVFFTLKVSAAAAA"  # 20 residues
    # Only the first 2 residues get a real SIFTS mapping; the rest are tags.
    exp_rows = [sifts_row(1, 1, comp_id=AA3[real[0]]), sifts_row(2, 2, comp_id=AA3[real[1]])]
    exp_rows += [unmapped_row(i, comp_id="HIS") for i in range(3, 21)]
    exp_path, af_path = _write_full_fixture(tmp_path, exp_rows, real, af_start=1)

    rows, report = ar.process_representative(_rep(exp_path, af_path, 1, len(real)))
    assert report["status"] == "excluded"
    assert report["exclusion_reason"].startswith(ar.CHAIN_EXCL_COVERAGE_BELOW_THRESHOLD)
    assert rows is None


def test_process_representative_keeps_engineered_mutation_chain(tmp_path):
    real = list("MSSTLHSVFFTLKVSAAAAA")  # 20 residues, all SIFTS-mapped
    exp_rows = [sifts_row(i, i, comp_id=AA3[c]) for i, c in enumerate(real, start=1)]
    af_sequence = list(real)
    af_sequence[9] = "G"  # AlphaFold differs at position 10 -> one engineered mutation
    exp_path, af_path = _write_full_fixture(tmp_path, exp_rows, "".join(af_sequence), af_start=1)

    rows, report = ar.process_representative(_rep(exp_path, af_path, 1, len(real)))
    assert report["status"] == "mapped"
    assert report["n_mutations"] == 1
    assert report["identity"] == pytest.approx(19 / 20)
    mutated_row = next(r for r in rows if r["auth_seq_id"] == 10)
    assert mutated_row["match_flag"] is False
    assert mutated_row["unmapped_reason"] == ""  # kept, not excluded


# --- write_parquet round-trip / schema -------------------------------------


def test_write_parquet_schema_round_trip(tmp_path):
    rows = [
        {
            "auth_seq_id": 1, "auth_ins_code": "", "label_seq_id": 1, "residue_name": "ALA",
            "is_observed": True, "uniprot_acc": "P00000", "uniprot_resnum": 1, "alphafold_resnum": 1,
            "alphafold_residue_name": "ALA", "match_flag": True, "method": "sifts", "unmapped_reason": "",
        }
    ]
    out_path = tmp_path / "TEST_A.parquet"
    ar.write_parquet(rows, out_path)
    df = pd.read_parquet(out_path)
    assert list(df.columns) == ar.PARQUET_FIELDS
    assert df.iloc[0]["auth_seq_id"] == 1


# --- attrition / coverage sweep --------------------------------------------


def test_compute_attrition_per_leakage_mode():
    report_rows = [
        {"eligible_homolog": True, "eligible_exact_train": True, "eligible_none": True, "status": "mapped"},
        {"eligible_homolog": False, "eligible_exact_train": True, "eligible_none": True, "status": "excluded"},
    ]
    attrition = ar.compute_attrition(report_rows)
    by_mode = {r["leakage_filter_mode"]: r for r in attrition}
    assert by_mode["homolog"]["n_eligible_representatives"] == 1
    assert by_mode["homolog"]["n_mapped"] == 1
    assert by_mode["exact_train"]["n_eligible_representatives"] == 2
    assert by_mode["exact_train"]["n_mapped"] == 1
    assert by_mode["exact_train"]["n_excluded"] == 1


def test_compute_coverage_sweep_counts_mapped_chains_at_each_threshold():
    report_rows = [
        {"status": "mapped", "coverage": 0.95},
        {"status": "mapped", "coverage": 0.85},
        {"status": "mapped", "coverage": 0.55},
        {"status": "excluded", "coverage": 0.0},
    ]
    sweep = ar.compute_coverage_sweep(report_rows)
    by_threshold = {r["coverage_threshold"]: r["n_chains_retained"] for r in sweep}
    assert by_threshold[0.50] == 3
    assert by_threshold[0.80] == 2
    assert by_threshold[0.90] == 1
    assert by_threshold[0.95] == 1


# --- priority ordering -----------------------------------------------------


def test_sort_by_priority_puts_primary_set_first():
    rows = [
        {"eligible_homolog": False, "eligible_exact_train": False, "cluster_id": "3CCC_A"},
        {"eligible_homolog": True, "eligible_exact_train": True, "cluster_id": "2BBB_A"},
        {"eligible_homolog": True, "eligible_exact_train": True, "cluster_id": "1AAA_A"},
    ]
    ordered = ar.sort_by_priority(rows)
    assert [r["cluster_id"] for r in ordered] == ["1AAA_A", "2BBB_A", "3CCC_A"]
