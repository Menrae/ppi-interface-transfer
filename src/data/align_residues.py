"""Phase 4: residue-level mapping between experimental PDB chains and
AlphaFold models.

The most correctness-critical step in the project -- a mapping error
silently corrupts every downstream metric. This module prefers failing
loudly and excluding a chain over guessing: every residue that isn't kept
carries an explicit reason, and every chain-level exclusion does too.

Primary mapping source: the per-residue SIFTS cross-references PDBe embeds
directly in ``_atom_site`` rows of its "updated" mmCIF (already cached by
Phase 3, see PLAN.md Phase 3/4) -- ``pdbx_sifts_xref_db_name/_acc/_num``.
The AlphaFold mmCIF (fetched from AlphaFold DB via the EBI/PDBe mirror)
carries the *same* embedded columns, so the same parsing/reduction code
handles both files; this is also what lets AlphaFold's own claim that its
numbering is UniProt-native be verified per chain, not assumed.

Runnable as: .venv/bin/python -m src.data.align_residues
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
from dataclasses import dataclass
from pathlib import Path

import gemmi
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

from src import config
from src.data.fetch_structures import ALPHAFOLD_API_CACHE_DIR

logger = logging.getLogger(__name__)

# --- Paths ---------------------------------------------------------------

FETCH_REPORT_PATH = config.INTERIM_DATA_DIR / "fetch_report.csv"
RESIDUE_MAPS_DIR = config.INTERIM_DATA_DIR / "residue_mappings"
MAPPING_REPORT_PATH = config.INTERIM_DATA_DIR / "mapping_report.csv"
ATTRITION_PATH = config.INTERIM_DATA_DIR / "phase4_attrition.csv"
COVERAGE_SWEEP_PATH = config.INTERIM_DATA_DIR / "phase4_coverage_sweep.csv"

# Thresholds swept for the "how many chains would different coverage
# cutoffs retain" report (PLAN.md Phase 4 item 6) -- purely a reporting
# sweep, does not change MIN_MAPPED_COVERAGE itself.
COVERAGE_SWEEP_THRESHOLDS = (0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# --- Per-residue and per-chain exclusion reasons (logged distinctly) -----

REASON_NO_SIFTS_MAPPING = "no_sifts_mapping"
REASON_SIFTS_WRONG_ACCESSION = "sifts_wrong_accession"
REASON_SIFTS_NUM_NOT_INTEGER = "sifts_num_not_integer"
REASON_NON_AMINO_ACID = "non_amino_acid_residue"
REASON_OUTSIDE_AF_RANGE = "outside_alphafold_model_range"
REASON_FALLBACK_NO_PARTNER = "fallback_no_alignment_partner"

CHAIN_EXCL_ALPHAFOLD_RANGE_MISMATCH = "alphafold_range_mismatch"
CHAIN_EXCL_ALPHAFOLD_NOT_UNIPROT_NATIVE = "alphafold_numbering_not_uniprot_native"
CHAIN_EXCL_MULTIPLE_MODELS = "multiple_models_unexpected"
CHAIN_EXCL_MULTICHAIN_ALPHAFOLD = "multichain_alphafold_model_unexpected"
CHAIN_EXCL_CHAIN_NOT_FOUND = "chain_not_found_in_mmcif"
CHAIN_EXCL_FALLBACK_FAILED = "fallback_failed_validation"
CHAIN_EXCL_COVERAGE_BELOW_THRESHOLD = "coverage_below_threshold"
CHAIN_EXCL_READ_ERROR = "structure_read_error"
CHAIN_EXCL_NOT_OFFICIAL_ALPHAFOLD = "alphafold_model_not_official_afdb"

# AlphaFold DB's own "providerId" for models computed by Google DeepMind
# with the standard AlphaFold2 pipeline. AlphaFold DB also hosts
# third-party "Community" submissions (e.g. ColabFold, at various
# versions) under other provider ids for the *same* accession, returned by
# the same prediction API Phase 3 queried -- these do not carry the
# training-cutoff/pipeline guarantees this project's AF2 comparison
# assumes (PLAN.md confound (a)), so they must never be silently accepted
# as "the AlphaFold model" just because they were the API's first/only hit
# for that accession.
OFFICIAL_ALPHAFOLD_PROVIDER_ID = "GDM"

PARQUET_FIELDS = [
    "auth_seq_id",
    "auth_ins_code",
    "label_seq_id",
    "residue_name",
    "is_observed",
    "uniprot_acc",
    "uniprot_resnum",
    "alphafold_resnum",
    "alphafold_residue_name",
    "match_flag",
    "method",
    "unmapped_reason",
]

MAPPING_REPORT_FIELDS = [
    "cluster_id",
    "pdb_id",
    "chain_id",
    "uniprot_acc",
    "eligible_homolog",
    "eligible_exact_train",
    "eligible_none",
    "method",
    "n_observed",
    "n_kept",
    "coverage",
    "n_compared",
    "n_mutations",
    "identity",
    "status",
    "exclusion_reason",
]

ATTRITION_FIELDS = [
    "leakage_filter_mode",
    "n_eligible_representatives",
    "n_mapped",
    "n_excluded",
    "note",
]

_ATOM_SITE_BASE_TAGS = [
    "_atom_site.group_PDB",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.label_seq_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.auth_comp_id",
    "_atom_site.occupancy",
    "_atom_site.pdbx_PDB_model_num",
]
_ATOM_SITE_SIFTS_TAGS = [
    "_atom_site.pdbx_sifts_xref_db_name",
    "_atom_site.pdbx_sifts_xref_db_acc",
    "_atom_site.pdbx_sifts_xref_db_num",
]


# --- Data model ------------------------------------------------------------


@dataclass(frozen=True)
class Residue:
    auth_seq_id: int
    ins_code: str
    label_seq_id: int
    comp_id: str
    sifts_db_name: str | None
    sifts_db_acc: str | None
    sifts_db_num: int | None


@dataclass
class MappedResidue:
    residue: Residue
    uniprot_resnum: int | None
    alphafold_resnum: int | None
    alphafold_residue_name: str | None
    match_flag: bool | None
    unmapped_reason: str


class ChainExcluded(Exception):
    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


# --- mmCIF parsing (shared by experimental and AlphaFold files) -----------
# The two files share the exact same _atom_site schema, including the
# embedded pdbx_sifts_xref_db_* columns -- one parser, used for both.


def _clean(value: str | None) -> str | None:
    if value is None or value in ("?", "."):
        return None
    return value


def _parse_int(value: str | None) -> int | None:
    value = _clean(value)
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def load_raw_atom_site_rows(cif_path: Path, auth_asym_id: str) -> list[dict]:
    """Return one dict per _atom_site row for auth_asym_id, all fields as
    raw strings (mmCIF '?'/'.' left as-is -- callers use _clean/_parse_int).

    Reads directly from the gzip-compressed cache file (gemmi handles the
    decompression). Missing sifts columns (absent from the file entirely,
    as opposed to present-but-'?') are treated as all-None, never crash.
    """
    doc = gemmi.cif.read(str(cif_path))
    block = doc.sole_block()

    has_sifts = all(len(block.find_loop(tag)) > 0 for tag in _ATOM_SITE_SIFTS_TAGS)
    tags = _ATOM_SITE_BASE_TAGS + (_ATOM_SITE_SIFTS_TAGS if has_sifts else [])
    table = block.find(tags)

    rows = []
    for row in table:
        values = {tag.split(".", 1)[1]: row[i] for i, tag in enumerate(tags)}
        if values["auth_asym_id"] != auth_asym_id:
            continue
        if _clean(values["group_PDB"]) is None:
            continue
        rows.append(values)
    return rows


def reduce_to_residues(rows: list[dict]) -> tuple[list[Residue], set[str]]:
    """Collapse raw atom rows to one Residue per (auth_seq_id, ins_code).

    Only rows with a real label_seq_id (i.e. part of a polymer entity, not
    a water/ligand that happens to share the chain's auth_asym_id -- see
    PLAN.md Phase 4) are considered. Alternate locations are resolved
    deterministically: highest occupancy wins, ties broken by ascending
    label_alt_id -- documented here, not silently arbitrary.

    Returns (residues sorted by (auth_seq_id, ins_code), set of distinct
    pdbx_PDB_model_num values seen) -- callers check the latter is a
    singleton rather than assuming a single-model file.
    """
    groups: dict[tuple[int, str], list[dict]] = {}
    model_nums: set[str] = set()
    for row in rows:
        model_nums.add(row["pdbx_PDB_model_num"])
        if _clean(row["label_seq_id"]) is None:
            continue
        auth_seq_id = _parse_int(row["auth_seq_id"])
        if auth_seq_id is None:
            continue
        ins_code = _clean(row["pdbx_PDB_ins_code"]) or ""
        groups.setdefault((auth_seq_id, ins_code), []).append(row)

    residues = []
    for (auth_seq_id, ins_code), group_rows in groups.items():

        def sort_key(r: dict) -> tuple[float, str]:
            occ = r["occupancy"]
            occupancy = float(occ) if _clean(occ) is not None else 1.0
            alt_id = _clean(r["label_alt_id"]) or ""
            return (-occupancy, alt_id)

        best = sorted(group_rows, key=sort_key)[0]
        residues.append(
            Residue(
                auth_seq_id=auth_seq_id,
                ins_code=ins_code,
                label_seq_id=_parse_int(best["label_seq_id"]),
                comp_id=best["auth_comp_id"],
                sifts_db_name=_clean(best.get("pdbx_sifts_xref_db_name")),
                sifts_db_acc=_clean(best.get("pdbx_sifts_xref_db_acc")),
                sifts_db_num=_parse_int(best.get("pdbx_sifts_xref_db_num")),
            )
        )
    residues.sort(key=lambda r: (r.auth_seq_id, r.ins_code))
    return residues, model_nums


def load_chain_residues(cif_path: Path, auth_asym_id: str) -> list[Residue]:
    """Load + reduce one chain's observed residues; raises ChainExcluded on
    any structural surprise (never silently guesses)."""
    try:
        rows = load_raw_atom_site_rows(cif_path, auth_asym_id)
    except Exception as exc:  # noqa: BLE001 - any parse failure excludes, not crashes the batch
        raise ChainExcluded(f"{CHAIN_EXCL_READ_ERROR}:{exc}") from exc

    if not rows:
        raise ChainExcluded(CHAIN_EXCL_CHAIN_NOT_FOUND)

    residues, model_nums = reduce_to_residues(rows)
    if len(model_nums) > 1:
        raise ChainExcluded(f"{CHAIN_EXCL_MULTIPLE_MODELS}:{sorted(model_nums)}")
    if not residues:
        raise ChainExcluded(CHAIN_EXCL_CHAIN_NOT_FOUND)
    return residues


def load_alphafold_residues(cif_path: Path) -> list[Residue]:
    """Load the AlphaFold model's residues, verifying it really is the
    single-chain, single-model file we assume (never guessed)."""
    try:
        doc = gemmi.cif.read(str(cif_path))
        block = doc.sole_block()
        chain_col = block.find_loop("_atom_site.auth_asym_id")
        chains = sorted({c for c in chain_col})
    except Exception as exc:  # noqa: BLE001
        raise ChainExcluded(f"{CHAIN_EXCL_READ_ERROR}:{exc}") from exc

    if len(chains) != 1:
        raise ChainExcluded(f"{CHAIN_EXCL_MULTICHAIN_ALPHAFOLD}:{chains}")
    return load_chain_residues(cif_path, chains[0])


# --- Residue identity (parent-residue resolution for modified residues) --


def one_letter_code(comp_id: str) -> str | None:
    """Parent amino-acid one-letter code for comp_id, or None if comp_id
    isn't an amino acid at all (e.g. a non-standard linker/crosslinker).

    Uses gemmi's tabulated chemical-component data, which already resolves
    modified residues (MSE, SEP, KCX, ...) to their parent identity -- e.g.
    MSE (selenomethionine) -> 'M', matching MET. This is how modified
    residues are "treated through their parent residue" per PLAN.md Phase 4.
    """
    info = gemmi.find_tabulated_residue(comp_id)
    if info is None or not info.is_amino_acid():
        return None
    return info.one_letter_code.upper()


def verify_official_alphafold_provider(uniprot_acc: str) -> None:
    """Raise ChainExcluded unless AlphaFold DB's own cached prediction-API
    response (Phase 3) says this accession's model was computed by Google
    DeepMind's standard pipeline (providerId "GDM"), not a third-party
    Community submission (e.g. ColabFold) that happens to also match the
    accession exactly. Checked directly from Phase 3's cached JSON --
    never inferred from the mmCIF alone, and never skipped."""
    path = ALPHAFOLD_API_CACHE_DIR / f"{uniprot_acc}.json"
    try:
        with open(path) as f:
            cached = json.load(f)
    except (OSError, ValueError) as exc:
        raise ChainExcluded(f"{CHAIN_EXCL_READ_ERROR}:{exc}") from exc

    body = cached.get("body") or []
    matches = [e for e in body if e.get("uniprotAccession") == uniprot_acc]
    provider_id = matches[0].get("providerId") if matches else None
    if not matches or provider_id != OFFICIAL_ALPHAFOLD_PROVIDER_ID:
        raise ChainExcluded(f"{CHAIN_EXCL_NOT_OFFICIAL_ALPHAFOLD}:providerId={provider_id}")


# --- AlphaFold numbering self-validation (never assumed) -------------------


def validate_alphafold_numbering(
    af_residues: list[Residue], expected_accession: str, start: int, end: int
) -> dict[int, Residue]:
    """Confirm the AlphaFold model's own numbering is UniProt-native over
    [start, end], using its own embedded SIFTS annotation -- never assumed
    from the recorded uniprotStart/uniprotEnd alone. Raises ChainExcluded
    on any disagreement. Returns {uniprot_position: Residue}."""
    by_seq_id = {r.auth_seq_id: r for r in af_residues}
    expected_range = set(range(start, end + 1))
    if set(by_seq_id) != expected_range:
        raise ChainExcluded(
            f"{CHAIN_EXCL_ALPHAFOLD_RANGE_MISMATCH}:expected={start}-{end},"
            f"got={min(by_seq_id) if by_seq_id else None}-{max(by_seq_id) if by_seq_id else None},"
            f"n={len(by_seq_id)}"
        )
    for auth_seq_id, residue in by_seq_id.items():
        if (
            residue.sifts_db_name != "UNP"
            or residue.sifts_db_acc != expected_accession
            or residue.sifts_db_num != auth_seq_id
        ):
            raise ChainExcluded(CHAIN_EXCL_ALPHAFOLD_NOT_UNIPROT_NATIVE)
    return by_seq_id


# --- Mapping methods ---------------------------------------------------


def _score_against_alphafold(
    uniprot_pos: int | None, residue: Residue, af_by_pos: dict[int, Residue], reason_if_unmapped: str
) -> MappedResidue:
    if uniprot_pos is None:
        return MappedResidue(residue, None, None, None, None, reason_if_unmapped)

    if one_letter_code(residue.comp_id) is None:
        return MappedResidue(residue, uniprot_pos, None, None, None, REASON_NON_AMINO_ACID)

    af_residue = af_by_pos.get(uniprot_pos)
    if af_residue is None:
        return MappedResidue(residue, uniprot_pos, None, None, None, REASON_OUTSIDE_AF_RANGE)

    exp_letter = one_letter_code(residue.comp_id)
    af_letter = one_letter_code(af_residue.comp_id)
    match = exp_letter is not None and af_letter is not None and exp_letter == af_letter
    return MappedResidue(residue, uniprot_pos, uniprot_pos, af_residue.comp_id, match, "")


def map_via_sifts(
    exp_residues: list[Residue], af_by_pos: dict[int, Residue], expected_accession: str
) -> list[MappedResidue]:
    mapped = []
    for residue in exp_residues:
        if residue.sifts_db_name != "UNP":
            mapped.append(_score_against_alphafold(None, residue, af_by_pos, REASON_NO_SIFTS_MAPPING))
        elif residue.sifts_db_acc != expected_accession:
            mapped.append(_score_against_alphafold(None, residue, af_by_pos, REASON_SIFTS_WRONG_ACCESSION))
        elif residue.sifts_db_num is None:
            mapped.append(_score_against_alphafold(None, residue, af_by_pos, REASON_SIFTS_NUM_NOT_INTEGER))
        else:
            mapped.append(_score_against_alphafold(residue.sifts_db_num, residue, af_by_pos, ""))
    return mapped


def _build_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    # Free end gaps on both sides -- lets an expression tag/linker at
    # either terminus align to nothing instead of forcing a bad-scoring
    # internal alignment, without weakening internal gap penalties.
    aligner.end_insertion_score = 0.0
    aligner.end_deletion_score = 0.0
    return aligner


def map_via_fallback(
    exp_residues: list[Residue], af_by_pos: dict[int, Residue]
) -> list[MappedResidue]:
    """Global (free-end-gap) alignment of the observed experimental
    sequence against the AlphaFold model's own sequence, used when SIFTS
    residue mapping is absent or fails validation (PLAN.md Phase 4)."""
    exp_index_map: list[int] = []
    exp_seq_chars: list[str] = []
    for i, residue in enumerate(exp_residues):
        letter = one_letter_code(residue.comp_id)
        if letter is None:
            continue
        exp_index_map.append(i)
        exp_seq_chars.append(letter)
    exp_seq = "".join(exp_seq_chars)

    af_positions = sorted(af_by_pos)
    af_seq = "".join(one_letter_code(af_by_pos[p].comp_id) or "X" for p in af_positions)

    mapped = [
        _score_against_alphafold(None, residue, af_by_pos, REASON_FALLBACK_NO_PARTNER)
        for residue in exp_residues
    ]
    # Non-amino-acid residues were excluded from the alignment sequence
    # entirely (can't be aligned as themselves); classify them properly.
    aligned_indices = set(exp_index_map)
    for i, residue in enumerate(exp_residues):
        if i not in aligned_indices:
            mapped[i] = _score_against_alphafold(None, residue, af_by_pos, REASON_NON_AMINO_ACID)

    if not exp_seq or not af_seq:
        return mapped

    aligner = _build_aligner()
    alignment = next(iter(aligner.align(exp_seq, af_seq)))
    target_blocks, query_blocks = alignment.aligned
    for (t_start, t_end), (q_start, q_end) in zip(target_blocks, query_blocks):
        for offset in range(t_end - t_start):
            exp_seq_idx = t_start + offset
            af_seq_idx = q_start + offset
            original_idx = exp_index_map[exp_seq_idx]
            uniprot_pos = af_positions[af_seq_idx]
            residue = exp_residues[original_idx]
            mapped[original_idx] = _score_against_alphafold(uniprot_pos, residue, af_by_pos, "")
    return mapped


# --- Per-chain scoring / method selection --------------------------------


def score_mapping(mapped: list[MappedResidue]) -> dict:
    n_observed = len(mapped)
    n_compared = sum(1 for m in mapped if m.alphafold_resnum is not None)
    n_matched = sum(1 for m in mapped if m.match_flag is True)
    return {
        "n_observed": n_observed,
        "n_kept": n_compared,
        "coverage": (n_compared / n_observed) if n_observed else 0.0,
        "n_compared": n_compared,
        "n_mutations": n_compared - n_matched,
        "identity": (n_matched / n_compared) if n_compared else 0.0,
    }


def choose_mapping(
    exp_residues: list[Residue], af_by_pos: dict[int, Residue], expected_accession: str
) -> tuple[list[MappedResidue], str]:
    """Pick SIFTS or fallback per PLAN.md Phase 4 item 4. Raises
    ChainExcluded if neither method validates."""
    sifts_mapped = map_via_sifts(exp_residues, af_by_pos, expected_accession)
    sifts_score = score_mapping(sifts_mapped)

    if sifts_score["n_compared"] > 0 and sifts_score["identity"] >= config.RESIDUE_MAPPING_VALIDATION_IDENTITY:
        return sifts_mapped, "sifts"

    fallback_mapped = map_via_fallback(exp_residues, af_by_pos)
    fallback_score = score_mapping(fallback_mapped)
    if fallback_score["n_compared"] > 0 and fallback_score["identity"] >= config.RESIDUE_MAPPING_VALIDATION_IDENTITY:
        return fallback_mapped, "fallback"

    raise ChainExcluded(
        f"{CHAIN_EXCL_FALLBACK_FAILED}:sifts_n_compared={sifts_score['n_compared']},"
        f"sifts_identity={sifts_score['identity']:.2f},"
        f"fallback_n_compared={fallback_score['n_compared']},"
        f"fallback_identity={fallback_score['identity']:.2f}"
    )


def build_parquet_rows(mapped: list[MappedResidue], uniprot_acc: str, method: str) -> list[dict]:
    rows = []
    for m in mapped:
        r = m.residue
        rows.append(
            {
                "auth_seq_id": r.auth_seq_id,
                "auth_ins_code": r.ins_code,
                "label_seq_id": r.label_seq_id,
                "residue_name": r.comp_id,
                "is_observed": True,
                "uniprot_acc": uniprot_acc,
                "uniprot_resnum": m.uniprot_resnum,
                "alphafold_resnum": m.alphafold_resnum,
                "alphafold_residue_name": m.alphafold_residue_name,
                "match_flag": m.match_flag,
                "method": method,
                "unmapped_reason": m.unmapped_reason,
            }
        )
    return rows


# --- One representative, start to finish --------------------------------


def process_representative(rep: dict) -> tuple[list[dict] | None, dict]:
    """Returns (parquet_rows_or_None, mapping_report_row)."""
    base_report = {
        "cluster_id": rep["cluster_id"],
        "pdb_id": rep["pdb_id"],
        "chain_id": rep["chain_id"],
        "uniprot_acc": rep["uniprot_acc"],
        "eligible_homolog": rep["eligible_homolog"],
        "eligible_exact_train": rep["eligible_exact_train"],
        "eligible_none": rep["eligible_none"],
        "method": "",
        "n_observed": 0,
        "n_kept": 0,
        "coverage": 0.0,
        "n_compared": 0,
        "n_mutations": 0,
        "identity": 0.0,
        "status": "excluded",
        "exclusion_reason": "",
    }
    try:
        verify_official_alphafold_provider(rep["uniprot_acc"])
        af_residues = load_alphafold_residues(Path(rep["alphafold_cif_path"]))
        af_by_pos = validate_alphafold_numbering(
            af_residues, rep["uniprot_acc"], rep["alphafold_uniprot_start"], rep["alphafold_uniprot_end"]
        )
        exp_residues = load_chain_residues(Path(rep["pdb_cif_path"]), rep["chain_id"])
        mapped, method = choose_mapping(exp_residues, af_by_pos, rep["uniprot_acc"])
        scores = score_mapping(mapped)
        if scores["coverage"] < config.MIN_MAPPED_COVERAGE:
            raise ChainExcluded(f"{CHAIN_EXCL_COVERAGE_BELOW_THRESHOLD}:{scores['coverage']:.2f}")
    except ChainExcluded as exc:
        base_report["exclusion_reason"] = exc.reason
        logger.warning("excluded %s_%s: %s", rep["pdb_id"], rep["chain_id"], exc.reason)
        return None, base_report

    base_report.update(
        {
            "method": method,
            "n_observed": scores["n_observed"],
            "n_kept": scores["n_kept"],
            "coverage": scores["coverage"],
            "n_compared": scores["n_compared"],
            "n_mutations": scores["n_mutations"],
            "identity": scores["identity"],
            "status": "mapped",
        }
    )
    return build_parquet_rows(mapped, rep["uniprot_acc"], method), base_report


# --- Loading fetch_report.csv / priority ordering -------------------------


def _csv_bool(value: str) -> bool:
    return value == "True"


def load_representatives(path: Path = FETCH_REPORT_PATH) -> list[dict]:
    """Load fetch_report.csv rows with both structures downloaded
    (overall_status == "complete"), the scope for Phase 4."""
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["overall_status"] != "complete":
                continue
            rows.append(
                {
                    "cluster_id": row["cluster_id"],
                    "pdb_id": row["pdb_id"],
                    "chain_id": row["chain_id"],
                    "uniprot_acc": row["uniprot_acc"],
                    "pdb_cif_path": row["pdb_cif_path"],
                    "alphafold_cif_path": row["alphafold_cif_path"],
                    "alphafold_uniprot_start": int(row["alphafold_uniprot_start"]),
                    "alphafold_uniprot_end": int(row["alphafold_uniprot_end"]),
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


# --- Output writers / attrition ------------------------------------------


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
        n_mapped = sum(1 for r in subset if r["status"] == "mapped")
        attrition.append(
            {
                "leakage_filter_mode": mode,
                "n_eligible_representatives": len(subset),
                "n_mapped": n_mapped,
                "n_excluded": len(subset) - n_mapped,
                "note": f"meets_target({config.MIN_TEST_CHAINS_TARGET})="
                f"{n_mapped >= config.MIN_TEST_CHAINS_TARGET}",
            }
        )
    return attrition


def compute_coverage_sweep(report_rows: list[dict]) -> list[dict]:
    """How many *mapped* chains would survive at each coverage threshold
    (PLAN.md Phase 4 item 6) -- reporting only, does not change
    MIN_MAPPED_COVERAGE or which chains are already excluded by it."""
    mapped_rows = [r for r in report_rows if r["status"] == "mapped"]
    sweep = []
    for threshold in COVERAGE_SWEEP_THRESHOLDS:
        n_survive = sum(1 for r in mapped_rows if r["coverage"] >= threshold)
        sweep.append(
            {
                "coverage_threshold": threshold,
                "n_chains_retained": n_survive,
                "note": f"of {len(mapped_rows)} chains mapped at MIN_MAPPED_COVERAGE="
                f"{config.MIN_MAPPED_COVERAGE}",
            }
        )
    return sweep


# --- Orchestration -------------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "align_residues.log"

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
    rows = sort_by_priority(load_representatives())
    logger.info(
        "Loaded %d representatives with both structures downloaded (fetch_report.csv "
        "overall_status == 'complete')",
        len(rows),
    )

    report_rows = []
    for i, rep in enumerate(rows):
        parquet_rows, report_row = process_representative(rep)
        report_rows.append(report_row)
        if parquet_rows is not None:
            out_path = RESIDUE_MAPS_DIR / f"{rep['pdb_id']}_{rep['chain_id']}.parquet"
            write_parquet(parquet_rows, out_path)
        if (i + 1) % 250 == 0 or (i + 1) == len(rows):
            n_mapped_so_far = sum(1 for r in report_rows if r["status"] == "mapped")
            logger.info("Processed %d/%d representatives (%d mapped)", i + 1, len(rows), n_mapped_so_far)

    write_csv(report_rows, MAPPING_REPORT_FIELDS, MAPPING_REPORT_PATH)
    logger.info("Wrote mapping report (%d rows) to %s", len(report_rows), MAPPING_REPORT_PATH)

    attrition = compute_attrition(report_rows)
    write_csv(attrition, ATTRITION_FIELDS, ATTRITION_PATH)
    logger.info("Wrote per-leakage-mode attrition to %s", ATTRITION_PATH)
    for row in attrition:
        logger.info(
            "mode=%-14s n_eligible=%-6d n_mapped=%-6d %s",
            row["leakage_filter_mode"],
            row["n_eligible_representatives"],
            row["n_mapped"],
            row["note"],
        )

    sweep = compute_coverage_sweep(report_rows)
    write_csv(sweep, ["coverage_threshold", "n_chains_retained", "note"], COVERAGE_SWEEP_PATH)
    logger.info("Wrote coverage-threshold sweep to %s", COVERAGE_SWEEP_PATH)

    exclusion_counts: dict[str, int] = {}
    for r in report_rows:
        if r["exclusion_reason"]:
            key = r["exclusion_reason"].split(":", 1)[0]
            exclusion_counts[key] = exclusion_counts.get(key, 0) + 1
    logger.info("Exclusion-reason breakdown: %s", exclusion_counts)

    method_counts: dict[str, int] = {}
    for r in report_rows:
        if r["method"]:
            method_counts[r["method"]] = method_counts.get(r["method"], 0) + 1
    logger.info("Mapping-method breakdown: %s", method_counts)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Phase 4: residue-numbering mapping")
    parser.parse_args(argv)
    run()


if __name__ == "__main__":
    main()
