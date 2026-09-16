"""Phase 2: redundancy reduction + PeSTo-overlap (leakage) flagging.

Clusters Phase 1's candidate chains with MMseqs2 at
``config.SEQUENCE_IDENTITY_CUTOFF`` / ``config.CLUSTER_MIN_COVERAGE``, picks
one representative per cluster deterministically, and searches those
representatives against the sequences of every chain in PeSTo's published
train/test/validation split files (see PLAN.md Phase 2 for why all three
files are treated as excluded, not just train+validation) to flag sequence-
level overlap with PeSTo's own training/model-selection data.

Runnable as: .venv/bin/python -m src.data.cluster_and_split
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import subprocess
from pathlib import Path

from src import config
from src.data.select_complexes import build_session

logger = logging.getLogger(__name__)

# --- Paths -----------------------------------------------------------------

CANDIDATES_PATH = config.INTERIM_DATA_DIR / "candidates.csv"
POLYMER_ENTITY_CACHE_DIR = config.RAW_DATA_DIR / "rcsb" / "polymer_entities"

PDB_SEQRES_URL = "https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
PDB_SEQRES_DIR = config.RAW_DATA_DIR / "pdb_seqres"
PDB_SEQRES_PATH = PDB_SEQRES_DIR / "pdb_seqres.txt.gz"

PESTO_SPLITS_DIR = config.RAW_DATA_DIR / "pesto_splits"
PESTO_SPLIT_FILES = {
    "train": PESTO_SPLITS_DIR / "subunits_train_set.txt",
    "test": PESTO_SPLITS_DIR / "subunits_test_set.txt",
    "validation": PESTO_SPLITS_DIR / "subunits_validation_set.txt",
}
PESTO_ALL_SPLITS_FASTA = PESTO_SPLITS_DIR / "pesto_all_splits.fasta"

MMSEQS_BIN = config.EXTERNAL_DIR / "mmseqs" / "bin" / "mmseqs"
MMSEQS_WORK_DIR = config.RAW_DATA_DIR / "mmseqs_work"
CLUSTER_TMP_DIR = MMSEQS_WORK_DIR / "cluster_tmp"
CLUSTER_OUT_PREFIX = MMSEQS_WORK_DIR / "candidates_cluster"
SEARCH_TMP_DIR = MMSEQS_WORK_DIR / "search_tmp"
SEARCH_OUT_PATH = MMSEQS_WORK_DIR / "pesto_homology_search_raw.tsv"
REPRESENTATIVES_FASTA = MMSEQS_WORK_DIR / "candidate_representatives.fasta"

SEQUENCES_FASTA_PATH = config.INTERIM_DATA_DIR / "sequences.fasta"
CLUSTERS_PATH = config.INTERIM_DATA_DIR / "clusters.tsv"
CANDIDATES_DEDUP_PATH = config.INTERIM_DATA_DIR / "candidates_dedup.csv"
ATTRITION_PATH = config.INTERIM_DATA_DIR / "phase2_attrition.csv"
LEAKAGE_SWEEP_PATH = config.INTERIM_DATA_DIR / "leakage_threshold_sweep.csv"
PESTO_HOMOLOGY_SEARCH_PATH = config.INTERIM_DATA_DIR / "pesto_homology_search.tsv"

CANDIDATE_FIELDS = [
    "pdb_id",
    "release_date",
    "resolution",
    "n_protein_entities",
    "entity_id",
    "chain_id",
    "seq_length",
    "polymer_type",
    "uniprot_ids",
    "n_uniprot_ids",
    "sifts_uniprot_ids",
    "sifts_agrees",
]

DEDUP_FIELDS = CANDIDATE_FIELDS + [
    "cluster_id",
    "cluster_size",
    "max_identity_train",
    "best_hit_train_id",
    "max_identity_test",
    "best_hit_test_id",
    "max_identity_validation",
    "best_hit_validation_id",
    "max_identity_any",
    "pesto_homolog_overlap",
    "pesto_exact_train_overlap",
]

CLUSTERS_FIELDS = ["cluster_id", "pdb_id", "chain_id", "is_representative"]
ATTRITION_FIELDS = ["stage", "n_before", "n_dropped", "n_remaining", "note"]
SWEEP_FIELDS = ["mode", "identity_threshold", "n_survive", "meets_target", "meets_floor"]
HOMOLOGY_FIELDS = [
    "rep_id",
    "max_identity_train",
    "best_hit_train_id",
    "max_identity_test",
    "best_hit_test_id",
    "max_identity_validation",
    "best_hit_validation_id",
    "max_identity_any",
    "pesto_homolog_overlap",
    "pesto_exact_train_overlap",
]

SEARCH_FORMAT_COLUMNS = ["query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"]


def chain_id_str(pdb_id: str, chain_id: str) -> str:
    return f"{pdb_id.upper()}_{chain_id}"


# --- Candidate sequence loading (reuses Phase 1's cached entity JSON) ------


def load_candidates() -> list[dict]:
    with open(CANDIDATES_PATH, newline="") as f:
        return list(csv.DictReader(f))


def load_candidate_sequences(rows: list[dict]) -> tuple[dict[str, str], list[dict]]:
    """Return ({id -> sequence}, dropped_rows) for candidate chains.

    Sequences are read from the polymer-entity JSON already cached by Phase 1
    (`entity_poly.pdbx_seq_one_letter_code_can`) -- no network calls. A row
    is dropped (and logged) only if its cached entity JSON is missing or has
    no canonical sequence, which should not happen for a clean Phase 1 run.
    """
    sequences: dict[str, str] = {}
    dropped: list[dict] = []
    entity_cache: dict[tuple[str, str], dict] = {}

    for row in rows:
        pdb_id, entity_id, chain_id = row["pdb_id"], row["entity_id"], row["chain_id"]
        key = (pdb_id, entity_id)
        if key not in entity_cache:
            path = POLYMER_ENTITY_CACHE_DIR / f"{pdb_id}_{entity_id}.json"
            if path.exists():
                with open(path) as f:
                    entity_cache[key] = json.load(f)
            else:
                entity_cache[key] = {}

        entity_data = entity_cache[key]
        seq = (entity_data.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
        cid = chain_id_str(pdb_id, chain_id)
        if seq:
            sequences[cid] = seq
        else:
            logger.warning(
                "dropped %s: no cached canonical sequence in polymer-entity JSON", cid
            )
            dropped.append({**row, "drop_reason": "missing_cached_sequence"})

    return sequences, dropped


def write_fasta(records: dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for seq_id, seq in records.items():
            f.write(f">{seq_id}\n{seq}\n")


# --- PeSTo split sequences (bulk pdb_seqres.txt.gz lookup) ------------------


def download_pdb_seqres(session) -> Path:
    if PDB_SEQRES_PATH.exists():
        return PDB_SEQRES_PATH
    logger.info("Downloading %s to %s", PDB_SEQRES_URL, PDB_SEQRES_PATH)
    PDB_SEQRES_DIR.mkdir(parents=True, exist_ok=True)
    response = session.get(PDB_SEQRES_URL, stream=True)
    response.raise_for_status()
    tmp_path = PDB_SEQRES_PATH.with_suffix(".tmp")
    with open(tmp_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=1 << 20):
            f.write(chunk)
    tmp_path.rename(PDB_SEQRES_PATH)
    return PDB_SEQRES_PATH


def load_pesto_split_membership() -> dict[str, set[str]]:
    """Return {"train"/"test"/"validation": {"PDBID_CHAINID", ...}}."""
    membership: dict[str, set[str]] = {}
    for name, path in PESTO_SPLIT_FILES.items():
        if not path.exists():
            raise FileNotFoundError(
                f"PeSTo split file missing: {path} -- expected from "
                "github.com/LBM-EPFL/PeSTo/data/datasets/ (see PLAN.md Phase 2 "
                "risk: pesto_overlap must fail loudly, not silently skip flagging)"
            )
        with open(path) as f:
            # Normalize only the PDB ID (case-insensitive by convention), not
            # the chain ID (case-SENSITIVE -- e.g. chain "C" and chain "c" can
            # be genuinely distinct chains in the same entry). Uppercasing the
            # whole "PDBID_CHAINID" string would silently collapse those.
            entries = set()
            for line in f:
                line = line.strip()
                if not line:
                    continue
                pdb_id, chain = line.split("_", 1)
                entries.add(f"{pdb_id.upper()}_{chain}")
            membership[name] = entries
        logger.info("Loaded %d chains from PeSTo %s split (%s)", len(membership[name]), name, path.name)
    return membership


def lookup_pdb_seqres_sequences(
    seqres_path: Path, wanted_ids: set[str]
) -> tuple[dict[str, str], set[str]]:
    """Stream-parse pdb_seqres.txt.gz, returning ({id -> seq}, missing_ids).

    ``wanted_ids`` and the returned dict keys are "PDBID_CHAINID" with an
    uppercase PDB ID (matching chain_id_str and PeSTo's split-file naming).
    Only wanted entries are kept in memory -- the full file has ~1.16M
    records.
    """
    found: dict[str, str] = {}
    remaining = set(wanted_ids)
    with gzip.open(seqres_path, "rt") as f:
        header = None
        for line in f:
            if line.startswith(">"):
                pdb_chain = line[1:].split(None, 1)[0]
                if "_" not in pdb_chain:
                    header = None
                    continue
                pdb_id, chain = pdb_chain.split("_", 1)
                cid = f"{pdb_id.upper()}_{chain}"
                header = cid if cid in remaining else None
            elif header is not None:
                found[header] = line.strip()
                remaining.discard(header)
                header = None
    return found, remaining


# --- MMseqs2 invocation ------------------------------------------------------


def _run_mmseqs(args: list[str]) -> None:
    cmd = [str(MMSEQS_BIN)] + [str(a) for a in args]
    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("mmseqs stdout:\n%s", result.stdout)
        logger.error("mmseqs stderr:\n%s", result.stderr)
        raise RuntimeError(f"mmseqs command failed (exit {result.returncode}): {' '.join(cmd)}")
    logger.debug("mmseqs stdout:\n%s", result.stdout)


def run_mmseqs_cluster(fasta_path: Path, out_prefix: Path, tmp_dir: Path) -> Path:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    _run_mmseqs(
        [
            "easy-cluster",
            fasta_path,
            out_prefix,
            tmp_dir,
            "--min-seq-id",
            config.SEQUENCE_IDENTITY_CUTOFF,
            "-c",
            config.CLUSTER_MIN_COVERAGE,
            "--cov-mode",
            0,
        ]
    )
    cluster_tsv = Path(f"{out_prefix}_cluster.tsv")
    if not cluster_tsv.exists():
        raise RuntimeError(f"mmseqs easy-cluster did not produce {cluster_tsv}")
    return cluster_tsv


def parse_cluster_tsv(path: Path) -> dict[str, list[str]]:
    """Return {mmseqs_assigned_rep_id: [member_id, ...]} (rep included)."""
    groups: dict[str, list[str]] = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rep, member = line.split("\t")
            groups.setdefault(rep, []).append(member)
    return groups


def choose_representative(members: list[str], candidate_index: dict[str, dict]) -> str:
    """Deterministically pick one member per cluster.

    Criteria, in order: best (lowest) resolution, then longest chain
    (seq_length), then ascending id string (PDB ID, then chain ID) as a
    final, fully deterministic tiebreak.
    """

    def sort_key(member_id: str) -> tuple[float, int, str]:
        row = candidate_index[member_id]
        resolution = float(row["resolution"])
        seq_length = -int(row["seq_length"])  # longest first -> negate for ascending sort
        return (resolution, seq_length, member_id)

    return min(members, key=sort_key)


def run_mmseqs_search(
    query_fasta: Path, target_fasta: Path, out_path: Path, tmp_dir: Path
) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    min_identity = min(config.LEAKAGE_SWEEP_IDENTITY_THRESHOLDS)
    _run_mmseqs(
        [
            "easy-search",
            query_fasta,
            target_fasta,
            out_path,
            tmp_dir,
            "--min-seq-id",
            min_identity,
            "-c",
            config.CLUSTER_MIN_COVERAGE,
            "--cov-mode",
            0,
            "--format-mode",
            0,
            "--format-output",
            ",".join(SEARCH_FORMAT_COLUMNS),
        ]
    )
    return out_path


def parse_search_tsv(path: Path) -> dict[str, list[dict]]:
    """Return {query_id: [{"target":..., "fident": float, ...}, ...]}."""
    hits: dict[str, list[dict]] = {}
    if not path.exists():
        return hits
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            record = dict(zip(SEARCH_FORMAT_COLUMNS, fields))
            record["fident"] = float(record["fident"])
            hits.setdefault(record["query"], []).append(record)
    return hits


# --- Leakage flag computation ------------------------------------------------


def compute_leakage_flags(
    representative_ids: list[str],
    search_hits: dict[str, list[dict]],
    split_membership: dict[str, set[str]],
) -> dict[str, dict]:
    """Return {rep_id: {max_identity_train, best_hit_train_id, ..., pesto_homolog_overlap, pesto_exact_train_overlap}}."""
    flags: dict[str, dict] = {}
    for rep_id in representative_ids:
        result = {}
        best_any = 0.0
        for split_name in ("train", "test", "validation"):
            best_identity = 0.0
            best_hit = ""
            for hit in search_hits.get(rep_id, []):
                if hit["target"] in split_membership[split_name] and hit["fident"] > best_identity:
                    best_identity = hit["fident"]
                    best_hit = hit["target"]
            result[f"max_identity_{split_name}"] = best_identity
            result[f"best_hit_{split_name}_id"] = best_hit
            best_any = max(best_any, best_identity)

        result["max_identity_any"] = best_any
        result["pesto_homolog_overlap"] = best_any >= config.SEQUENCE_IDENTITY_CUTOFF
        result["pesto_exact_train_overlap"] = rep_id in split_membership["train"]
        flags[rep_id] = result
    return flags


def sweep_thresholds(flags: dict[str, dict]) -> list[dict]:
    """Report survivor counts at each swept identity threshold + exact-ID-only.

    A representative "survives" a threshold if its max_identity_any is
    strictly below that threshold (i.e. not flagged as overlap at that
    cutoff); it survives the exact-ID mode if pesto_exact_train_overlap is
    False. This never loosens SEQUENCE_IDENTITY_CUTOFF for the primary
    analysis -- it only reports how the count *would* change, per PLAN.md
    confound (e).
    """
    rows: list[dict] = []
    for threshold in config.LEAKAGE_SWEEP_IDENTITY_THRESHOLDS:
        n_survive = sum(1 for f in flags.values() if f["max_identity_any"] < threshold)
        rows.append(
            {
                "mode": f"identity<{threshold:.2f}",
                "identity_threshold": threshold,
                "n_survive": n_survive,
                "meets_target": n_survive >= config.MIN_TEST_CHAINS_TARGET,
                "meets_floor": n_survive >= config.MIN_TEST_CHAINS_FLOOR,
            }
        )
    n_survive_exact = sum(1 for f in flags.values() if not f["pesto_exact_train_overlap"])
    rows.append(
        {
            "mode": "exact_id_train_only",
            "identity_threshold": "",
            "n_survive": n_survive_exact,
            "meets_target": n_survive_exact >= config.MIN_TEST_CHAINS_TARGET,
            "meets_floor": n_survive_exact >= config.MIN_TEST_CHAINS_FLOOR,
        }
    )
    return rows


# --- Output writers -----------------------------------------------------


def write_csv(rows: list[dict], fields: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_tsv(rows: list[dict], fields: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# --- Orchestration -----------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "cluster_and_split.log"

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


def require_mmseqs() -> None:
    if not MMSEQS_BIN.exists():
        raise FileNotFoundError(
            f"MMseqs2 binary not found at {MMSEQS_BIN}. Per PLAN.md/CLAUDE.md, leakage "
            "control must be sequence-level -- do not fall back to exact-ID matching. "
            "Install the static binary (see PLAN.md Sec 6 Risks) before rerunning."
        )


def run() -> None:
    _setup_logging()
    require_mmseqs()
    attrition: list[dict] = []

    candidate_rows = load_candidates()
    candidate_index = {chain_id_str(r["pdb_id"], r["chain_id"]): r for r in candidate_rows}
    n_initial = len(candidate_rows)
    attrition.append(
        {"stage": "initial_candidate_chains", "n_before": n_initial, "n_dropped": 0, "n_remaining": n_initial, "note": ""}
    )

    sequences, dropped_missing_seq = load_candidate_sequences(candidate_rows)
    attrition.append(
        {
            "stage": "missing_cached_sequence",
            "n_before": n_initial,
            "n_dropped": len(dropped_missing_seq),
            "n_remaining": len(sequences),
            "note": "chain had no cached polymer-entity JSON / canonical sequence",
        }
    )
    write_fasta(sequences, SEQUENCES_FASTA_PATH)
    logger.info("Wrote %d candidate sequences to %s", len(sequences), SEQUENCES_FASTA_PATH)

    # --- Redundancy clustering ---
    cluster_tsv = run_mmseqs_cluster(SEQUENCES_FASTA_PATH, CLUSTER_OUT_PREFIX, CLUSTER_TMP_DIR)
    groups = parse_cluster_tsv(cluster_tsv)
    logger.info("MMseqs2 clustering produced %d clusters from %d sequences", len(groups), len(sequences))

    cluster_rows = []
    representatives: list[str] = []
    for members in groups.values():
        members = [m for m in members if m in candidate_index]
        if not members:
            continue
        rep_id = choose_representative(members, candidate_index)
        representatives.append(rep_id)
        for member_id in members:
            row = candidate_index[member_id]
            cluster_rows.append(
                {
                    "cluster_id": rep_id,
                    "pdb_id": row["pdb_id"],
                    "chain_id": row["chain_id"],
                    "is_representative": member_id == rep_id,
                }
            )
    write_tsv(cluster_rows, CLUSTERS_FIELDS, CLUSTERS_PATH)
    logger.info("Wrote %d cluster-membership rows to %s", len(cluster_rows), CLUSTERS_PATH)

    attrition.append(
        {
            "stage": "redundancy_clustering",
            "n_before": len(sequences),
            "n_dropped": len(sequences) - len(representatives),
            "n_remaining": len(representatives),
            "note": f"MMseqs2 easy-cluster --min-seq-id {config.SEQUENCE_IDENTITY_CUTOFF} -c {config.CLUSTER_MIN_COVERAGE} --cov-mode 0; "
            "one deterministic representative kept per cluster (best resolution, then longest chain, then ID sort)",
        }
    )

    # --- PeSTo split sequences + homology search ---
    split_membership = load_pesto_split_membership()
    split_union = set().union(*split_membership.values())
    logger.info("PeSTo train+test+validation union: %d distinct chains", len(split_union))

    session = build_session()
    seqres_path = download_pdb_seqres(session)
    pesto_sequences, missing_pesto = lookup_pdb_seqres_sequences(seqres_path, split_union)
    if missing_pesto:
        logger.warning(
            "%d/%d PeSTo split chains had no sequence in pdb_seqres.txt.gz (obsolete/superseded "
            "entries) -- excluded from the homology search target set, not silently merged in",
            len(missing_pesto),
            len(split_union),
        )
    write_fasta(pesto_sequences, PESTO_ALL_SPLITS_FASTA)
    logger.info("Wrote %d PeSTo split sequences to %s", len(pesto_sequences), PESTO_ALL_SPLITS_FASTA)

    rep_sequences = {rep_id: sequences[rep_id] for rep_id in representatives}
    write_fasta(rep_sequences, REPRESENTATIVES_FASTA)

    search_out = run_mmseqs_search(REPRESENTATIVES_FASTA, PESTO_ALL_SPLITS_FASTA, SEARCH_OUT_PATH, SEARCH_TMP_DIR)
    search_hits = parse_search_tsv(search_out)
    n_reps_with_hits = sum(1 for r in representatives if search_hits.get(r))
    logger.info(
        "MMseqs2 homology search: %d/%d representatives had >=1 hit against the PeSTo split union at identity>=%.2f",
        n_reps_with_hits,
        len(representatives),
        min(config.LEAKAGE_SWEEP_IDENTITY_THRESHOLDS),
    )

    flags = compute_leakage_flags(representatives, search_hits, split_membership)
    n_homolog_overlap = sum(1 for f in flags.values() if f["pesto_homolog_overlap"])
    n_exact_overlap = sum(1 for f in flags.values() if f["pesto_exact_train_overlap"])
    logger.info(
        "pesto_homolog_overlap=True for %d/%d representatives (>=%.2f identity to train+test+validation union); "
        "pesto_exact_train_overlap=True for %d/%d",
        n_homolog_overlap,
        len(representatives),
        config.SEQUENCE_IDENTITY_CUTOFF,
        n_exact_overlap,
        len(representatives),
    )

    homology_rows = [{"rep_id": rep_id, **flags[rep_id]} for rep_id in representatives]
    write_tsv(homology_rows, HOMOLOGY_FIELDS, PESTO_HOMOLOGY_SEARCH_PATH)
    logger.info("Wrote per-representative homology search results to %s", PESTO_HOMOLOGY_SEARCH_PATH)

    # --- Combined dedup + leakage output ---
    cluster_size = {}
    for row in cluster_rows:
        cluster_size[row["cluster_id"]] = cluster_size.get(row["cluster_id"], 0) + 1

    dedup_rows = []
    for rep_id in representatives:
        row = dict(candidate_index[rep_id])
        row["cluster_id"] = rep_id
        row["cluster_size"] = cluster_size[rep_id]
        row.update(flags[rep_id])
        dedup_rows.append(row)
    write_csv(dedup_rows, DEDUP_FIELDS, CANDIDATES_DEDUP_PATH)
    logger.info("Wrote %d deduplicated+flagged representative chains to %s", len(dedup_rows), CANDIDATES_DEDUP_PATH)

    # --- Threshold sweep ---
    sweep_rows = sweep_thresholds(flags)
    write_csv(sweep_rows, SWEEP_FIELDS, LEAKAGE_SWEEP_PATH)
    logger.info("Wrote leakage threshold sweep to %s", LEAKAGE_SWEEP_PATH)
    for row in sweep_rows:
        logger.info(
            "sweep %-20s n_survive=%-6d meets_target(%d)=%-5s meets_floor(%d)=%s",
            row["mode"],
            row["n_survive"],
            config.MIN_TEST_CHAINS_TARGET,
            row["meets_target"],
            config.MIN_TEST_CHAINS_FLOOR,
            row["meets_floor"],
        )

    write_csv(attrition, ATTRITION_FIELDS, ATTRITION_PATH)
    logger.info("Wrote attrition table to %s", ATTRITION_PATH)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Phase 2: redundancy reduction + PeSTo-overlap flagging")
    parser.parse_args(argv)
    run()


if __name__ == "__main__":
    main()
