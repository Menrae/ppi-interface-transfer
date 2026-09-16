"""Phase 3: download experimental (PDBe updated mmCIF) and AlphaFold DB
structures for every Phase 2 candidate representative.

Downloads all 3,009 representatives in ``data/interim/candidates_dedup.csv``
(so the leakage sensitivity analyses have inputs too), processing the
primary leakage-filtered set (``eligible_homolog``) first so the run is
usable if interrupted. For each representative: the PDBe "updated" mmCIF
for its PDB entry (embeds SIFTS residue-level UniProt cross-references,
reused directly by Phase 4 -- see PLAN.md Phase 3/4) and the AlphaFold DB
model for its UniProt accession, queried via the AlphaFold DB prediction
API (never a hardcoded model version/URL pattern).

Before the bulk download, a small real (not synthetic) calibration sample
is fetched to empirically project total wall-clock time and disk usage; if
that projection exceeds ``config.PHASE3_TIME_BUDGET_SECONDS`` or
``config.PHASE3_DISK_BUDGET_BYTES`` the run stops before downloading the
rest, per PLAN.md Phase 3.

Runnable as: .venv/bin/python -m src.data.fetch_structures
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import time
from pathlib import Path

import requests

from src import config
from src.data.select_complexes import build_session

logger = logging.getLogger(__name__)

# --- URL templates -----------------------------------------------------

PDBE_UPDATED_CIF_URL_TEMPLATE = (
    "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id_lower}_updated.cif"
)
ALPHAFOLD_API_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
UNIPROT_FASTA_URL_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"

# --- Cache/output paths --------------------------------------------------

PDB_UPDATED_DIR = config.RAW_DATA_DIR / "pdb"
ALPHAFOLD_DIR = config.RAW_DATA_DIR / "alphafold"
ALPHAFOLD_API_CACHE_DIR = ALPHAFOLD_DIR / "api"
UNIPROT_CACHE_DIR = config.RAW_DATA_DIR / "uniprot"

CANDIDATES_DEDUP_PATH = config.INTERIM_DATA_DIR / "candidates_dedup.csv"
FETCH_REPORT_PATH = config.INTERIM_DATA_DIR / "fetch_report.csv"
ATTRITION_PATH = config.INTERIM_DATA_DIR / "phase3_attrition.csv"

# --- Failure reasons (logged distinctly, never silently dropped) --------

FAILURE_ACCESSION_ABSENT = "accession_absent_from_alphafold_db"
FAILURE_FRAGMENTED = "fragmented_in_alphafold_db"
FAILURE_SEQUENCE_MISMATCH = "sequence_mismatch"
FAILURE_NETWORK_ERROR = "network_error"
FAILURE_PDB_DOWNLOAD_ERROR = "pdb_download_error"

FETCH_REPORT_FIELDS = [
    "cluster_id",
    "pdb_id",
    "chain_id",
    "uniprot_acc",
    "release_date",
    "eligible_homolog",
    "eligible_exact_train",
    "eligible_none",
    "pdb_cif_path",
    "pdb_status",
    "pdb_failure_reason",
    "alphafold_cif_path",
    "alphafold_model_version",
    "alphafold_uniprot_start",
    "alphafold_uniprot_end",
    "alphafold_status",
    "alphafold_failure_reason",
    "overall_status",
]

ATTRITION_FIELDS = [
    "leakage_filter_mode",
    "n_eligible_representatives",
    "n_pdb_success",
    "n_pdb_failed",
    "n_alphafold_success",
    "n_alphafold_failed",
    "n_both_success",
    "note",
]


# --- Pure functions (URL construction, response classification) --------
# Kept free of I/O so tests/test_fetch_structures.py can exercise them
# directly on saved fixtures with no network access.


def pdbe_updated_cif_url(pdb_id: str) -> str:
    return PDBE_UPDATED_CIF_URL_TEMPLATE.format(pdb_id_lower=pdb_id.lower())


def alphafold_api_url(accession: str) -> str:
    return ALPHAFOLD_API_URL_TEMPLATE.format(accession=accession)


def uniprot_fasta_url(accession: str) -> str:
    return UNIPROT_FASTA_URL_TEMPLATE.format(accession=accession)


def classify_pdb_http_status(status_code: int) -> str | None:
    """Return a failure reason string, or None if status_code means success."""
    if status_code == 200:
        return None
    return f"{FAILURE_PDB_DOWNLOAD_ERROR}:{status_code}"


def classify_alphafold_api_response(
    status_code: int, body: list[dict] | None, accession: str
) -> tuple[str | None, dict | None]:
    """Classify one AlphaFold DB prediction-API response.

    Returns (failure_reason, matched_entry): exactly one is truthy/None+dict
    on success.

    The API can return entries for *isoforms* of the queried accession (e.g.
    querying "Q9UPN9" can return both "AF-Q9UPN9-F1" (uniprotAccession
    "Q9UPN9") and "AF-Q9UPN9-2-F1" (uniprotAccession "Q9UPN9-2")) -- only an
    entry whose own ``uniprotAccession`` is an exact match to the queried
    accession counts. If AFDB has genuinely split the queried accession's
    own sequence into multiple fragment models (>1 exact match), that is
    ``FAILURE_FRAGMENTED``, not success -- fragments are not stitched here.
    """
    if status_code != 200 or not body:
        return FAILURE_ACCESSION_ABSENT, None
    matches = [entry for entry in body if entry.get("uniprotAccession") == accession]
    if not matches:
        return FAILURE_ACCESSION_ABSENT, None
    if len(matches) > 1:
        return FAILURE_FRAGMENTED, None
    return None, matches[0]


def parse_fasta_sequence(text: str) -> str:
    """Return the concatenated sequence (no header) from single-record FASTA text."""
    lines = text.strip().splitlines()
    return "".join(line.strip() for line in lines if line and not line.startswith(">"))


def sequences_match(model_sequence: str, current_sequence: str) -> bool:
    return bool(model_sequence) and model_sequence == current_sequence


def leakage_priority_rank(row: dict) -> tuple[int, str]:
    """Sort key so the primary (eligible_homolog) set downloads first.

    Rank 0 = primary set (eligible_homolog); rank 1 = eligible_exact_train
    but not primary; rank 2 = everything else (eligible_none only). Ties
    broken by cluster_id for a fully deterministic, resumable order.
    """
    if row["eligible_homolog"]:
        rank = 0
    elif row["eligible_exact_train"]:
        rank = 1
    else:
        rank = 2
    return (rank, row["cluster_id"])


def build_fetch_report_row(rep: dict, pdb_result: dict, af_result: dict) -> dict:
    if pdb_result["status"] == "success" and af_result["status"] == "success":
        overall = "complete"
    elif pdb_result["status"] == "success" or af_result["status"] == "success":
        overall = "partial"
    else:
        overall = "failed"
    return {
        "cluster_id": rep["cluster_id"],
        "pdb_id": rep["pdb_id"],
        "chain_id": rep["chain_id"],
        "uniprot_acc": rep["uniprot_acc"],
        "release_date": rep["release_date"],
        "eligible_homolog": rep["eligible_homolog"],
        "eligible_exact_train": rep["eligible_exact_train"],
        "eligible_none": rep["eligible_none"],
        "pdb_cif_path": pdb_result["path"],
        "pdb_status": pdb_result["status"],
        "pdb_failure_reason": pdb_result["failure_reason"],
        "alphafold_cif_path": af_result["cif_path"],
        "alphafold_model_version": af_result["model_version"],
        "alphafold_uniprot_start": af_result["uniprot_start"],
        "alphafold_uniprot_end": af_result["uniprot_end"],
        "alphafold_status": af_result["status"],
        "alphafold_failure_reason": af_result["failure_reason"],
        "overall_status": overall,
    }


# --- I/O: cached fetch functions (single source of truth for both the ---
# --- calibration sample and the full run) -------------------------------


def fetch_pdb_structure(session: requests.Session, pdb_id: str) -> dict:
    """Ensure data/raw/pdb/{PDB_ID}_updated.cif.gz exists; return a status dict.

    Cache hit -> zero network calls. Only ever caches on HTTP 200; a
    transient network exception or non-200 response leaves no cache file,
    so a rerun retries it (resumable, consistent with Phase 1/2).
    """
    cache_path = PDB_UPDATED_DIR / f"{pdb_id.upper()}_updated.cif.gz"
    if cache_path.exists():
        return {
            "path": str(cache_path),
            "status": "success",
            "failure_reason": "",
            "bytes": cache_path.stat().st_size,
            "network_call": False,
        }

    time.sleep(config.PHASE3_REQUEST_DELAY_SECONDS)
    url = pdbe_updated_cif_url(pdb_id)
    try:
        response = session.get(url, timeout=60)
    except requests.RequestException as exc:
        logger.warning("network error fetching PDB %s: %s", pdb_id, exc)
        return {
            "path": "",
            "status": "failed",
            "failure_reason": FAILURE_NETWORK_ERROR,
            "bytes": 0,
            "network_call": True,
        }

    failure = classify_pdb_http_status(response.status_code)
    if failure is not None:
        logger.warning("PDB download error for %s: HTTP %d", pdb_id, response.status_code)
        return {
            "path": "",
            "status": "failed",
            "failure_reason": failure,
            "bytes": 0,
            "network_call": True,
        }

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
    with gzip.open(tmp_path, "wb") as f:
        f.write(response.content)
    tmp_path.rename(cache_path)
    return {
        "path": str(cache_path),
        "status": "success",
        "failure_reason": "",
        "bytes": len(response.content),
        "network_call": True,
    }


def _cached_alphafold_api_response(session: requests.Session, accession: str) -> tuple[int, list | None, bool]:
    """Return (status_code, body, made_network_call).

    Caches both success (200) and permanent-failure (non-200) outcomes --
    "accession absent from AlphaFold DB" is a common, stable outcome here,
    not a transient error, so it's worth memoizing to keep reruns from
    re-querying every permanently-absent accession. A network exception is
    never cached (propagates to the caller, retried on the next run).
    """
    cache_path = ALPHAFOLD_API_CACHE_DIR / f"{accession}.json"
    if cache_path.exists():
        with open(cache_path) as f:
            cached = json.load(f)
        return cached["status_code"], cached["body"], False

    time.sleep(config.PHASE3_REQUEST_DELAY_SECONDS)
    response = session.get(alphafold_api_url(accession), timeout=30)
    status_code = response.status_code
    body = response.json() if status_code == 200 else None
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump({"status_code": status_code, "body": body}, f)
    return status_code, body, True


def _cached_current_uniprot_sequence(session: requests.Session, accession: str) -> tuple[str, bool]:
    """Return (sequence, made_network_call); caches the raw FASTA text on success."""
    cache_path = UNIPROT_CACHE_DIR / f"{accession}.fasta"
    if cache_path.exists():
        with open(cache_path) as f:
            return parse_fasta_sequence(f.read()), False

    time.sleep(config.PHASE3_REQUEST_DELAY_SECONDS)
    response = session.get(uniprot_fasta_url(accession), timeout=30)
    if response.status_code != 200:
        raise requests.RequestException(
            f"UniProt fetch for {accession} returned HTTP {response.status_code}"
        )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        f.write(response.text)
    return parse_fasta_sequence(response.text), True


def fetch_alphafold_model(session: requests.Session, accession: str) -> dict:
    """Ensure the AlphaFold model (+ metadata) for accession is cached.

    Three steps, each of which can fail distinctly (see module-level
    FAILURE_* constants): (1) query the prediction API and match the exact
    accession (not an isoform) with exactly one fragment; (2) compare the
    matched model's UniProt sequence against the accession's *current*
    UniProt sequence (a mismatch means the model may be numbered against an
    outdated sequence -- unsafe for Phase 4's residue mapping, so the mmCIF
    is not downloaded); (3) download the mmCIF itself.
    """
    result = {
        "cif_path": "",
        "status": "failed",
        "failure_reason": "",
        "model_version": "",
        "uniprot_start": "",
        "uniprot_end": "",
        "bytes": 0,
        "network_call": False,
    }

    try:
        api_status, api_body, made_call = _cached_alphafold_api_response(session, accession)
    except requests.RequestException as exc:
        logger.warning("network error fetching AlphaFold API for %s: %s", accession, exc)
        result["failure_reason"] = FAILURE_NETWORK_ERROR
        result["network_call"] = True
        return result
    result["network_call"] = result["network_call"] or made_call

    failure, entry = classify_alphafold_api_response(api_status, api_body, accession)
    if failure is not None:
        logger.warning("AlphaFold DB lookup failed for %s: %s", accession, failure)
        result["failure_reason"] = failure
        return result

    result["model_version"] = entry.get("latestVersion", "")
    result["uniprot_start"] = entry.get("uniprotStart", "")
    result["uniprot_end"] = entry.get("uniprotEnd", "")

    try:
        current_seq, made_call = _cached_current_uniprot_sequence(session, accession)
    except requests.RequestException as exc:
        logger.warning("network error fetching current UniProt sequence for %s: %s", accession, exc)
        result["failure_reason"] = FAILURE_NETWORK_ERROR
        result["network_call"] = True
        return result
    result["network_call"] = result["network_call"] or made_call

    model_seq = entry.get("uniprotSequence", "")
    if not sequences_match(model_seq, current_seq):
        logger.warning(
            "sequence mismatch for %s: AlphaFold model UniProt sequence (%d aa) != "
            "current UniProt sequence (%d aa)",
            accession,
            len(model_seq),
            len(current_seq),
        )
        result["failure_reason"] = FAILURE_SEQUENCE_MISMATCH
        return result

    cif_cache_path = ALPHAFOLD_DIR / f"{accession}.cif.gz"
    if cif_cache_path.exists():
        result["cif_path"] = str(cif_cache_path)
        result["status"] = "success"
        result["bytes"] = cif_cache_path.stat().st_size
        return result

    result["network_call"] = True
    time.sleep(config.PHASE3_REQUEST_DELAY_SECONDS)
    try:
        response = session.get(entry["cifUrl"], timeout=60)
    except requests.RequestException as exc:
        logger.warning("network error downloading AlphaFold model for %s: %s", accession, exc)
        result["failure_reason"] = FAILURE_NETWORK_ERROR
        return result
    if response.status_code != 200:
        logger.warning(
            "AlphaFold model download error for %s: HTTP %d", accession, response.status_code
        )
        result["failure_reason"] = FAILURE_NETWORK_ERROR
        return result

    cif_cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cif_cache_path.with_suffix(cif_cache_path.suffix + ".tmp")
    with gzip.open(tmp_path, "wb") as f:
        f.write(response.content)
    tmp_path.rename(cif_cache_path)
    result["cif_path"] = str(cif_cache_path)
    result["status"] = "success"
    result["bytes"] = len(response.content)
    return result


# --- Loading candidates_dedup.csv ---------------------------------------


def _csv_bool(value: str) -> bool:
    return value == "True"


def load_representatives(path: Path = CANDIDATES_DEDUP_PATH) -> list[dict]:
    """Load candidates_dedup.csv rows with the fields fetch_structures needs.

    Booleans are parsed immediately (not left as "True"/"False" strings) so
    downstream code never re-derives them.
    """
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append(
                {
                    "cluster_id": row["cluster_id"],
                    "pdb_id": row["pdb_id"],
                    "chain_id": row["chain_id"],
                    "uniprot_acc": row["uniprot_ids"],
                    "release_date": row["release_date"],
                    "eligible_homolog": _csv_bool(row["eligible_homolog"]),
                    "eligible_exact_train": _csv_bool(row["eligible_exact_train"]),
                    "eligible_none": _csv_bool(row["eligible_none"]),
                }
            )
    return rows


def sort_by_priority(rows: list[dict]) -> list[dict]:
    return sorted(rows, key=leakage_priority_rank)


def unique_ids_in_order(rows: list[dict], key: str) -> list[str]:
    """Unique values of rows[i][key], first-seen order preserved."""
    seen: set[str] = set()
    ordered: list[str] = []
    for row in rows:
        value = row[key]
        if value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


# --- Pre-flight calibration/estimate -------------------------------------


def calibrate_and_estimate(
    session: requests.Session,
    pdb_ids_priority_order: list[str],
    uniprot_accs_priority_order: list[str],
    sample_size: int = None,
) -> dict:
    """Fetch a small real, uncached sample to project full-run time/disk.

    Only *uncached* ids are sampled and measured -- already-downloaded ids
    cost ~0 additional time, so including them would bias the estimate low
    on a rerun with most files already in place. The sampled fetches are
    real (not thrown away): they populate the same cache the full run uses.
    """
    if sample_size is None:
        sample_size = config.PHASE3_ESTIMATE_SAMPLE_SIZE

    uncached_pdb = [
        pid for pid in pdb_ids_priority_order
        if not (PDB_UPDATED_DIR / f"{pid.upper()}_updated.cif.gz").exists()
    ]
    uncached_uniprot = [
        acc for acc in uniprot_accs_priority_order
        if not (ALPHAFOLD_DIR / f"{acc}.cif.gz").exists()
    ]
    n_pdb_remaining = len(uncached_pdb)
    n_uniprot_remaining = len(uncached_uniprot)

    pdb_sample = uncached_pdb[:sample_size]
    uniprot_sample = uncached_uniprot[:sample_size]

    t0 = time.time()
    pdb_bytes = 0
    for pdb_id in pdb_sample:
        pdb_bytes += fetch_pdb_structure(session, pdb_id)["bytes"]
    pdb_elapsed = time.time() - t0
    pdb_seconds_per_item = pdb_elapsed / len(pdb_sample) if pdb_sample else 0.0
    pdb_bytes_per_item = pdb_bytes / len(pdb_sample) if pdb_sample else 0.0

    t0 = time.time()
    af_bytes = 0
    for acc in uniprot_sample:
        af_bytes += fetch_alphafold_model(session, acc)["bytes"]
    af_elapsed = time.time() - t0
    af_seconds_per_item = af_elapsed / len(uniprot_sample) if uniprot_sample else 0.0
    af_bytes_per_item = af_bytes / len(uniprot_sample) if uniprot_sample else 0.0

    estimated_seconds = (
        pdb_seconds_per_item * n_pdb_remaining + af_seconds_per_item * n_uniprot_remaining
    )
    estimated_bytes = (
        pdb_bytes_per_item * n_pdb_remaining + af_bytes_per_item * n_uniprot_remaining
    )

    return {
        "n_unique_pdb_total": len(pdb_ids_priority_order),
        "n_unique_uniprot_total": len(uniprot_accs_priority_order),
        "n_pdb_remaining": n_pdb_remaining,
        "n_uniprot_remaining": n_uniprot_remaining,
        "n_pdb_calibrated": len(pdb_sample),
        "n_uniprot_calibrated": len(uniprot_sample),
        "pdb_seconds_per_item": pdb_seconds_per_item,
        "af_seconds_per_item": af_seconds_per_item,
        "pdb_bytes_per_item": pdb_bytes_per_item,
        "af_bytes_per_item": af_bytes_per_item,
        "estimated_seconds": estimated_seconds,
        "estimated_bytes": estimated_bytes,
        "over_time_budget": estimated_seconds > config.PHASE3_TIME_BUDGET_SECONDS,
        "over_disk_budget": estimated_bytes > config.PHASE3_DISK_BUDGET_BYTES,
    }


def log_estimate(estimate: dict) -> None:
    logger.info(
        "Phase 3 scope: %d unique PDB entries (%d not yet cached), %d unique UniProt "
        "accessions (%d not yet cached)",
        estimate["n_unique_pdb_total"],
        estimate["n_pdb_remaining"],
        estimate["n_unique_uniprot_total"],
        estimate["n_uniprot_remaining"],
    )
    logger.info(
        "Calibration sample: %d PDB fetches (%.2fs/item), %d AlphaFold fetches (%.2fs/item)",
        estimate["n_pdb_calibrated"],
        estimate["pdb_seconds_per_item"],
        estimate["n_uniprot_calibrated"],
        estimate["af_seconds_per_item"],
    )
    logger.info(
        "Projected remaining work: %.1f minutes, %.2f GB (raw, pre-gzip) -- budget is "
        "%.1f minutes / %.1f GB",
        estimate["estimated_seconds"] / 60,
        estimate["estimated_bytes"] / 1024**3,
        config.PHASE3_TIME_BUDGET_SECONDS / 60,
        config.PHASE3_DISK_BUDGET_BYTES / 1024**3,
    )
    if estimate["over_time_budget"] or estimate["over_disk_budget"]:
        logger.error(
            "Projected time/disk usage exceeds the configured budget -- stopping before "
            "the bulk download. Re-run with a smaller scope or raise "
            "PHASE3_TIME_BUDGET_SECONDS/PHASE3_DISK_BUDGET_BYTES in src/config.py after "
            "confirming with the user."
        )


# --- Orchestration -------------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "fetch_structures.log"

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
        n_pdb_success = sum(1 for r in subset if r["pdb_status"] == "success")
        n_af_success = sum(1 for r in subset if r["alphafold_status"] == "success")
        n_both = sum(
            1 for r in subset if r["pdb_status"] == "success" and r["alphafold_status"] == "success"
        )
        attrition.append(
            {
                "leakage_filter_mode": mode,
                "n_eligible_representatives": len(subset),
                "n_pdb_success": n_pdb_success,
                "n_pdb_failed": len(subset) - n_pdb_success,
                "n_alphafold_success": n_af_success,
                "n_alphafold_failed": len(subset) - n_af_success,
                "n_both_success": n_both,
                "note": f"meets_target({config.MIN_TEST_CHAINS_TARGET})="
                f"{n_both >= config.MIN_TEST_CHAINS_TARGET}",
            }
        )
    return attrition


def run(estimate_only: bool = False, force: bool = False) -> None:
    _setup_logging()
    session = build_session()

    rows = sort_by_priority(load_representatives())
    logger.info("Loaded %d Phase 2 representatives from %s", len(rows), CANDIDATES_DEDUP_PATH)

    pdb_ids = unique_ids_in_order(rows, "pdb_id")
    uniprot_accs = unique_ids_in_order(rows, "uniprot_acc")

    estimate = calibrate_and_estimate(session, pdb_ids, uniprot_accs)
    log_estimate(estimate)

    if estimate_only:
        return
    if (estimate["over_time_budget"] or estimate["over_disk_budget"]) and not force:
        raise SystemExit(
            "Phase 3 pre-flight estimate exceeded the configured time/disk budget -- "
            "stopping before the bulk download. See the logged estimate above."
        )

    report_rows = []
    for i, rep in enumerate(rows):
        pdb_result = fetch_pdb_structure(session, rep["pdb_id"])
        af_result = fetch_alphafold_model(session, rep["uniprot_acc"])
        report_rows.append(build_fetch_report_row(rep, pdb_result, af_result))
        if (i + 1) % 250 == 0 or (i + 1) == len(rows):
            logger.info("Processed %d/%d representatives", i + 1, len(rows))

    write_csv(report_rows, FETCH_REPORT_FIELDS, FETCH_REPORT_PATH)
    logger.info("Wrote fetch report (%d rows) to %s", len(report_rows), FETCH_REPORT_PATH)

    attrition = compute_attrition(report_rows)
    write_csv(attrition, ATTRITION_FIELDS, ATTRITION_PATH)
    logger.info("Wrote per-leakage-mode attrition to %s", ATTRITION_PATH)
    for row in attrition:
        logger.info(
            "mode=%-14s n_eligible=%-6d n_both_success=%-6d %s",
            row["leakage_filter_mode"],
            row["n_eligible_representatives"],
            row["n_both_success"],
            row["note"],
        )

    failure_counts: dict[str, int] = {}
    for r in report_rows:
        for reason in (r["pdb_failure_reason"], r["alphafold_failure_reason"]):
            if reason:
                key = reason.split(":", 1)[0]
                failure_counts[key] = failure_counts.get(key, 0) + 1
    logger.info("Failure-reason breakdown: %s", failure_counts)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Phase 3: download PDB (updated mmCIF) and AlphaFold DB structures"
    )
    parser.add_argument(
        "--estimate-only",
        action="store_true",
        help="Run the pre-flight calibration/estimate and stop (no bulk download)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Proceed with the bulk download even if the pre-flight estimate exceeds budget",
    )
    args = parser.parse_args(argv)
    run(estimate_only=args.estimate_only, force=args.force)


if __name__ == "__main__":
    main()
