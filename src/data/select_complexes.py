"""Phase 1: candidate complex selection.

Queries the RCSB Search API for X-ray entries meeting the resolution,
release-date, protein-entity-count, and nucleic-acid-free filters from
``src.config``, fetches per-entry and per-polymer-entity detail from the
RCSB Data API, cross-checks chain-to-UniProt mapping against the bulk SIFTS
file, and writes one row per surviving protein chain to
``data/interim/candidates.csv``.

Runnable as: .venv/bin/python -m src.data.select_complexes --max-entries 500
"""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import logging
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

from src import config

logger = logging.getLogger(__name__)

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{entry_id}"
POLYMER_ENTITY_URL_TEMPLATE = (
    "https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}"
)
SIFTS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/"
    "pdb_chain_uniprot.tsv.gz"
)

RCSB_RAW_DIR = config.RAW_DATA_DIR / "rcsb"
SEARCH_CACHE_DIR = RCSB_RAW_DIR / "search"
ENTRY_CACHE_DIR = RCSB_RAW_DIR / "entries"
POLYMER_ENTITY_CACHE_DIR = RCSB_RAW_DIR / "polymer_entities"

SIFTS_DIR = config.RAW_DATA_DIR / "sifts"
SIFTS_CACHE_PATH = SIFTS_DIR / "pdb_chain_uniprot.tsv.gz"

CANDIDATES_PATH = config.INTERIM_DATA_DIR / "candidates.csv"
ATTRITION_PATH = config.INTERIM_DATA_DIR / "phase1_attrition.csv"

SEARCH_PAGE_SIZE = 100
REQUEST_DELAY_SECONDS = 0.1

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

ATTRITION_FIELDS = ["stage", "n_before", "n_dropped", "n_remaining", "note"]


# --- HTTP session with retry/backoff -------------------------------------


def build_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=5,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def _cached_request(
    session: requests.Session,
    cache_path: Path,
    url: str,
    method: str = "GET",
    json_body: dict | None = None,
) -> dict:
    """GET/POST url, caching the parsed JSON response at cache_path.

    A cache hit means zero network calls, so reruns work offline.
    """
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    time.sleep(REQUEST_DELAY_SECONDS)
    if method == "POST":
        response = session.post(
            url, json=json_body, headers={"Content-Type": "application/json"}
        )
    else:
        response = session.get(url)
    response.raise_for_status()
    data = response.json()

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(data, f)
    return data


def _query_hash(query: dict) -> str:
    payload = json.dumps(query, sort_keys=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]


# --- RCSB Search API -------------------------------------------------------


def build_search_query(max_protein_entities: int | None = None) -> dict:
    """Build the RCSB Search API v2 query body from config thresholds.

    Filters: X-ray method, resolution <= RESOLUTION_CUTOFF_ANGSTROM, release
    date > PDB_RELEASE_DATE_CUTOFF, 2..max_protein_entities distinct protein
    entities, zero nucleic-acid and zero nucleic-acid-hybrid entities.
    """
    if max_protein_entities is None:
        max_protein_entities = config.MAX_PROTEIN_ENTITIES

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": config.RESOLUTION_CUTOFF_ANGSTROM,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "greater",
                        "value": config.PDB_RELEASE_DATE_CUTOFF,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "range",
                        "value": {
                            "from": 2,
                            "to": max_protein_entities,
                            "include_lower": True,
                            "include_upper": True,
                        },
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_nucleic_acid",
                        "operator": "equals",
                        "value": 0,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_nucleic_acid_hybrid",
                        "operator": "equals",
                        "value": 0,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_id", "direction": "asc"}],
        },
    }


def search_entries(
    session: requests.Session, max_entries: int | None = None
) -> tuple[list[str], int]:
    """Run the paginated search query, returning (entry_ids, total_count).

    Page size is fixed regardless of max_entries so cached pages are reused
    across runs with different --max-entries values.
    """
    base_query = build_search_query()
    qhash = _query_hash(base_query)
    cache_dir = SEARCH_CACHE_DIR / qhash

    entry_ids: list[str] = []
    total_count = 0
    start = 0

    while True:
        query = copy.deepcopy(base_query)
        query["request_options"]["paginate"] = {
            "start": start,
            "rows": SEARCH_PAGE_SIZE,
        }
        cache_path = cache_dir / f"page_start{start:06d}.json"
        data = _cached_request(
            session, cache_path, SEARCH_URL, method="POST", json_body=query
        )

        total_count = data.get("total_count", 0)
        hits = data.get("result_set", [])
        entry_ids.extend(hit["identifier"] for hit in hits)
        logger.info(
            "search page start=%d -> %d hits (total_count=%d, collected=%d)",
            start,
            len(hits),
            total_count,
            len(entry_ids),
        )

        if not hits or start + SEARCH_PAGE_SIZE >= total_count:
            break
        if max_entries is not None and len(entry_ids) >= max_entries:
            break
        start += SEARCH_PAGE_SIZE

    if max_entries is not None:
        entry_ids = entry_ids[:max_entries]
    return entry_ids, total_count


# --- RCSB Data API ---------------------------------------------------------


def fetch_entry(session: requests.Session, entry_id: str) -> dict:
    cache_path = ENTRY_CACHE_DIR / f"{entry_id}.json"
    return _cached_request(
        session, cache_path, ENTRY_URL_TEMPLATE.format(entry_id=entry_id)
    )


def fetch_polymer_entity(
    session: requests.Session, entry_id: str, entity_id: str
) -> dict:
    cache_path = POLYMER_ENTITY_CACHE_DIR / f"{entry_id}_{entity_id}.json"
    return _cached_request(
        session,
        cache_path,
        POLYMER_ENTITY_URL_TEMPLATE.format(entry_id=entry_id, entity_id=entity_id),
    )


def parse_entry_summary(entry_id: str, entry_data: dict) -> dict:
    accession = entry_data.get("rcsb_accession_info", {}) or {}
    info = entry_data.get("rcsb_entry_info", {}) or {}
    identifiers = entry_data.get("rcsb_entry_container_identifiers", {}) or {}
    resolution_list = info.get("resolution_combined") or []

    release_date = accession.get("initial_release_date") or ""
    return {
        "pdb_id": entry_id,
        "release_date": release_date[:10] if release_date else "",
        "resolution": resolution_list[0] if resolution_list else None,
        "polymer_entity_count_protein": info.get("polymer_entity_count_protein"),
        "polymer_entity_count_nucleic_acid": info.get(
            "polymer_entity_count_nucleic_acid"
        ),
        "polymer_entity_ids": identifiers.get("polymer_entity_ids", []) or [],
    }


def parse_polymer_entity(entry_id: str, entity_id: str, entity_data: dict) -> list[dict]:
    """Return one dict per auth_asym_id (chain) for this polymer entity."""
    entity_poly = entity_data.get("entity_poly", {}) or {}
    container = entity_data.get("rcsb_polymer_entity_container_identifiers", {}) or {}

    polymer_type = entity_poly.get("rcsb_entity_polymer_type")
    seq_length = entity_poly.get("rcsb_sample_sequence_length")
    auth_asym_ids = container.get("auth_asym_ids", []) or []
    uniprot_ids = sorted(set(container.get("uniprot_ids") or []))

    rows = []
    for chain_id in auth_asym_ids:
        rows.append(
            {
                "pdb_id": entry_id,
                "entity_id": entity_id,
                "chain_id": chain_id,
                "polymer_type": polymer_type,
                "seq_length": seq_length,
                "uniprot_ids": list(uniprot_ids),
            }
        )
    return rows


# --- SIFTS bulk mapping ----------------------------------------------------


def download_sifts_mapping(session: requests.Session) -> Path:
    if SIFTS_CACHE_PATH.exists():
        return SIFTS_CACHE_PATH

    logger.info("Downloading SIFTS pdb_chain_uniprot.tsv.gz to %s", SIFTS_CACHE_PATH)
    SIFTS_DIR.mkdir(parents=True, exist_ok=True)
    time.sleep(REQUEST_DELAY_SECONDS)
    response = session.get(SIFTS_URL, stream=True)
    response.raise_for_status()

    tmp_path = SIFTS_CACHE_PATH.with_suffix(".tmp")
    with open(tmp_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=1 << 20):
            f.write(chunk)
    tmp_path.rename(SIFTS_CACHE_PATH)
    return SIFTS_CACHE_PATH


def load_sifts_mapping(path: Path) -> dict[tuple[str, str], tuple[str, ...]]:
    """Return {(PDB_ID upper, chain_id): (uniprot_acc, ...)} from the bulk file.

    The file has a '#'-prefixed version-comment line, a header line, then
    tab-separated PDB/CHAIN/SP_PRIMARY/... rows (possibly several rows per
    chain, one per aligned segment); accessions are deduped per chain.
    """
    mapping: dict[tuple[str, str], set[str]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        first_line = f.readline()
        if not first_line.startswith("#"):
            f.seek(0)
        f.readline()  # column header line (PDB CHAIN SP_PRIMARY ...)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            pdb_id, chain_id, sp_primary = parts[0], parts[1], parts[2]
            key = (pdb_id.upper(), chain_id)
            mapping.setdefault(key, set()).add(sp_primary)
    return {k: tuple(sorted(v)) for k, v in mapping.items()}


def cross_check_uniprot(
    pdb_id: str,
    chain_id: str,
    rcsb_uniprot_ids: list[str],
    sifts_mapping: dict[tuple[str, str], tuple[str, ...]],
) -> tuple[tuple[str, ...], bool]:
    """Return (sifts_uniprot_ids, agrees_with_rcsb)."""
    sifts_ids = sifts_mapping.get((pdb_id.upper(), chain_id), ())
    agrees = set(sifts_ids) == set(rcsb_uniprot_ids)
    return sifts_ids, agrees


# --- Chain-level filters (pure functions, no I/O) --------------------------


def filter_non_protein_entities(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    kept, dropped = [], []
    for row in rows:
        if row.get("polymer_type") == "Protein":
            kept.append(row)
        else:
            dropped.append(
                {**row, "drop_reason": f"non_protein_entity_type:{row.get('polymer_type')}"}
            )
    return kept, dropped


def filter_min_chain_length(rows: list[dict], min_length: int) -> tuple[list[dict], list[dict]]:
    kept, dropped = [], []
    for row in rows:
        length = row.get("seq_length")
        if length is not None and length >= min_length:
            kept.append(row)
        else:
            dropped.append(
                {**row, "drop_reason": f"chain_length<{min_length} (length={length})"}
            )
    return kept, dropped


def filter_uniprot_mapping(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    kept, dropped = [], []
    for row in rows:
        if row.get("uniprot_ids"):
            kept.append(row)
        else:
            dropped.append({**row, "drop_reason": "no_uniprot_mapping"})
    return kept, dropped


def filter_chimera(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    kept, dropped = [], []
    for row in rows:
        ids = set(row.get("uniprot_ids") or [])
        if len(ids) <= 1:
            kept.append(row)
        else:
            dropped.append({**row, "drop_reason": f"chimera:{','.join(sorted(ids))}"})
    return kept, dropped


def run_filter_pipeline(
    chain_rows: list[dict], min_chain_length: int
) -> tuple[list[dict], list[dict]]:
    """Apply all Phase 1 chain-level filters in order.

    Returns (surviving_rows, attrition_records). Every dropped row is
    logged at DEBUG with its reason (never silently dropped).
    """
    attrition: list[dict] = []
    n_before = len(chain_rows)
    attrition.append(
        {
            "stage": "initial_protein_chains",
            "n_before": n_before,
            "n_dropped": 0,
            "n_remaining": n_before,
            "note": "protein-type chains extracted from fetched entities",
        }
    )

    stages = [
        ("non_protein_entity_type", filter_non_protein_entities),
        (f"min_chain_length>={min_chain_length}", lambda rows: filter_min_chain_length(rows, min_chain_length)),
        ("uniprot_mapping", filter_uniprot_mapping),
        ("chimera", filter_chimera),
    ]

    rows = chain_rows
    for stage_name, filter_fn in stages:
        n_before = len(rows)
        rows, dropped = filter_fn(rows)
        for d in dropped:
            logger.debug(
                "dropped %s_%s (entity %s): %s",
                d["pdb_id"],
                d["chain_id"],
                d.get("entity_id"),
                d["drop_reason"],
            )
        attrition.append(
            {
                "stage": stage_name,
                "n_before": n_before,
                "n_dropped": len(dropped),
                "n_remaining": len(rows),
                "note": "",
            }
        )
        logger.info(
            "filter %-28s n_before=%-7d n_dropped=%-7d n_remaining=%-7d",
            stage_name,
            n_before,
            len(dropped),
            len(rows),
        )

    return rows, attrition


# --- Orchestration -----------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "select_complexes.log"

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)-8s %(name)s: %(message)s")
    )

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))

    root.handlers.clear()
    root.addHandler(file_handler)
    root.addHandler(console_handler)


def write_candidates(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CANDIDATE_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_attrition(records: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ATTRITION_FIELDS)
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def run(max_entries: int | None = None, max_protein_entities: int | None = None) -> None:
    _setup_logging()
    session = build_session()

    if max_protein_entities is not None:
        # build_search_query reads config directly; honor an override by
        # monkeypatching for the duration of this run rather than adding a
        # second source of truth for the threshold.
        original = config.MAX_PROTEIN_ENTITIES
        config.MAX_PROTEIN_ENTITIES = max_protein_entities
    try:
        logger.info(
            "Phase 1 search: resolution<=%.2f release_date>%s protein_entities=2..%d nucleic_acid=0",
            config.RESOLUTION_CUTOFF_ANGSTROM,
            config.PDB_RELEASE_DATE_CUTOFF,
            config.MAX_PROTEIN_ENTITIES,
        )
        entry_ids, total_count = search_entries(session, max_entries=max_entries)
        logger.info(
            "Search returned %d total hits; processing %d entries",
            total_count,
            len(entry_ids),
        )
    finally:
        if max_protein_entities is not None:
            config.MAX_PROTEIN_ENTITIES = original

    sifts_path = download_sifts_mapping(session)
    logger.info("Loading SIFTS mapping from %s", sifts_path)
    sifts_mapping = load_sifts_mapping(sifts_path)
    logger.info("Loaded %d SIFTS chain mappings", len(sifts_mapping))

    entry_summaries: dict[str, dict] = {}
    chain_rows: list[dict] = []
    n_entry_fetch_errors = 0
    n_sanity_check_failures = 0

    for entry_id in entry_ids:
        try:
            entry_data = fetch_entry(session, entry_id)
        except requests.RequestException as exc:
            logger.warning("failed to fetch entry %s: %s", entry_id, exc)
            n_entry_fetch_errors += 1
            continue

        summary = parse_entry_summary(entry_id, entry_data)

        # Defensive re-check of the search API's own filters (PLAN.md Phase 1
        # success criterion: 100% of retained rows satisfy the cutoffs).
        resolution = summary["resolution"]
        release_date = summary["release_date"]
        if (
            resolution is None
            or resolution > config.RESOLUTION_CUTOFF_ANGSTROM
            or not release_date
            or release_date <= config.PDB_RELEASE_DATE_CUTOFF
        ):
            logger.warning(
                "entry %s failed local sanity check (resolution=%s release_date=%s) - excluded",
                entry_id,
                resolution,
                release_date,
            )
            n_sanity_check_failures += 1
            continue

        entry_summaries[entry_id] = summary

        for entity_id in summary["polymer_entity_ids"]:
            try:
                entity_data = fetch_polymer_entity(session, entry_id, entity_id)
            except requests.RequestException as exc:
                logger.warning(
                    "failed to fetch polymer entity %s/%s: %s", entry_id, entity_id, exc
                )
                continue

            for chain_row in parse_polymer_entity(entry_id, entity_id, entity_data):
                sifts_ids, agrees = cross_check_uniprot(
                    entry_id, chain_row["chain_id"], chain_row["uniprot_ids"], sifts_mapping
                )
                if not agrees:
                    logger.debug(
                        "SIFTS disagreement %s_%s: rcsb=%s sifts=%s",
                        entry_id,
                        chain_row["chain_id"],
                        chain_row["uniprot_ids"],
                        sifts_ids,
                    )
                chain_row["sifts_uniprot_ids"] = sifts_ids
                chain_row["sifts_agrees"] = agrees
                chain_row["release_date"] = summary["release_date"]
                chain_row["resolution"] = summary["resolution"]
                chain_row["n_protein_entities"] = summary["polymer_entity_count_protein"]
                chain_rows.append(chain_row)

    n_sifts_disagreements = sum(1 for r in chain_rows if not r["sifts_agrees"])
    logger.info(
        "Fetched %d/%d entries successfully (%d fetch errors, %d sanity-check failures)",
        len(entry_summaries),
        len(entry_ids),
        n_entry_fetch_errors,
        n_sanity_check_failures,
    )
    logger.info(
        "Extracted %d raw chain rows; SIFTS disagreed with RCSB on %d/%d chains",
        len(chain_rows),
        n_sifts_disagreements,
        len(chain_rows),
    )

    surviving_rows, attrition = run_filter_pipeline(chain_rows, config.MIN_CHAIN_LENGTH)

    # Finalize output rows (flatten list fields to ';'-joined strings for CSV).
    output_rows = []
    for row in surviving_rows:
        output_rows.append(
            {
                "pdb_id": row["pdb_id"],
                "release_date": row["release_date"],
                "resolution": row["resolution"],
                "n_protein_entities": row["n_protein_entities"],
                "entity_id": row["entity_id"],
                "chain_id": row["chain_id"],
                "seq_length": row["seq_length"],
                "polymer_type": row["polymer_type"],
                "uniprot_ids": ";".join(row["uniprot_ids"]),
                "n_uniprot_ids": len(row["uniprot_ids"]),
                "sifts_uniprot_ids": ";".join(row["sifts_uniprot_ids"]),
                "sifts_agrees": row["sifts_agrees"],
            }
        )

    write_candidates(output_rows, CANDIDATES_PATH)
    write_attrition(attrition, ATTRITION_PATH)

    n_complexes = len({r["pdb_id"] for r in output_rows})
    n_uniprot = len({acc for r in surviving_rows for acc in r["uniprot_ids"]})
    logger.info(
        "Wrote %d candidate chains across %d entries (%d unique UniProt accessions) to %s",
        len(output_rows),
        n_complexes,
        n_uniprot,
        CANDIDATES_PATH,
    )
    logger.info("Wrote attrition table to %s", ATTRITION_PATH)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Phase 1: candidate complex selection")
    parser.add_argument(
        "--max-entries",
        type=int,
        default=None,
        help="Cap on number of search hits to process (omit for a full run)",
    )
    parser.add_argument(
        "--max-protein-entities",
        type=int,
        default=None,
        help="Override config.MAX_PROTEIN_ENTITIES for this run",
    )
    args = parser.parse_args(argv)
    run(max_entries=args.max_entries, max_protein_entities=args.max_protein_entities)


if __name__ == "__main__":
    main()
