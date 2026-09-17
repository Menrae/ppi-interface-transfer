"""Phase 6 (lean): PeSTo inference on isolated single chains.

Scope this pass (see PLAN.md "Scope decision (2026-09-17): lean path for
Phases 6-9"): the primary set (`eligible_homolog`) only, with two input
variants -- the isolated experimental chain (`exp`) and the AlphaFold
model trimmed to the experimental chain's mapped UniProt range
(`af_trimmed`). The `af_full` code path stays working (`--inputs
exp,af_full,af_trimmed`) but isn't part of the default invocation.

Each (chain, input variant) job runs PeSTo in its own OS subprocess
(`pesto_worker.py`, invoked as a plain script -- see that file's docstring
for why), both to isolate an OOM/crash to just that job and to avoid a
`sys.modules["src"]` collision between this project's own `src` package
and PeSTo's own absolutely-imported `src` package.

Runnable as: .venv/bin/python -m src.models.run_pesto --subset primary --inputs exp,af_trimmed
"""

from __future__ import annotations

import argparse
import csv
import logging
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import pandas as pd

from src import config
from src.data.interface_labels import HeavyAtom, resolve_heavy_atoms

logger = logging.getLogger(__name__)

# --- Paths -----------------------------------------------------------------

PESTO_ROOT = config.EXTERNAL_DIR / "PeSTo"
PESTO_CHECKPOINT_DIR = PESTO_ROOT / "model" / "save" / "i_v4_1_2021-09-07_11-21"
PESTO_WORKER_SCRIPT = Path(__file__).resolve().parent / "pesto_worker.py"

LABELS_REPORT_PATH = config.INTERIM_DATA_DIR / "labels_report.csv"
FETCH_REPORT_PATH = config.INTERIM_DATA_DIR / "fetch_report.csv"
RESIDUE_MAPS_DIR = config.INTERIM_DATA_DIR / "residue_mappings"
PDB_UPDATED_DIR = config.RAW_DATA_DIR / "pdb"
ALPHAFOLD_DIR = config.RAW_DATA_DIR / "alphafold"

PESTO_INPUTS_DIR = config.INTERIM_DATA_DIR / "pesto_inputs"
PREDICTIONS_DIR = config.PROCESSED_DATA_DIR / "predictions"
INFERENCE_REPORT_PATH = config.INTERIM_DATA_DIR / "inference_report.csv"
ATTRITION_PATH = config.INTERIM_DATA_DIR / "phase6_attrition.csv"

INPUT_VARIANTS = ("exp", "af_full", "af_trimmed")
DEFAULT_INPUTS = ("exp", "af_trimmed")
SUBSET_TO_MODE = {"primary": "homolog", "homolog": "homolog", "exact_train": "exact_train", "none": "none"}

# --- Failure/skip reasons ---------------------------------------------------

REASON_EXCEEDS_MEMORY = "exceeds_memory_ceiling"
REASON_WORKER_CRASHED = "worker_crashed_or_oom"
REASON_INPUT_PREP_ERROR = "input_prep_error"
REASON_JOIN_INCOMPLETE = "join_incomplete"

# Ideal (fully-modeled) heavy-atom count per standard amino acid --
# textbook values (backbone N/CA/C/O + sidechain). Modified residues are
# resolved to their parent one-letter code (same convention as Phases 4/5)
# and use the parent's count as an approximation (documented as such --
# e.g. MSE is scored against MET's count, since Se simply replaces S).
IDEAL_HEAVY_ATOM_COUNTS = {
    "G": 4, "A": 5, "S": 6, "C": 6, "P": 7, "T": 7, "V": 7,
    "N": 8, "D": 8, "I": 8, "L": 8, "M": 8,
    "E": 9, "Q": 9, "K": 9,
    "H": 10, "F": 11, "R": 11, "Y": 12, "W": 14,
}


def one_letter_code(comp_id: str) -> str | None:
    info = gemmi.find_tabulated_residue(comp_id)
    if info is None or not info.is_amino_acid():
        return None
    return info.one_letter_code.upper()


# --- Data model --------------------------------------------------------


@dataclass
class Job:
    cluster_id: str
    pdb_id: str
    chain_id: str
    uniprot_acc: str
    variant: str
    pdb_path: Path
    keys_path: Path
    output_path: Path
    n_atoms: int = 0
    n_residues: int = 0
    missing_heavy_atoms: dict = field(default_factory=dict)  # {(auth_seq_id, ins_code): n_missing}, exp only
    status: str = "pending"  # pending|success|skipped_resumed|skipped_too_large|failed
    reason: str = ""
    runtime_s: float = 0.0


# --- Input preparation -------------------------------------------------


def _write_pdb(atoms: list[tuple], path: Path) -> None:
    """atoms: list of (atom_name, element, comp_id, chain_letter, resnum, x, y, z)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for i, (name, element, comp_id, chain_letter, resnum, x, y, z) in enumerate(atoms, start=1):
            line = "{:<6s}{:>5d} {:<4s} {:>3s} {:1s}{:>4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:<2s}  ".format(
                "ATOM", i, name, comp_id, chain_letter, resnum, x, y, z, 0.0, 0.0, element
            )
            f.write(line + "\n")
        f.write("TER\nEND\n")


def _write_keys(rows: list[tuple], header: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def prepare_experimental_input(
    pdb_cif_path: Path, chain_id: str, phase4_keys: list[tuple[int, str]], out_pdb: Path, out_keys: Path
) -> tuple[int, dict]:
    """Isolated single chain: same alt-loc rule as Phase 4/5, restricted to
    exactly the Phase 4 observed-residue set (sorted ascending), all other
    chains/waters/ligands/ions removed. Returns (n_atoms, missing_heavy_atoms).

    Raises ValueError if the chain isn't found or its residue set doesn't
    match phase4_keys exactly (never silently proceeds on a mismatch).
    """
    st = gemmi.read_structure(str(pdb_cif_path))
    st.setup_entities()
    chain_names = [c.name for c in st[0]]
    if chain_id not in chain_names:
        raise ValueError(f"chain {chain_id} not found in {pdb_cif_path}")
    chain = st[0][chain_id]

    by_key = {}
    for res in chain:
        if res.label_seq is None:
            continue
        key = (res.seqid.num, res.seqid.icode.strip())
        by_key[key] = res

    wanted = set(phase4_keys)
    found = set(by_key)
    if wanted != found:
        raise ValueError(
            f"residue key mismatch for {pdb_cif_path.name}/{chain_id}: "
            f"missing={len(wanted - found)} extra={len(found - wanted)}"
        )

    atoms = []
    missing_heavy_atoms = {}
    for key in sorted(phase4_keys):
        res = by_key[key]
        heavy = resolve_heavy_atoms(res)
        for atom in heavy:
            atoms.append((atom.name, atom.element, atom.comp_id, chain_id[:1], key[0], atom.pos[0], atom.pos[1], atom.pos[2]))
        letter = one_letter_code(res.name)
        ideal = IDEAL_HEAVY_ATOM_COUNTS.get(letter) if letter else None
        missing_heavy_atoms[key] = max(ideal - len(heavy), 0) if ideal is not None else None

    _write_pdb(atoms, out_pdb)
    _write_keys([[k[0], k[1]] for k in sorted(phase4_keys)], ["auth_seq_id", "auth_ins_code"], out_keys)
    return len(atoms), missing_heavy_atoms


def prepare_alphafold_input(
    af_cif_path: Path, uniprot_range: tuple[int, int] | None, out_pdb: Path, out_keys: Path
) -> int:
    """AlphaFold model, format-converted to PDB. If uniprot_range is given
    (af_trimmed), restrict to that inclusive range (AF's own numbering is
    UniProt-native -- verified in Phase 4); otherwise (af_full) keep the
    whole model. AlphaFold models have no alt-locs/waters/ligands to strip.
    """
    st = gemmi.read_structure(str(af_cif_path))
    st.setup_entities()
    chain_names = [c.name for c in st[0]]
    if len(chain_names) != 1:
        raise ValueError(f"expected a single-chain AlphaFold model, got {chain_names} in {af_cif_path}")
    chain = st[0][chain_names[0]]

    residues = [res for res in chain if res.label_seq is not None]
    if uniprot_range is not None:
        lo, hi = uniprot_range
        residues = [res for res in residues if lo <= res.seqid.num <= hi]
    residues.sort(key=lambda r: r.seqid.num)

    atoms = []
    keys = []
    for res in residues:
        heavy = resolve_heavy_atoms(res)
        for atom in heavy:
            atoms.append((atom.name, atom.element, atom.comp_id, "A", res.seqid.num, atom.pos[0], atom.pos[1], atom.pos[2]))
        keys.append([res.seqid.num])

    _write_pdb(atoms, out_pdb)
    _write_keys(keys, ["uniprot_resnum"], out_keys)
    return len(atoms)


# --- Loading scope -------------------------------------------------------


def _csv_bool(value: str) -> bool:
    return value == "True"


def load_representatives(subset: str) -> list[dict]:
    """Chains labeled in Phase 5, eligible under the requested leakage mode
    (`subset` in {"primary","homolog","exact_train","none"} -- "primary"
    is an alias for "homolog", this lean pass's in-scope mode)."""
    mode = SUBSET_TO_MODE[subset]
    labels = pd.read_csv(LABELS_REPORT_PATH)
    labels = labels[(labels["status"] == "labeled") & (labels[f"eligible_{mode}"])]

    fetch = pd.read_csv(FETCH_REPORT_PATH)[
        ["pdb_id", "chain_id", "uniprot_acc", "pdb_cif_path", "alphafold_cif_path"]
    ]
    merged = labels.merge(fetch, on=["pdb_id", "chain_id", "uniprot_acc"], how="left")
    if merged["pdb_cif_path"].isna().any() or merged["alphafold_cif_path"].isna().any():
        raise ValueError("some labeled representatives are missing fetch_report structure paths")

    return merged[
        ["cluster_id", "pdb_id", "chain_id", "uniprot_acc", "pdb_cif_path", "alphafold_cif_path",
         "eligible_homolog", "eligible_exact_train", "eligible_none"]
    ].to_dict("records")


def load_phase4_keys(pdb_id: str, chain_id: str) -> tuple[list[tuple[int, str]], tuple[int, int]]:
    """Returns (observed residue keys, (min_uniprot_resnum, max_uniprot_resnum))."""
    df = pd.read_parquet(RESIDUE_MAPS_DIR / f"{pdb_id}_{chain_id}.parquet")
    keys = [(int(r.auth_seq_id), r.auth_ins_code or "") for r in df.itertuples()]
    mapped = df[df["uniprot_resnum"].notna()]
    uniprot_range = (int(mapped["uniprot_resnum"].min()), int(mapped["uniprot_resnum"].max()))
    return keys, uniprot_range


# --- Job construction ----------------------------------------------------


def output_path_for(variant: str, pdb_id: str, chain_id: str) -> Path:
    return PREDICTIONS_DIR / variant / f"{pdb_id}_{chain_id}.parquet"


def _count_keys_rows(keys_path: Path) -> int:
    with open(keys_path, newline="") as f:
        return sum(1 for _ in f) - 1  # minus header


def is_valid_output(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        df = pd.read_parquet(path)
    except Exception:  # noqa: BLE001
        return False
    return "pesto_interface_prob" in df.columns and len(df) > 0


def build_jobs(reps: list[dict], inputs: tuple[str, ...]) -> list[Job]:
    jobs = []
    for rep in reps:
        try:
            phase4_keys, uniprot_range = load_phase4_keys(rep["pdb_id"], rep["chain_id"])
        except Exception as exc:  # noqa: BLE001
            logger.warning("skipping %s_%s: could not load Phase 4 keys: %s", rep["pdb_id"], rep["chain_id"], exc)
            continue

        for variant in inputs:
            base = f"{rep['pdb_id']}_{rep['chain_id']}_{variant}"
            job = Job(
                cluster_id=rep["cluster_id"], pdb_id=rep["pdb_id"], chain_id=rep["chain_id"],
                uniprot_acc=rep["uniprot_acc"], variant=variant,
                pdb_path=PESTO_INPUTS_DIR / f"{base}.pdb", keys_path=PESTO_INPUTS_DIR / f"{base}.keys.csv",
                output_path=output_path_for(variant, rep["pdb_id"], rep["chain_id"]),
            )

            if is_valid_output(job.output_path):
                job.status, job.reason = "skipped_resumed", ""
                job.n_residues = len(pd.read_parquet(job.output_path))
                jobs.append(job)
                continue

            try:
                if variant == "exp":
                    n_atoms, missing = prepare_experimental_input(
                        Path(rep["pdb_cif_path"]), rep["chain_id"], phase4_keys, job.pdb_path, job.keys_path
                    )
                    job.missing_heavy_atoms = missing
                    job.n_residues = len(missing)
                elif variant == "af_full":
                    n_atoms = prepare_alphafold_input(Path(rep["alphafold_cif_path"]), None, job.pdb_path, job.keys_path)
                    job.n_residues = _count_keys_rows(job.keys_path)
                elif variant == "af_trimmed":
                    n_atoms = prepare_alphafold_input(
                        Path(rep["alphafold_cif_path"]), uniprot_range, job.pdb_path, job.keys_path
                    )
                    job.n_residues = _count_keys_rows(job.keys_path)
                else:
                    raise ValueError(f"unknown variant {variant}")
            except Exception as exc:  # noqa: BLE001
                job.status, job.reason = "failed", f"{REASON_INPUT_PREP_ERROR}:{exc}"
                jobs.append(job)
                continue

            job.n_atoms = n_atoms
            if n_atoms > config.PESTO_MAX_ATOMS:
                job.status, job.reason = "skipped_too_large", REASON_EXCEEDS_MEMORY
            jobs.append(job)
    return jobs


# --- Execution -----------------------------------------------------------


def run_one_job(job: Job, n_threads: int) -> None:
    """Never raises -- any failure (subprocess crash/OOM/timeout, or an
    unexpected error in this function itself, e.g. augmenting the output)
    is caught and recorded on the job, so one bad job can never take down
    run_tiered's ThreadPoolExecutor or the rest of the batch."""
    t0 = time.time()
    try:
        try:
            result = subprocess.run(
                [
                    sys.executable, str(PESTO_WORKER_SCRIPT), str(PESTO_ROOT), str(PESTO_CHECKPOINT_DIR),
                    str(job.pdb_path), str(job.keys_path), str(job.output_path), str(n_threads),
                ],
                capture_output=True, text=True, timeout=1800,
            )
        except subprocess.TimeoutExpired:
            job.status, job.reason = "failed", f"{REASON_WORKER_CRASHED}:timeout"
            return

        if result.returncode == 0 and is_valid_output(job.output_path):
            job.status = "success"
            if job.variant == "exp" and job.missing_heavy_atoms:
                _augment_exp_output(job)
        else:
            job.status = "failed"
            job.reason = f"{REASON_WORKER_CRASHED}:returncode={result.returncode}"
            logger.warning(
                "job failed %s_%s/%s (returncode=%d): %s",
                job.pdb_id, job.chain_id, job.variant, result.returncode,
                result.stderr[-500:] if result.stderr else "",
            )
    except Exception as exc:  # noqa: BLE001
        job.status, job.reason = "failed", f"{REASON_WORKER_CRASHED}:unexpected_error:{exc}"
        logger.warning("unexpected error running job %s_%s/%s: %s", job.pdb_id, job.chain_id, job.variant, exc)
    finally:
        job.runtime_s = time.time() - t0


def _augment_exp_output(job: Job) -> None:
    df = pd.read_parquet(job.output_path)
    df["n_missing_heavy_atoms"] = [
        job.missing_heavy_atoms.get((int(r.auth_seq_id), r.auth_ins_code or ""))
        for r in df.itertuples()
    ]
    df.to_parquet(job.output_path, index=False)


def run_tiered(jobs: list[Job]) -> None:
    """Runs `pending` jobs in ascending-atom-count tiers, one tier fully
    before the next, per config.PESTO_CONCURRENCY_TIERS -- so memory from
    different tiers is never concurrent, and each tier's own concurrency
    was chosen to stay within this container's memory budget (see
    config.py). A crashed/OOM-killed job is logged and does not stop the
    run -- other jobs in the same tier, and all later tiers, still run.
    """
    pending = sorted((j for j in jobs if j.status == "pending"), key=lambda j: j.n_atoms)
    tier_bounds = config.PESTO_CONCURRENCY_TIERS

    idx = 0
    for max_atoms, n_workers, n_threads in tier_bounds:
        tier_jobs = []
        while idx < len(pending) and pending[idx].n_atoms <= max_atoms:
            tier_jobs.append(pending[idx])
            idx += 1
        if not tier_jobs:
            continue
        logger.info(
            "tier <=%d atoms: %d jobs, %d workers x %d threads", max_atoms, len(tier_jobs), n_workers, n_threads
        )
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(run_one_job, job, n_threads): job for job in tier_jobs}
            for i, future in enumerate(as_completed(futures)):
                future.result()  # run_one_job never raises; mutates job in place
                if (i + 1) % 50 == 0 or (i + 1) == len(tier_jobs):
                    logger.info("tier progress: %d/%d", i + 1, len(tier_jobs))

    if idx < len(pending):
        # Should not happen: every job's n_atoms <= PESTO_MAX_ATOMS <= the
        # last tier's bound by construction (build_jobs already filters
        # anything larger out as skipped_too_large).
        leftover = pending[idx:]
        logger.error("%d jobs did not fit any concurrency tier -- marking failed", len(leftover))
        for job in leftover:
            job.status, job.reason = "failed", "no_matching_concurrency_tier"


# --- Reporting -------------------------------------------------------------


def write_inference_report(jobs: list[Job], path: Path = INFERENCE_REPORT_PATH) -> None:
    by_chain: dict[tuple[str, str], dict] = {}
    for job in jobs:
        key = (job.pdb_id, job.chain_id)
        row = by_chain.setdefault(
            key,
            {"cluster_id": job.cluster_id, "pdb_id": job.pdb_id, "chain_id": job.chain_id, "uniprot_acc": job.uniprot_acc},
        )
        row[f"{job.variant}_status"] = job.status
        row[f"{job.variant}_reason"] = job.reason
        row[f"{job.variant}_n_atoms"] = job.n_atoms
        row[f"{job.variant}_n_residues"] = job.n_residues
        row[f"{job.variant}_runtime_s"] = round(job.runtime_s, 2)

    variants_present = sorted({job.variant for job in jobs})
    fields = ["cluster_id", "pdb_id", "chain_id", "uniprot_acc"]
    for variant in variants_present:
        fields += [f"{variant}_status", f"{variant}_reason", f"{variant}_n_atoms", f"{variant}_n_residues", f"{variant}_runtime_s"]

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, restval="")
        writer.writeheader()
        for row in by_chain.values():
            writer.writerow(row)
    logger.info("Wrote inference report (%d chains) to %s", len(by_chain), path)


def write_attrition(jobs: list[Job], subset: str, path: Path = ATTRITION_PATH) -> None:
    labels = pd.read_csv(LABELS_REPORT_PATH)
    labeled = labels[labels["status"] == "labeled"]
    chains_run = {(j.pdb_id, j.chain_id) for j in jobs}

    rows = []
    for mode in config.LEAKAGE_FILTER_MODES:
        eligible = labeled[labeled[f"eligible_{mode}"]]
        n_eligible = len(eligible)
        eligible_keys = set(zip(eligible["pdb_id"], eligible["chain_id"]))
        run_keys = eligible_keys & chains_run
        n_succeeded = sum(
            1 for key in run_keys
            if all(j.status in ("success", "skipped_resumed") for j in jobs if (j.pdb_id, j.chain_id) == key)
        )
        note = (
            f"lean scope: only the '{subset}' subset was run this pass"
            if mode != SUBSET_TO_MODE[subset]
            else f"meets_target({config.MIN_TEST_CHAINS_TARGET})={n_succeeded >= config.MIN_TEST_CHAINS_TARGET}"
        )
        rows.append(
            {
                "leakage_filter_mode": mode, "n_eligible_representatives": n_eligible,
                "n_run": len(run_keys), "n_succeeded": n_succeeded, "note": note,
            }
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["leakage_filter_mode", "n_eligible_representatives", "n_run", "n_succeeded", "note"])
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote attrition table to %s", path)


def verify_joins(jobs: list[Job]) -> list[tuple[str, str, str]]:
    """Item 2's assert: every Phase-4-mapped residue must have a prediction
    in every successfully-run variant. Returns a list of (pdb_id, chain_id,
    problem) for any chain that fails -- expected to be empty."""
    problems = []
    by_chain: dict[tuple[str, str], list[Job]] = {}
    for job in jobs:
        by_chain.setdefault((job.pdb_id, job.chain_id), []).append(job)

    for (pdb_id, chain_id), chain_jobs in by_chain.items():
        successful = [j for j in chain_jobs if j.status in ("success", "skipped_resumed")]
        if not successful:
            continue
        mapped = pd.read_parquet(RESIDUE_MAPS_DIR / f"{pdb_id}_{chain_id}.parquet")
        mapped = mapped[mapped["uniprot_resnum"].notna()]
        for job in successful:
            df = pd.read_parquet(job.output_path)
            if job.variant == "exp":
                have = set(zip(df["auth_seq_id"], df["auth_ins_code"].fillna("")))
                want = set(zip(mapped["auth_seq_id"], mapped["auth_ins_code"].fillna("")))
            else:
                have = set(df["uniprot_resnum"])
                want = set(mapped["uniprot_resnum"].astype(int))
            if not want.issubset(have):
                problems.append((pdb_id, chain_id, f"{job.variant}: missing {len(want - have)} mapped residues"))
    return problems


# --- Orchestration -----------------------------------------------------


def _setup_logging() -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = config.LOGS_DIR / "run_pesto.log"

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


def run(subset: str = "primary", inputs: tuple[str, ...] = DEFAULT_INPUTS, limit: int | None = None) -> list[Job]:
    _setup_logging()
    reps = load_representatives(subset)
    if limit is not None:
        reps = reps[:limit]
    logger.info("Scope: subset=%s inputs=%s -> %d representatives", subset, inputs, len(reps))

    jobs = build_jobs(reps, inputs)
    n_by_status = {}
    for job in jobs:
        n_by_status[job.status] = n_by_status.get(job.status, 0) + 1
    logger.info("Built %d jobs: %s", len(jobs), n_by_status)

    run_tiered(jobs)

    n_by_status = {}
    for job in jobs:
        n_by_status[job.status] = n_by_status.get(job.status, 0) + 1
    logger.info("Final job status counts: %s", n_by_status)

    problems = verify_joins(jobs)
    if problems:
        logger.error("Join verification found %d problems: %s", len(problems), problems[:10])
    else:
        logger.info("Join verification: all successful outputs join exactly to Phase 4's mapped residues")

    write_inference_report(jobs)
    write_attrition(jobs, subset)
    return jobs


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Phase 6: PeSTo inference on isolated single chains")
    parser.add_argument("--subset", choices=list(SUBSET_TO_MODE), default="primary")
    parser.add_argument("--inputs", default=",".join(DEFAULT_INPUTS), help="comma-separated subset of exp,af_full,af_trimmed")
    parser.add_argument("--limit", type=int, default=None, help="process only the first N representatives (smoke testing)")
    args = parser.parse_args(argv)
    inputs = tuple(args.inputs.split(","))
    for v in inputs:
        if v not in INPUT_VARIANTS:
            parser.error(f"unknown input variant {v!r}, must be one of {INPUT_VARIANTS}")
    run(subset=args.subset, inputs=inputs, limit=args.limit)


if __name__ == "__main__":
    main()
