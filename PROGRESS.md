# Progress

Standing status file. Read this first in any new session on this project —
it says what's actually done, what's decided, and what to do next. Update
it whenever you finish a chunk of work or make a decision that changes the
plan. `PLAN.md` is the design document (10 phases, confounds, deviations);
this file tracks execution against it.

Last updated: 2026-09-16.

## Where things stand

**No analysis code has been written yet.** Everything so far is
environment setup, planning, and one piece of investigation (PeSTo split
resolution). Phase 1 of `PLAN.md` (candidate complex selection) has not
been started.

**Nothing has been committed to git.** The user reviews and commits
themselves — do not `git commit` or `git push` unless explicitly asked to
in that session, even if a task appears complete.

## What exists on disk right now

- `README.md`, `CLAUDE.md` — project description and dev conventions
  (environment facts, no-magic-numbers rule, mandatory tests, no-silent-drop
  logging rule, offline-caching rule). Read `CLAUDE.md` before writing any code.
- `PLAN.md` — the 10-phase implementation plan, confounds/mitigations,
  deviations from the proposal, open questions, and risks. This is the
  authoritative design doc; **implement in phase order**, Phase 1 first.
- `src/config.py` — all paths and tunable thresholds as constants:
  `RESOLUTION_CUTOFF_ANGSTROM=2.5`, `PDB_RELEASE_DATE_CUTOFF="2018-04-30"`
  (release date, not deposit date), `PLDDT_BANDS=(50,70,90)`,
  `INTERFACE_DISTANCE_CUTOFF_ANGSTROM=5.0`, `SEQUENCE_IDENTITY_CUTOFF=0.30`,
  `LEAKAGE_FILTER_MODES=("homolog","exact_train","none")`. Only
  `src/__init__.py` and `config.py` exist under `src/` so far — no other
  modules have been written.
- `docs/Project_Proposal.pdf` + `docs/proposal.txt` — the original proposal
  and its extracted text (via `pypdf`, now in `requirements.txt`).
- `data/raw/pesto_splits/` — **already downloaded and cached**: PeSTo's
  three dataset-split files from `github.com/LBM-EPFL/PeSTo/data/datasets/`
  (`subunits_train_set.txt` 376,216 lines, `subunits_test_set.txt` 97,424
  lines, `subunits_validation_set.txt` 101,700 lines, `bc-30.out`,
  `README.md`). These are real cached files, not placeholders — Phase 2 can
  read them directly instead of re-downloading.
- `.venv/` — Python 3.11 venv with everything in `requirements.txt`
  installed (biotite, biopython, gemmi, numpy/scipy/pandas/pyarrow,
  scikit-learn, matplotlib, seaborn, requests, tqdm, pytest, torch+cpu,
  pypdf). **Not yet installed:** `statsmodels` (needed for Phase 8's VIF
  check — install and re-`pip freeze` when starting that phase), MMseqs2
  (needed for Phase 2's homology search — no system package manager access
  in this container; plan is a static binary in `external/mmseqs/`, not
  yet fetched).
- Empty scaffold directories: `data/{raw,interim,processed}` (beyond
  `pesto_splits`), `results/`, `notebooks/`, `tests/` (just `__init__.py`),
  `logs/`, `external/` (empty — PeSTo repo itself has not been checked out
  here yet, only its split-list files were fetched for investigation).

## Decisions already made (don't re-derive these — see `PLAN.md` for full reasoning)

- **AF2 memorization control:** candidates are filtered to
  `initial_release_date > PDB_RELEASE_DATE_CUTOFF` (RCSB release date, not
  deposit date).
- **PeSTo leakage control is sequence-level, not exact-ID.** Verified by
  reading PeSTo's actual training code
  (`model/main.py`/`model/config.py` on GitHub) against the paper's Methods
  text: the training file is unambiguous
  (`subunits_train_set.txt`), but the paper and the code **disagree** about
  which of `subunits_test_set.txt`/`subunits_validation_set.txt` was the
  true untouched holdout (code shows `test_set.txt` was actually used for
  checkpoint selection during training, contradicting the paper). Resolution:
  both are treated as excluded; the homology search runs against the union
  of all three files. Don't re-litigate this without new evidence — see
  `PLAN.md` Phase 2 for the full evidence table.
- **Three-way leakage sensitivity reporting** in Phase 7, controlled by
  `LEAKAGE_FILTER_MODES`: primary = `"homolog"` (exclude sequence-identity
  hits against train+test+validation union), sensitivity A =
  `"exact_train"` (exact-ID match to training file only), sensitivity B =
  `"none"`.
- **Sample-size floor:** target ≥100 non-overlap chains for the primary
  benchmark; 50 is a hard floor. Below 50, stratify by overlap status
  instead of loosening the filter.
- **Apo/unbound arm is best-effort**, reported with its own sample size,
  never blocks the main experimental-vs-AlphaFold benchmark.
- **Phase 9 (fine-tuning) is gated**: only runs if Phase 7 shows a
  significant AlphaFold-vs-experimental gap (primary mode) AND Phase 8 ties
  it to pLDDT. Training script must stay environment-agnostic
  (`--device` flag; Colab vs. HPC undecided).
- **Interfaces computed directly from coordinates** (5 Å heavy-atom +
  ΔSASA), not via PISA — deviation from the proposal, documented in
  `PLAN.md` §4.
- **PeSTo's "no architectural changes" claim is wrong** — the input
  embedding layer must be widened for the pLDDT channel (Phase 9). The
  Phase 9 CPU smoke test must use PeSTo's real released checkpoint, not a
  random init, to verify the widened layer reproduces original outputs.
- **DSSP/freesasa/MMseqs2 are not installed** and there's no root access in
  this container — fallbacks are `biotite.structure.annotate_sse` (DSSP),
  `biotite.structure.sasa` (freesasa), and a static MMseqs2 binary in
  `external/` (MMseqs2 fallback not yet fetched).

## Next step

**Start implementing Phase 1** (`PLAN.md` "Phase 1 — Candidate complex
selection"): build `src/data/rcsb_search.py`, `src/data/sifts.py`,
`src/pipeline/select_candidates.py`, plus their tests, to query the RCSB
Search API for X-ray structures meeting `RESOLUTION_CUTOFF_ANGSTROM` and
`PDB_RELEASE_DATE_CUTOFF`, with ≥2 protein entities and no nucleic acids,
joined against the bulk SIFTS `pdb_chain_uniprot.tsv.gz` for chain-level
UniProt mapping. Follow `PLAN.md`'s exact inputs/outputs/file paths for
that phase, and `CLAUDE.md`'s logging/testing conventions.

Before starting, it's worth resolving open question #3 from `PLAN.md` §5
(the Phase 9 significance threshold) since it's cheap to decide now and
avoids picking a threshold after seeing results later.

## How to update this file

When you finish a phase, a sub-step, or make a new decision that isn't yet
in `PLAN.md`: update "What exists on disk," move finished items out of
"Next step" and into "Decisions already made" or a new "Completed phases"
section, and bump "Last updated." Keep entries about *state and decisions*,
not a narrated log of tool calls.
