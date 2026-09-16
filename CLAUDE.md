# CLAUDE.md

Conventions for working in this repo. Read this before writing code here.
Read `PROGRESS.md` first, though — it says what's actually been done and
what to do next; this file is just standing conventions.

## Environment

- Python 3.11 virtualenv already exists at `.venv` and is expected to be
  active in the shell. Do not create another venv, and do not use `uv` or
  `conda` for this project. Install packages with `.venv/bin/pip`; run
  Python and pytest via `.venv/bin/python` (e.g. `.venv/bin/python -m
  pytest`) so commands never land in another environment.
- The shell also has a conda env called `astro` underneath, and unrelated
  projects live alongside this one in `/workspace` (e.g.
  `imbh-galactic-nuclei`, `jwst-sed-anomaly`). Never read, modify, or
  install into anything outside `/workspace/ppi-interface-transfer`.
- git version in this environment is 2.25.1, older than 2.28 — avoid
  newer git options such as `git init -b` or `git switch`; use
  `git checkout -b` etc. instead.
- There is no GPU. Only install/use the CPU-only build of torch (from
  `https://download.pytorch.org/whl/cpu`).
- `.gitignore` should be appended to, not overwritten, when new ignore
  patterns are needed.

## Code conventions

- All paths and tunable thresholds (resolution cutoffs, date cutoffs,
  pLDDT bands, distance cutoffs, sequence identity cutoffs, etc.) live in
  `src/config.py` as module-level constants. No magic numbers or
  hardcoded paths anywhere else in the codebase — import from
  `src.config` instead.
- Every data-transforming module (parsing, filtering, feature
  extraction, dataset construction, etc.) gets corresponding pytest
  tests under `tests/`.
- Never silently drop residues, chains, structures, or dataset entries.
  Any filtering step must log the count dropped and the reason (e.g.
  "dropped 12 structures: resolution > RESOLUTION_CUTOFF_ANGSTROM").
- Cache all downloaded data (PDB/AlphaFold structures, sequence
  databases, etc.) under `data/raw/` or `external/` so that reruns work
  offline without re-fetching.
