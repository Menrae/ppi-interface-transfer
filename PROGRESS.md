# Progress

Standing status file. Read this first in any new session on this project —
it says what's actually done, what's decided, and what to do next. Update
it whenever you finish a chunk of work or make a decision that changes the
plan. `PLAN.md` is the design document (10 phases, confounds, deviations);
this file tracks execution against it.

Last updated: 2026-09-16.

## Where things stand

**Phase 1 (candidate complex selection) is done and has run successfully.**
Phases 2-10 are not started. See "Completed phases" below for what Phase 1
actually produced and what it means for Phase 2.

**Everything has been committed and pushed** to
`https://github.com/Menrae/ppi-interface-transfer.git` (`main`), through
the commit that set up scaffolding/environment/plan. Phase 1's code and
outputs (this session's work) are **not yet committed** — confirm with the
user before committing/pushing, per their standing preference to review
first; don't assume a green light carries over between sessions.

## What exists on disk right now

- `README.md`, `CLAUDE.md` — project description and dev conventions
  (environment facts, no-magic-numbers rule, mandatory tests, no-silent-drop
  logging rule, offline-caching rule). Read `CLAUDE.md` before writing any code.
- `PLAN.md` — the 10-phase implementation plan, confounds/mitigations,
  deviations from the proposal, open questions, and risks. This is the
  authoritative design doc; **implement in phase order**, Phase 1 first.
- `src/config.py` — all paths and tunable thresholds as constants:
  `RESOLUTION_CUTOFF_ANGSTROM=2.5`, `PDB_RELEASE_DATE_CUTOFF="2018-04-30"`
  (release date, not deposit date), `MAX_PROTEIN_ENTITIES=10`,
  `MIN_CHAIN_LENGTH=40`, `PLDDT_BANDS=(50,70,90)`,
  `INTERFACE_DISTANCE_CUTOFF_ANGSTROM=5.0`, `SEQUENCE_IDENTITY_CUTOFF=0.30`,
  `LEAKAGE_FILTER_MODES=("homolog","exact_train","none")`.
- `src/data/select_complexes.py` — Phase 1 implementation (see "Completed
  phases" below). `tests/test_select_complexes.py` (11 tests, passing) +
  `tests/fixtures/` (small saved RCSB/SIFTS JSON/TSV fixtures, no network
  needed to run the tests).
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
- `data/interim/candidates.csv` (1275 rows), `data/interim/phase1_attrition.csv`
  — Phase 1 outputs, see "Completed phases" below.
- `data/raw/rcsb/{search,entries,polymer_entities}/` — cached RCSB API
  responses from the Phase 1 run (500-entry cap); reruns of
  `select_complexes` against the same query/entries are fully offline.
- `data/raw/sifts/pdb_chain_uniprot.tsv.gz` — cached bulk SIFTS file
  (~986k chain mappings loaded from it).
- `logs/select_complexes.log` — full DEBUG log of the Phase 1 run (every
  dropped chain + every SIFTS disagreement, with reasons).
- Empty scaffold directories: `data/processed`, `results/`, `notebooks/`,
  `external/` (PeSTo repo itself has not been checked out here yet, only
  its split-list files were fetched for investigation).

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
- **`.gitignore` gotcha, fixed 2026-09-16:** the original patterns
  (`data/`, `results/`, `logs/`, `external/`, no leading slash) matched
  those directory names at *any* depth, not just repo root — so
  `src/data/` was being silently git-ignored the moment it was created,
  and `select_complexes.py` never even showed up in `git status`. Fixed by
  anchoring to `/data/`, `/results/`, `/logs/`, `/external/`. If a new
  top-level `src/<name>/` package ever collides with a root-level ignored
  dir name again, check `git check-ignore -v <path>` before assuming
  `git status` is complete — it can be silently empty.

## Completed phases

### Phase 1 — Candidate complex selection (2026-09-16)

`src/data/select_complexes.py`, runnable as
`.venv/bin/python -m src.data.select_complexes --max-entries N`. First live
run used `--max-entries 500` (of 10,774 total hits for the full filter —
the full run has not been done yet). See `PLAN.md` Phase 1 for the
file-path/module reorganization vs. the original sketch.

**Attrition (chain-level, from `data/interim/phase1_attrition.csv`):**

| stage | n_before | n_dropped | n_remaining |
| --- | --- | --- | --- |
| initial_protein_chains | 1655 | 0 | 1655 |
| non_protein_entity_type | 1655 | 0 | 1655 |
| min_chain_length>=40 | 1655 | 245 | 1410 |
| uniprot_mapping | 1410 | 132 | 1278 |
| chimera | 1278 | 3 | 1275 |

Entry-level: 500/500 entries fetched successfully (0 fetch errors, 0 local
sanity-check failures — every fetched entry genuinely satisfied
resolution/release-date cutoffs). Final: **1275 candidate chains across 473
entries, 160 unique UniProt accessions**.

**Distributions:** resolution mean 1.94 Å (range 1.00-2.50, as expected
given the cutoff); release years 2018 (196), 2019 (23), 2020 (275), 2021
(300), 2022 (144), 2026 (337) — lumpy because this 500-entry slice is the
alphabetically-first `rcsb_id`-sorted subset of hits, not a random/temporal
sample; a full run will smooth this out. `sifts_agrees` was True for
1265/1275 chains (99.2%); **all 10 disagreements were SIFTS-bulk-file
coverage gaps** (SIFTS had no mapping at all for that chain, not a
conflicting accession) — mostly newer TrEMBL-style accessions
(`A0A...`) not yet in the static bulk snapshot, plus one real antibody
chain (P01854, IgE heavy constant region, PDB 30AF). This validates using
RCSB's own `uniprot_ids` (SIFTS-derived, but live-computed) as the
authoritative source, with the bulk file as a secondary check, as
implemented.

**Important finding for Phase 2:** UniProt diversity is much lower than
entry count suggests (160 accessions / 473 entries) because this slice
contains at least two large **crystallographic fragment-screening
campaigns** deposited as many near-identical entries of the same complex:
184 chains (92 entries) of yeast Prp8–Aar2 (PDB `5QY*`, UniProt
P33334/P32357) and 102 chains (51 entries) of bovine tubulin α/β (PDB
`5S4*`, UniProt P81947/Q6B856) — together ~30% of all entries in this
slice. These will cluster hard in Phase 2's 30%-identity redundancy
reduction (near-100% identity within each campaign), which is the correct
behavior, but it means **raw entry/chain counts from Phase 1 substantially
overstate biological diversity** — don't use them as a proxy for Phase 7's
target sample size without going through Phase 2 first. No antibody-chain
dominance was found in this slice (only the one incidental IgE hit above).

## Next step

**Implement Phase 2** (`PLAN.md` "Phase 2 — Redundancy reduction +
PeSTo-overlap flagging"): build `src/data/sequences.py`,
`src/pipeline/redundancy_reduction.py`, `src/data/pesto_overlap.py` plus
tests. Needs: (1) sequences for the 1275 candidate chains (RCSB
polymer-entity API, or reuse `entity_poly.pdbx_seq_one_letter_code_can`
already present in the cached `data/raw/rcsb/polymer_entities/*.json` from
Phase 1 — check there before re-fetching); (2) a static MMseqs2 binary
fetched into `external/mmseqs/` (not yet done — no root/apt in this
container, see `PLAN.md` §6 Risks); (3) sequences for PeSTo's
train+test+validation chains (the split files are already cached at
`data/raw/pesto_splits/`, but their sequences still need to be fetched).

Before running Phase 2 at scale, consider running Phase 1 without
`--max-entries` (or with a larger cap) first — the 500-entry test run is
a small, non-random slice (see fragment-screening finding above), so
Phase 2's redundancy numbers on just this slice won't be representative.

Open question #3 from `PLAN.md` §5 (the Phase 9 significance threshold)
is still unresolved and still cheap to decide now.

## How to update this file

When you finish a phase, a sub-step, or make a new decision that isn't yet
in `PLAN.md`: update "What exists on disk," move finished items out of
"Next step" and into "Decisions already made" or a new "Completed phases"
section, and bump "Last updated." Keep entries about *state and decisions*,
not a narrated log of tool calls.
