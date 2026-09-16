# Progress

Standing status file. Read this first in any new session on this project —
it says what's actually done, what's decided, and what to do next. Update
it whenever you finish a chunk of work or make a decision that changes the
plan. `PLAN.md` is the design document (10 phases, confounds, deviations);
this file tracks execution against it.

Last updated: 2026-09-16.

## Where things stand

**Phases 1 and 2 are done and have run successfully.** Phases 3-10 are not
started. See "Completed phases" below for what each produced and what it
means for the next phase. Phase 2's leakage-sweep result (see below) is
important context before starting Phase 3: the current candidate pool is
too small and too leakage-heavy for the primary benchmark, and a full-scale
Phase 1 rerun is recommended before investing further downstream.

**Phase 1's code/outputs were committed and pushed** to
`https://github.com/Menrae/ppi-interface-transfer.git` (`main`, commit
`bd28a1b`). **Phase 2's code and outputs (this session's work) are not yet
committed** — confirm with the user before committing/pushing, per their
standing preference to review first; don't assume a green light carries
over between sessions.

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
  `CLUSTER_MIN_COVERAGE=0.80`, `LEAKAGE_FILTER_MODES=("homolog","exact_train","none")`,
  `MIN_TEST_CHAINS_TARGET=100`, `MIN_TEST_CHAINS_FLOOR=50`,
  `LEAKAGE_SWEEP_IDENTITY_THRESHOLDS=(0.30,0.50,0.70,0.95)`.
- `src/data/select_complexes.py` — Phase 1 implementation (see "Completed
  phases" below). `tests/test_select_complexes.py` (11 tests, passing) +
  `tests/fixtures/` (small saved RCSB/SIFTS JSON/TSV fixtures, no network
  needed to run the tests).
- `src/data/cluster_and_split.py` — Phase 2 implementation (see "Completed
  phases" below). `tests/test_cluster_and_split.py` (14 tests, passing,
  including one real end-to-end MMseqs2 invocation on tiny synthetic
  sequences — skipped automatically if the binary is absent).
- `docs/Project_Proposal.pdf` + `docs/proposal.txt` — the original proposal
  and its extracted text (via `pypdf`, now in `requirements.txt`).
- `README.md` — rewritten for a general-science-background public audience
  (plain-language PPI/AlphaFold/PeSTo/pLDDT background, phase-by-phase
  approach, current status, repro instructions, key references w/ DOIs).
  Update its "Current status" section after each phase, same as this file.
- `external/mmseqs/bin/mmseqs` — **installed 2026-09-16**: official static
  AVX2 Linux build (commit `d401e78c2d18a822cdb1527d7464a043f6035a15`),
  gitignored. The AVX2 build was used (not SSE4.1) since `/proc/cpuinfo`
  confirmed AVX2 support in this container.
- `data/raw/pesto_splits/` — **already downloaded and cached**: PeSTo's
  three dataset-split files from `github.com/LBM-EPFL/PeSTo/data/datasets/`
  (`subunits_train_set.txt` 376,216 lines, `subunits_test_set.txt` 97,424
  lines, `subunits_validation_set.txt` 101,700 lines, `bc-30.out`,
  `README.md`). Note: this `README.md` itself states "the definition of
  test and validation set is swapped in this source code compared to the
  commonly used definition" — independent confirmation of the train-code-
  vs-paper discrepancy already recorded in `PLAN.md`'s evidence table.
  `data/raw/pesto_splits/pesto_all_splits.fasta` (560,829 sequences, Phase 2
  output) is also cached here now.
- `data/raw/pdb_seqres/pdb_seqres.txt.gz` — **new, Phase 2**: bulk
  PDBID_CHAIN→sequence FASTA for the entire PDB (~1.16M records, ~67 MB),
  from `files.wwpdb.org`. Used instead of per-chain RCSB calls to resolve
  sequences for PeSTo's ~575k split-file chains in one download.
- `.venv/` — Python 3.11 venv with everything in `requirements.txt`
  installed (biotite, biopython, gemmi, numpy/scipy/pandas/pyarrow,
  scikit-learn, matplotlib, seaborn, requests, tqdm, pytest, torch+cpu,
  pypdf). **Not yet installed:** `statsmodels` (needed for Phase 8's VIF
  check — install and re-`pip freeze` when starting that phase).
- `data/interim/candidates.csv` (1275 rows), `data/interim/phase1_attrition.csv`
  — Phase 1 outputs, see "Completed phases" below.
- `data/interim/{sequences.fasta,clusters.tsv,candidates_dedup.csv,
  pesto_homology_search.tsv,leakage_threshold_sweep.csv,phase2_attrition.csv}`
  — Phase 2 outputs, see "Completed phases" below.
- `data/raw/rcsb/{search,entries,polymer_entities}/` — cached RCSB API
  responses from the Phase 1 run (500-entry cap); reruns of
  `select_complexes` against the same query/entries are fully offline.
- `data/raw/sifts/pdb_chain_uniprot.tsv.gz` — cached bulk SIFTS file
  (~986k chain mappings loaded from it).
- `data/raw/mmseqs_work/` — MMseqs2 intermediate cluster/search DBs and tmp
  dirs from the Phase 2 run (cache, not a curated output).
- `logs/select_complexes.log`, `logs/cluster_and_split.log` — full DEBUG
  logs of the Phase 1 and Phase 2 runs respectively.
- Empty scaffold directories: `data/processed`, `results/`, `notebooks/`
  (PeSTo's model repo itself has not been checked out here yet — only its
  split-list files were fetched for investigation; that checkout is a
  Phase 6 input).

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
- **DSSP/freesasa are not installed** and there's no root access in this
  container — fallbacks are `biotite.structure.annotate_sse` (DSSP) and
  `biotite.structure.sasa` (freesasa). **MMseqs2 is now installed**
  (`external/mmseqs/bin/mmseqs`, static AVX2 build) — see "What exists on
  disk" above; the Biopython-pairwise fallback was not needed.
- **Cluster representative selection (Phase 2) is deterministic and
  independent of MMseqs2's own internal representative choice:** best
  (lowest) resolution, then longest chain, then ascending `PDBID_CHAINID`
  string. MMseqs2's `easy-cluster` is used only to assign cluster
  membership; `src/data/cluster_and_split.py` picks the actual
  representative itself so the choice is reproducible and documented.
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

### Phase 2 — Redundancy reduction + PeSTo-overlap flagging (2026-09-16)

`src/data/cluster_and_split.py`, runnable as
`.venv/bin/python -m src.data.cluster_and_split`. Run against Phase 1's
500-entry-slice candidate pool (1275 chains) — **not yet run at full
scale**, see "Next step" below. See `PLAN.md` Phase 2 for the
module/output-path reorganization vs. the original sketch.

**Attrition (chain-level, from `data/interim/phase2_attrition.csv`):**

| stage | n_before | n_dropped | n_remaining |
| --- | --- | --- | --- |
| initial_candidate_chains | 1275 | 0 | 1275 |
| missing_cached_sequence | 1275 | 0 | 1275 |
| redundancy_clustering | 1275 | 1129 | 146 |

All 1275 candidate chains had a cached canonical sequence (Phase 1's
polymer-entity JSON cache was complete — 0 dropped). MMseqs2 `easy-cluster`
(`--min-seq-id 0.30 -c 0.80 --cov-mode 0`) collapsed 1275 chains into **146
clusters**; the largest clusters (204 and 184 members) are exactly the two
fragment-screening campaigns flagged as a Phase 1 finding above (yeast
Prp8–Aar2, bovine tubulin), confirming that finding. 135 of 146
representatives' UniProt accessions are unique (close to the full 146,
since each cluster is essentially one biological pair).

**PeSTo split sequence resolution:** loaded all 376,216 train / 97,425 test
/ 101,701 validation chains (line counts match `PLAN.md`'s evidence table);
union = 575,342 distinct chains, of which 560,829 (97.5%) had a sequence in
`pdb_seqres.txt.gz` (the other 14,513 are presumably obsolete/superseded
PDB IDs — logged, excluded from the search target, not silently merged in).
**Bug found and fixed during this run:** the first pass through
`load_pesto_split_membership` uppercased the whole `PDBID_CHAINID` string
before deduplicating, which silently collapsed case-distinct chains (PDB
chain IDs are case-sensitive — e.g. chain `C` and chain `c` can be
genuinely different chains in the same entry) — this undercounted the
union by ~27k entries (548,654 vs. the correct 575,342). Fixed to only
uppercase the PDB ID half of each line; `tests/test_cluster_and_split.py`
now has a regression test
(`test_load_pesto_split_membership_preserves_chain_id_case`). The fix did
**not** change the leakage-flag results below (same 131/146 and 39/146,
same sweep table) — the case-collapsed entries happened not to affect any
representative's best hit — but don't assume that'll always be true; this
is the kind of bug that stays invisible until the target file is large
enough for it to matter.

**Leakage search result — the important finding:** of the 146
representatives, **131 (89.7%) are flagged `pesto_homolog_overlap`** (≥30%
identity to the train+test+validation union) and 39 (26.7%) are flagged
`pesto_exact_train_overlap` (literal ID match to the training file alone).
The gap between those two numbers is exactly what motivated using
sequence-level rather than exact-ID leakage control in the first place —
most of the true overlap would have been invisible to ID matching alone.

**Leakage threshold sweep (`data/interim/leakage_threshold_sweep.csv`),
against `MIN_TEST_CHAINS_TARGET=100` / `MIN_TEST_CHAINS_FLOOR=50`:**

| mode | n_survive | meets target (100) | meets floor (50) |
| --- | --- | --- | --- |
| identity < 0.30 (primary) | 15 | No | No |
| identity < 0.50 | 25 | No | No |
| identity < 0.70 | 27 | No | No |
| identity < 0.95 | 30 | No | No |
| exact-ID (train only) | 107 | Yes | Yes |

**Implication, per `PLAN.md` confound (e):** on this candidate pool, the
primary (`"homolog"`) mode is well under even the 50-chain hard floor, and
stays under it all the way out to 95% identity — this is not a borderline
case fixable by loosening the cutoff slightly. Per the pre-agreed
mitigation, the leakage filter should **not** be loosened to hit the
target; Phase 7 should stratify by `pesto_homolog_overlap` status and lean
on the `"exact_train"`/`"none"` sensitivity modes, which do clear the
target. However, the more likely explanation is simply that **this
candidate pool is a small, non-random 500-entry test slice** (see Phase 1's
fragment-screening-campaign finding) — a full-scale Phase 1 run should
produce enough non-redundant, non-overlapping clusters to reassess this
before concluding the primary analysis is structurally underpowered.

## Next step

**Rerun Phase 1 at full scale first**, then rerun Phase 2 on the result.
`select_complexes` reported 10,774 total hits for the full filter when last
run with `--max-entries 500`; rerunning with a much larger `--max-entries`
(or none) will take longer (more RCSB API calls, all cached/resumable) but
is needed before the Phase 2 leakage-sweep numbers above can be trusted as
representative — the current 15-chain primary-mode survivor count is
plausibly an artifact of the small, campaign-heavy 500-entry slice rather
than a true property of the full candidate pool. Re-run
`.venv/bin/python -m src.data.cluster_and_split` afterward (it will
re-cluster from scratch since `data/interim/candidates.csv` will have
changed; the PeSTo-splits/`pdb_seqres.txt.gz` caches are reusable as-is).

Once a representative-scale Phase 2 result exists, proceed to **Phase 3 —
Download PDB mmCIFs and AlphaFold DB models** (`PLAN.md` Phase 3):
`src/data/download_pdb.py`, `src/data/download_alphafold.py`,
`src/pipeline/download_structures.py`, reading from
`data/interim/candidates_dedup.csv`.

Open question #3 from `PLAN.md` §5 (the Phase 9 significance threshold)
is still unresolved and still cheap to decide now.

## How to update this file

When you finish a phase, a sub-step, or make a new decision that isn't yet
in `PLAN.md`: update "What exists on disk," move finished items out of
"Next step" and into "Decisions already made" or a new "Completed phases"
section, and bump "Last updated." Keep entries about *state and decisions*,
not a narrated log of tool calls.
