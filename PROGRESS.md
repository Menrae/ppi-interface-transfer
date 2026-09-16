# Progress

Standing status file. Read this first in any new session on this project —
it says what's actually done, what's decided, and what to do next. Update
it whenever you finish a chunk of work or make a decision that changes the
plan. `PLAN.md` is the design document (10 phases, confounds, deviations);
this file tracks execution against it.

Last updated: 2026-09-16.

## Where things stand

**Phases 1 and 2 are done, and both have now been run at full scale**
(uncapped Phase 1 search, all 10,774 hits; Phase 2 unchanged on the result).
Phases 3-10 are not started -- **do not start Phase 3 without an explicit
go-ahead**: the leakage-threshold decision below is the user's to make, not
pre-decided by this pipeline. See "Completed phases" for the full-scale
numbers and "Next step" for exactly what's pending.

The earlier 500-entry pilot run (Phase 1 + Phase 2) is preserved at
`data/interim/pilot_500/` for before/after comparison, not overwritten.
The pilot's leakage-sweep shortfall (primary mode well under
`MIN_TEST_CHAINS_FLOOR`) turned out to be a small/non-representative-sample
artifact, not a structural problem -- the full-scale run clears both
`MIN_TEST_CHAINS_TARGET` and `MIN_TEST_CHAINS_FLOOR` at every swept
threshold.

**Phase 1 and Phase 2's first-run code/outputs were committed and pushed**
to `https://github.com/Menrae/ppi-interface-transfer.git` (`main`, commits
`bd28a1b` and `f7b5c81`). **This session's changes (full-scale rerun code +
outputs) are not yet committed** — confirm with the user before
committing/pushing, per their standing preference to review first; don't
assume a green light carries over between sessions.

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
  phases" below). Now includes RCSB GraphQL-batched entry/entity fetch (with
  REST single-item fallback) and deposition-group collapse, added this
  session for the full-scale rerun. `tests/test_select_complexes.py` (21
  tests, passing) + `tests/fixtures/` (small saved RCSB/SIFTS/GraphQL
  JSON/TSV fixtures, no network needed to run the tests).
- `src/data/cluster_and_split.py` — Phase 2 implementation, **unchanged
  this session** (rerun as-is on the full candidate set, per this session's
  brief). `tests/test_cluster_and_split.py` (15 tests, passing, including
  one real end-to-end MMseqs2 invocation on tiny synthetic sequences —
  skipped automatically if the binary is absent).
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
- `data/interim/candidates.csv` (22,887 rows), `data/interim/phase1_attrition.csv`
  — **full-scale** Phase 1 outputs (this session), see "Completed phases"
  below.
- `data/interim/{sequences.fasta,clusters.tsv,candidates_dedup.csv,
  pesto_homology_search.tsv,leakage_threshold_sweep.csv,phase2_attrition.csv}`
  — **full-scale** Phase 2 outputs (this session), see "Completed phases"
  below.
- `data/interim/leakage_threshold_sweep_by_source.csv`,
  `data/interim/leakage_by_release_year.csv` — **new, this session**: the
  two extra leakage breakdowns requested to inform the (still-pending)
  threshold decision. Produced by a one-off analysis script (not part of
  the Phase 2 pipeline module, since `cluster_and_split.py` was run
  unchanged) reading `candidates_dedup.csv`. See "Completed phases" below
  for the numbers.
- `data/interim/pilot_500/` — the original 500-entry pilot run's Phase 1 +
  Phase 2 outputs (`candidates.csv`, `phase1_attrition.csv`,
  `candidates_dedup.csv`, `clusters.tsv`, `pesto_homology_search.tsv`,
  `leakage_threshold_sweep.csv`, `phase2_attrition.csv`, `sequences.fasta`),
  moved here (not overwritten) before the full-scale rerun, for
  before/after comparison.
- `data/raw/rcsb/{search,entries,polymer_entities}/` — cached RCSB API
  responses; now covers the full 10,774-hit search (previously just the
  500-entry pilot slice), fetched primarily via the new GraphQL batching
  path. `data/raw/rcsb/graphql_group_batches/`,
  `data/raw/rcsb/graphql_entry_batches/` — new cached raw GraphQL batch
  responses.
- `data/raw/sifts/pdb_chain_uniprot.tsv.gz` — cached bulk SIFTS file
  (~986k chain mappings loaded from it).
- `data/raw/mmseqs_work/` — MMseqs2 intermediate cluster/search DBs and tmp
  dirs from the (now full-scale) Phase 2 run (cache, not a curated output).
- `logs/select_complexes.log`, `logs/cluster_and_split.log` — full DEBUG
  logs of the full-scale Phase 1 and Phase 2 runs respectively.
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
- **RCSB metadata fetch uses batched GraphQL, not per-entry REST.** For the
  full-scale run, `data.rcsb.org/graphql` fetches an entry plus all its
  polymer entities in one call, batched `GRAPHQL_BATCH_SIZE=200` IDs per
  request (~125 entries/sec measured) instead of one REST call per entry
  plus one per entity. The GraphQL response is split into the exact same
  per-entry/per-entity cache-file shape the old REST fetchers wrote, so
  every downstream parser is unchanged. IDs GraphQL doesn't return fall
  back to the original single-item REST calls. Resumable at per-entry
  cache-file granularity regardless of batch boundaries.
- **Deposition groups (e.g. PanDDA fragment-screening campaigns) are
  collapsed to one best-resolution entry per group before the expensive
  fetch**, using RCSB's real `rcsb_entry_group_membership` field filtered
  to `aggregation_method == "matching_deposit_group_id"` specifically (not
  RCSB's separate, unrelated sequence-similarity browsing groups). This is
  a new, logged, entry-level Phase 1 attrition stage
  (`deposition_group_collapse`) — added once this authoritative field was
  found; it's complementary to, not a replacement for, Phase 2's chain-level
  MMseqs2 clustering, which still independently collapses near-duplicate
  sequences from campaigns not tagged with this field.
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

### Phase 1 — Candidate complex selection

**Pilot run (2026-09-16, `--max-entries 500`):** 1275 candidate chains
across 473 entries, 160 unique UniProt accessions. Preserved at
`data/interim/pilot_500/` — full details in git history / earlier revisions
of this file. Notable finding, carried forward below: two large
crystallographic fragment-screening campaigns (yeast Prp8–Aar2, bovine
tubulin) dominated this slice, inflating entry counts relative to true
biological diversity.

**Full-scale run (2026-09-16, no `--max-entries` cap, same day as the
pilot):** `src/data/select_complexes.py`, runnable as
`.venv/bin/python -m src.data.select_complexes`. See `PLAN.md` Phase 1 for
the GraphQL-batching and deposition-group-collapse mechanisms added for
this run.

**Wall-clock time: ~100 seconds** (measured, `time` around the full run) —
far under the 2-hour stop-and-ask threshold, so the run proceeded without
pausing (the pre-run estimate, based on the empirically measured ~125
entries/sec GraphQL fetch rate plus search-pagination overhead, was a few
minutes; the actual run beat that). Two transient connection resets during
search pagination were retried transparently by the existing Retry adapter
— no manual intervention, no data loss.

**Attrition (from `data/interim/phase1_attrition.csv`; first two stages are
entry-level, the rest are chain-level):**

| stage | n_before | n_dropped | n_remaining |
| --- | --- | --- | --- |
| entries_from_search | 10,774 | 0 | 10,774 |
| deposition_group_collapse | 10,774 | 786 | 9,988 |
| initial_protein_chains | 37,356 | 0 | 37,356 |
| non_protein_entity_type | 37,356 | 0 | 37,356 |
| min_chain_length>=40 | 37,356 | 7,551 | 29,805 |
| uniprot_mapping | 29,805 | 6,673 | 23,132 |
| chimera | 23,132 | 245 | 22,887 |

786 duplicate entries were dropped across 26 true RCSB deposit groups
(PanDDA-style campaigns); 9988/9988 remaining entries were fetched
successfully (0 fetch errors, 0 local sanity-check failures). Final:
**22,887 candidate chains across 8,808 entries, 4,029 unique UniProt
accessions** — roughly 18x the pilot's chain count and 25x its UniProt
diversity, confirming the pilot was a small, non-representative slice.

**Distributions:** resolution mean 2.00 Å (range 0.78-2.50). Release-year
distribution is now smooth (not lumpy like the alphabetically-sorted
pilot): 2018 (1907), 2019 (2754), 2020 (3012), 2021 (3078), 2022 (2667),
2023 (2736), 2024 (2440), 2025 (2599), 2026 (1694, partial year).
`sifts_agrees` True for 22,648/22,887 chains (99.0%), consistent with the
pilot's finding that disagreements are bulk-file coverage gaps, not real
conflicts.

### Phase 2 — Redundancy reduction + PeSTo-overlap flagging

**Pilot run (2026-09-16, on the 1275-chain pilot slice):** 146 clusters;
primary-mode (30% identity, train+test+validation union) leakage sweep
survivor count was 15 — well under `MIN_TEST_CHAINS_FLOOR` (50). Preserved
at `data/interim/pilot_500/`. This shortfall motivated the full-scale
Phase 1 rerun below rather than accepting the pilot's numbers as final.

**Full-scale run (2026-09-16, same day, `src/data/cluster_and_split.py`
run *unchanged* on the new 22,887-chain candidate pool):**

**Attrition (from `data/interim/phase2_attrition.csv`):**

| stage | n_before | n_dropped | n_remaining |
| --- | --- | --- | --- |
| initial_candidate_chains | 22,887 | 0 | 22,887 |
| missing_cached_sequence | 22,887 | 0 | 22,887 |
| redundancy_clustering | 22,887 | 19,878 | 3,009 |

MMseqs2 `easy-cluster` (`--min-seq-id 0.30 -c 0.80 --cov-mode 0`) collapsed
22,887 chains into **3,009 clusters**. Cluster-size distribution is heavily
right-skewed: top 10 cluster sizes are 786, 702, 627, 428, 411, 306, 270,
268, 208, 207 (these ten alone account for ~4,213 of the 22,887 candidate
chains); 762/3,009 clusters (25%) are singletons; 413 clusters have ≥10
members. This is the expected shape given how much redundancy Phase 1
already removed via deposition-group collapse — the remaining large
clusters are independent redeposits/homologs that aren't the same RCSB
deposit group.

**Leakage search result:** of the 3,009 representatives, **2,313 (76.9%)
are flagged `pesto_homolog_overlap`** (≥30% identity to the
train+test+validation union) and 834 (27.7%) are flagged
`pesto_exact_train_overlap` (literal ID match to the training file alone).

**Leakage threshold sweep (`data/interim/leakage_threshold_sweep.csv`),
against `MIN_TEST_CHAINS_TARGET=100` / `MIN_TEST_CHAINS_FLOOR=50`:**

| mode | n_survive | meets target (100) | meets floor (50) |
| --- | --- | --- | --- |
| identity < 0.30 (primary) | 696 | Yes | Yes |
| identity < 0.50 | 834 | Yes | Yes |
| identity < 0.70 | 945 | Yes | Yes |
| identity < 0.95 | 1,039 | Yes | Yes |
| exact-ID (train only) | 2,175 | Yes | Yes |

**Every mode, including the strictest (primary, 30% identity against the
full union), now clears both the target and the floor by a wide margin.**
The pilot's shortfall is fully resolved by the full-scale candidate pool —
it was a small-sample artifact, not evidence that the primary leakage
definition is unworkable.

**Extra breakdown (a) — train-only vs. the recorded train+test+validation
union, at each swept threshold** (`data/interim/leakage_threshold_sweep_by_source.csv`;
this does not change the recorded primary definition, which stays the
union):

| identity threshold | union (train+test+val, primary) | train-only |
| --- | --- | --- |
| < 0.30 | 696 | 1,105 |
| < 0.50 | 834 | 1,336 |
| < 0.70 | 945 | 1,451 |
| < 0.95 | 1,039 | 1,542 |

Including test+validation in the exclusion set (the recorded decision, per
`PLAN.md` confound (b) — both files' true roles couldn't be confidently
distinguished, see the evidence table) costs roughly 35-40% of the
survivors relative to train-only at every threshold (e.g. at 30%: 696 vs.
1105, a 409-chain/37% reduction). Both are comfortably above the floor
either way at full scale.

**Extra breakdown (b) — survivors by release year, primary (union, 30%)
definition** (`data/interim/leakage_by_release_year.csv`):

| release year | n representatives | n survive (union, primary) | % flagged (union) | n survive (train-only) | % flagged (train-only) |
| --- | --- | --- | --- | --- | --- |
| 2018 | 320 | 0 | 100.0% | 55 | 82.8% |
| 2019 | 416 | 1 | 99.8% | 67 | 83.9% |
| 2020 | 385 | 3 | 99.2% | 79 | 79.5% |
| 2021 | 328 | 96 | 70.7% | 134 | 59.1% |
| 2022 | 341 | 128 | 62.5% | 158 | 53.7% |
| 2023 | 335 | 127 | 62.1% | 166 | 50.4% |
| 2024 | 333 | 130 | 61.0% | 171 | 48.6% |
| 2025 | 331 | 122 | 63.1% | 161 | 51.4% |
| 2026 | 220 | 89 | 59.5% | 114 | 48.2% |

**This is the most surprising and most decision-relevant finding of this
session.** Chains released 2018-2020 (i.e. right after the AF2 training
cutoff) are almost entirely (99-100%) flagged as PeSTo-overlapping under
the primary definition — essentially none of them would survive the
primary leakage filter. The flag rate drops sharply and plateaus around
59-63% (union) / 48-59% (train-only) from 2021 onward. This is consistent
with PeSTo's own training-data snapshot extending to roughly 2020-2021
(PeSTo was published in 2023, so a ~2020-2021 cutoff is plausible), well
past the 2018-04-30 AF2 cutoff this project uses for its own candidate
pool — see `PLAN.md` §5 open question #2, now partially answered by this
evidence. **Practical implication:** almost all of the primary mode's 696
survivors are necessarily drawn from 2021+ releases; there are essentially
zero usable non-overlapping chains from the 2018-2020 window under the
strict/primary definition. If the eventual benchmark specifically wants
coverage close to the AF2 cutoff boundary, the primary leakage definition
as recorded will not provide it — that tension is inherent to the data
and not resolved by picking a different identity threshold within the
30-95% sweep range.

**Tradeoffs for the (still open) leakage-threshold decision:**

- Loosening the identity threshold from 30% to 95% only grows the primary
  survivor pool from 696 to 1,039 (a ~49% relative increase) — a modest
  gain for a large weakening of what "no leakage" actually means, since
  50-95%-identity homologs can still share fold/interface-pattern
  information with PeSTo's training data.
- Dropping test+validation from the exclusion set (train-only) grows the
  pool further (e.g. 696 → 1,105 at 30%) but reintroduces exactly the
  ambiguity `PLAN.md` confound (b) chose to resolve conservatively (which
  file was the true untouched holdout is unresolved from paper vs. code).
- All four swept thresholds and both source scopes (union/train-only)
  clear `MIN_TEST_CHAINS_TARGET` on their own now — sample size is no
  longer the constraint it was against the pilot slice. The remaining
  choice is a leakage-strictness/coverage tradeoff, and specifically a
  *temporal* coverage tradeoff (recent vs. right-after-cutoff releases),
  not a raw-sample-size one.

## Next step

**A human threshold decision is pending** — this session deliberately did
not choose one (see the tradeoffs above and the full sweep/breakdown CSVs
in `data/interim/`). Once decided, the choice should be recorded as the
final piece of `PLAN.md` Phase 7's leakage-sensitivity design (it's already
using `LEAKAGE_FILTER_MODES`/`SEQUENCE_IDENTITY_CUTOFF` from `src/config.py`
— no code changes needed regardless of which threshold is chosen, since
Phase 7 already re-derives eligibility from `candidates_dedup.csv`'s
per-representative identity columns rather than hardcoding 30%).

**Do not start Phase 3** until that decision is made (explicit instruction
this session). Once ready, Phase 3 is **Download PDB mmCIFs and AlphaFold
DB models** (`PLAN.md` Phase 3): `src/data/download_pdb.py`,
`src/data/download_alphafold.py`, `src/pipeline/download_structures.py`,
reading from `data/interim/candidates_dedup.csv`.

Open question #3 from `PLAN.md` §5 (the Phase 9 significance threshold)
is still unresolved and still cheap to decide now.

## How to update this file

When you finish a phase, a sub-step, or make a new decision that isn't yet
in `PLAN.md`: update "What exists on disk," move finished items out of
"Next step" and into "Decisions already made" or a new "Completed phases"
section, and bump "Last updated." Keep entries about *state and decisions*,
not a narrated log of tool calls.
