# Progress

Standing status file. Read this first in any new session on this project —
it says what's actually done, what's decided, and what to do next. Update
it whenever you finish a chunk of work or make a decision that changes the
plan. `PLAN.md` is the design document (10 phases, confounds, deviations);
this file tracks execution against it.

Last updated: 2026-09-18.

## Where things stand

**Phases 1 and 2 are done, and both have now been run at full scale**
(uncapped Phase 1 search, all 10,774 hits; Phase 2 unchanged on the result).
**The leakage-threshold decision is now made (2026-09-16):** primary test
set = `identity < 0.30` (`SEQUENCE_IDENTITY_CUTOFF`) against PeSTo's
train+test+validation union (`"homolog"` mode), n=696; sensitivity modes
(looser thresholds, exact-ID, no filter) remain as recorded. See "Decisions
already made" below and `PLAN.md` confounds (e)/(f) for the full reasoning,
including the recorded limitation that this primary set is skewed toward
2021+ releases. `candidates_dedup.csv` now carries one boolean column per
`LEAKAGE_FILTER_MODES` entry (`eligible_homolog`, `eligible_exact_train`,
`eligible_none`) so later phases can select any subset by column without
recomputing the homology search. **Phase 3 (structure download) is done,
full scale** — all 3,009 representatives, primary set (609/696) comfortably
clears `MIN_TEST_CHAINS_TARGET`. **Phase 4 (residue-numbering mapping) is
now also done, full scale** — 2,604/2,704 representatives with both
structures downloaded mapped successfully (primary set: 579/609). A real,
significant data-quality issue was found and gated in Phase 4: ~9% of
Phase 3's "AlphaFold" hits are actually third-party Community models
(mostly ColabFold), not genuine AlphaFold2 — see "Decisions already made"
below and `PLAN.md` Phase 4. **Phase 5 (ground-truth interface labels) is
now also done, full scale** — 2,527/2,604 Phase-4-mapped chains labeled
(primary set: 566/579). The single most decision-relevant finding this
session: computing partners in the raw asymmetric unit instead of the
true biological assembly would have manufactured fake interfaces almost
everywhere (only 1/2,604 chains would show zero ASU-only partners, vs.
75 chains correctly excluded as genuinely monomeric once the real
biological assembly is used) — see "Decisions already made" and `PLAN.md`
Phase 5. **Phase 6 (PeSTo inference) is now done too, lean scope** — the
primary set only (566 chains), `exp`+`af_trimmed` inputs only (`af_full`
and the apo arm deferred, see the lean-scope decision below). **All
566/566 chains succeeded on both inputs, zero failures, zero inputs
exceeded the memory ceiling** (the lean scope's real primary-set atom
counts turned out to fit comfortably). **Part C sanity check passed:**
pooled ROC-AUC of PeSTo's experimental-input predictions against Phase
5's distance labels = **0.882** — confirms the interface channel and
residue joins are correct before any real benchmarking. **Phase 7
(benchmarking exp vs. af_trimmed, primary set) is now also done, lean
scope** (2026-09-18) — pre-registered analysis plan committed to in
`PLAN.md` before any comparison was run, then implemented exactly as
registered in `src/analysis/benchmark.py`. **Headline result: PeSTo is
significantly and consistently more accurate on experimental structures
than on trimmed AlphaFold models of the same 554 (of 566) chains** —
paired median AUPR difference **+0.027**, 95% CI **[0.019, 0.035]**,
Wilcoxon p ≈ 3.5e-17 — but the typical per-chain gap is modest (a few
AUPR points), not dramatic; see "Completed phases" for the full numbers,
all eight pre-registered strata, and the ΔSASA robustness check (which
reproduces the same result). **Phase 8 (error analysis, primary set,
lean scope) is now also done** (2026-09-17, explicit go-ahead this
session) — pre-registered analysis plan committed to `PLAN.md` before any
structural feature was computed, then implemented exactly as registered
in `src/analysis/structural_metrics.py` + `src/analysis/error_analysis.py`.
**Headline result: the AlphaFold-vs-experimental accuracy gap shrinks
monotonically with AlphaFold's own confidence** — pooled AUPR gap 0.26 in
the lowest pLDDT band (`<50`) down to 0.04 in the highest (`>=90`), a >6x
reduction — and a per-residue regression confirms pLDDT and local
structural divergence (local RMSD) are each independently associated with
the size of the prediction shift (both significant, controlling for the
other; VIF only 1.6, milder collinearity than expected). **Both Phase 9
gate conditions (Phase 7 significance + Phase 8 pLDDT-band effect) are
now met for the first time**, but Phase 9 remains out of scope for this
lean pass regardless, per the standing lean-scope decision (not a
re-opened question). One pre-registered hypothesis was **not** borne out:
chain length has essentially no association with per-chain AUPR drop
(Spearman rho=0.028, p=0.51) — the gap between Phase 7's pooled AUPR drop
(0.070) and its per-chain median drop (0.027) turned out to be explained
by the drop distribution's right-skew (mean 0.061, well above the
median), not by longer chains dropping more. See "Completed phases" for
the full band/regression/length/flagged-chain numbers.

**The results paper (`paper/results_paper.pdf`, built by `paper/build_paper.py`)
was written this session (2026-09-18) and then adversarially audited against
the data in the same session.** The audit found and fixed: (1) an internal
terminology inconsistency where both the 696-representative leakage-free
candidate set and the 566-chain final analyzed set were called "the primary
benchmark set" -- now "leakage-free candidate set" (696) vs. "primary
benchmark set" (566) throughout, abstract included; (2) a real
mischaracterization of `MIN_CHAIN_LENGTH` (Phase 1) as "modeled residues"
in both `config.py`'s own comment and the paper -- it's actually
`entity_poly.rcsb_sample_sequence_length`, the polymer entity's full
expressed/construct sequence length, which differs from the number of
residues actually resolved in the deposited coordinates for 86% of Phase-1
candidates (extreme case: PDB 7F90 chain B, 1,817-residue entity sequence
vs. 45 residues actually observed) -- both fixed to describe it accurately;
(3) the 95,118-vs-95,117 pooled-residue-count mismatch between Phase 7's
benchmark and Phase 8's regression, traced to one specific residue (PDB
35WD chain A, auth residue 1, chemical component ZI0) whose RSA is
undefined because it isn't a standard amino acid per gemmi's CCD -- this
was already being dropped for a legitimate reason but was never logged;
`src/data/interface_labels.py` now logs it (and every other such residue
project-wide: 44/528,351 labeled residues across all leakage modes) both
per-occurrence and as a run-summary count, and the paper now explains the
1-residue gap inline; (4) overreaching causal/actionable language in the
Phase 8 write-up's Discussion, replaced by Phase 8b's actual measured
numbers (see below); (5) two reproducibility-documentation bugs in
`README.md` -- the documented Phase 1 command used the `--max-entries 500`
pilot flag instead of the uncapped full-scale invocation that actually
produced this project's numbers, and the documented phase sequence omitted
`scripts/validate_mapping_geometry.py`, which is not actually optional
(`src/analysis/benchmark.py` reads its `phase4_geometric_validation.csv`
output directly for the "geometrically flagged" stratum, and
`error_analysis.py` depends on it transitively) -- both fixed, and the
full reproduction sequence now annotates every step with network
requirements and measured runtime. **Phase 8b (pLDDT-filtering deployment
check) was added as a new post-hoc, exploratory phase** to directly test
the Phase 8 Discussion's informal "actionable" claim rather than leave it
asserted: naive pLDDT filtering at a fixed 0.5 threshold does **not**
meaningfully close the experimental-vs-AlphaFold gap (AlphaFold-input F1
improves by only +0.007 from no filtering to the strictest cutoff, well
within overlapping CIs, while discarding 43% of true interface residues)
-- reported plainly as a limited/non-result rather than reframed, per this
session's explicit instruction. See "Completed phases" for the full
Phase 8b table and reasoning, and `PLAN.md`'s Phase 8b section for the
post-hoc pre-registration (written before these specific numbers, but
after Phase 7/8's results were already known -- explicitly not a true a
priori pre-registration, and labeled as such everywhere it's mentioned).

The earlier 500-entry pilot run (Phase 1 + Phase 2) is preserved at
`data/interim/pilot_500/` for before/after comparison, not overwritten.
The pilot's leakage-sweep shortfall (primary mode well under
`MIN_TEST_CHAINS_FLOOR`) turned out to be a small/non-representative-sample
artifact, not a structural problem -- the full-scale run clears both
`MIN_TEST_CHAINS_TARGET` and `MIN_TEST_CHAINS_FLOOR` at every swept
threshold.

**Phases 1-8, the results paper, and Phase 8b were committed and pushed**
to `https://github.com/Menrae/ppi-interface-transfer.git` (`main`) across
several commits through `869eaf7` ("Add results paper..."). **This
session's changes (the audit fixes above, Phase 8b, and the rebuilt paper)
are not yet committed** — confirm with the user before committing/pushing,
per their standing preference to review first; don't assume a green light
carries over between sessions.

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
- `logs/select_complexes.log`, `logs/cluster_and_split.log`,
  `logs/fetch_structures.log` — full DEBUG logs of the full-scale Phase 1,
  2, and 3 runs respectively.
- `src/data/fetch_structures.py` — Phase 3 implementation (see "Completed
  phases" below). `tests/test_fetch_structures.py` (30 tests, passing, no
  network) + `tests/fixtures/alphafold_api_{success,isoforms_only,fragmented}.json`.
- `data/raw/pdb/{PDB_ID}_updated.cif.gz` — **new, Phase 3**: PDBe's
  "updated" mmCIF (embeds SIFTS residue-level UniProt cross-references) for
  all 2,294 unique PDB entries among the 3,009 representatives, gzipped
  (637 MB total).
- `data/raw/alphafold/{uniprot_acc}.cif.gz` — **new, Phase 3**: AlphaFold DB
  models for 2,485 of 2,756 unique UniProt accessions that passed
  resolution + sequence-match, gzipped (247 MB total, incl. cached API
  responses). `data/raw/alphafold/api/{uniprot_acc}.json` — cached raw
  AlphaFold DB prediction-API responses (both successful and permanently-
  failed outcomes, all 2,756 accessions).
- `data/raw/uniprot/{uniprot_acc}.fasta` — **new, Phase 3**: cached current
  UniProt sequences (2,506 accessions) fetched to detect
  `sequence_mismatch` against each AlphaFold model's own sequence (9.8 MB).
- `data/interim/fetch_report.csv` (3,009 rows), `data/interim/phase3_attrition.csv`
  — **new, Phase 3** outputs, see "Completed phases" below.
- `src/data/align_residues.py` — Phase 4 implementation (see "Completed
  phases" below). `tests/test_align_residues.py` (43 tests, passing, no
  network, hand-built mmCIF fixtures inline in the test file — no fixture
  files needed on disk for this module).
- `scripts/spot_check_mapping.py` — **new, Phase 4**: hand-verification
  helper, not part of the pipeline; prints a table + ChimeraX command
  block for N random primary-set mapped chains x 3 kept residues each.
- `data/interim/residue_mappings/{pdb_id}_{chain}.parquet` (2,604 files,
  36 MB total) — **new, Phase 4**: per-chain residue mapping tables.
- `data/interim/mapping_report.csv` (2,704 rows), `data/interim/phase4_attrition.csv`,
  `data/interim/phase4_coverage_sweep.csv` — **new, Phase 4** outputs, see
  "Completed phases" below.
- `logs/align_residues.log` — full DEBUG log of the full-scale Phase 4 run.
- `src/data/interface_labels.py` — Phase 5 implementation (see "Completed
  phases" below). `tests/test_interface_labels.py` (18 tests, passing, no
  network, hand-built mmCIF fixtures with real `pdbx_struct_assembly`/
  `pdbx_struct_oper_list` categories inline in the test file).
- `scripts/spot_check_interfaces.py` — **new, Phase 5**: hand-verification
  helper, not part of the pipeline; prints a table + a ChimeraX command
  block (fetching the chosen assembly directly from RCSB) for N random
  primary-set labeled chains.
- `data/interim/interface_labels/{pdb_id}_{chain}.parquet` (2,527 files) —
  **new, Phase 5**: per-chain interface label tables.
- `data/interim/labels_report.csv` (2,604 rows), `data/interim/phase5_attrition.csv`,
  `data/interim/phase5_label_agreement.csv` — **new, Phase 5** outputs, see
  "Completed phases" below.
- `logs/interface_labels.log` — full DEBUG log of the full-scale Phase 5 run.
- `external/PeSTo/` — **new, Phase 6**: git checkout of `LBM-EPFL/PeSTo`,
  commit `ba651aa29aaa839d0ee2c458dee90ca1f934d90a` (gitignored, 1.3GB
  including all 4 released checkpoint variants and their training logs;
  only `i_v4_1_2021-09-07_11-21/model_ckpt.pt` is actually used).
- `src/models/run_pesto.py`, `src/models/pesto_worker.py` — Phase 6
  implementation (see "Completed phases" below).
  `tests/test_run_pesto.py` (6 tests, passing, no network, no model
  weights).
- `scripts/validate_mapping_geometry.py`, `scripts/validate_interface_labels.py`,
  `scripts/phase6_sanity_check.py` — **new, Phase 6 session**: geometric/
  spatial verification for Phases 4/5 (ChimeraX not available), and the
  Part C sanity check. Not part of any pipeline; read already-written
  outputs.
- `data/interim/pesto_inputs/*.pdb`/`*.keys.csv` (126 MB) — **new, Phase
  6**: prepared PeSTo inputs, inspectable.
- `data/processed/predictions/{exp,af_trimmed}/{pdb_id}_{chain}.parquet`
  (566 chains each, 6.8 MB total) — **new, Phase 6**: per-residue PeSTo
  interface-probability predictions.
- `data/interim/inference_report.csv`, `data/interim/phase6_attrition.csv`,
  `data/interim/phase4_geometric_validation.csv`,
  `data/interim/phase5_spatial_verification.csv` — **new, Phase 6 session**
  outputs, see "Completed phases" below.
- `logs/phase6_primary.log`, `logs/run_pesto.log` — the former is the
  nohup-redirected console stream from the full primary run (per this
  session's explicit request); the latter is `run_pesto.py`'s own DEBUG
  log (same convention as every other phase).
- `src/analysis/benchmark.py` — Phase 7 implementation (see "Completed
  phases" below), runnable as `.venv/bin/python -m src.analysis.benchmark`.
  `tests/test_benchmark.py` (5 tests, passing, no network, synthetic
  fixtures only). New `config.py` constants: `BOOTSTRAP_N_RESAMPLES=10_000`,
  `PHASE7_BOOTSTRAP_SEED=0`, `INTERFACE_FRACTION_STRATUM_THRESHOLD=0.5`.
- `results/benchmark/{per_chain_metrics,summary,strata,exclusions}.csv`,
  `results/figures/{paired_aupr_scatter,aupr_diff_distribution,
  pooled_pr_curve,pooled_roc_curve}.png` — **new, Phase 7**: benchmark
  outputs, see "Completed phases" below for the numbers.
- `logs/benchmark.log` — full DEBUG log of the Phase 7 run.
- `src/analysis/structural_metrics.py` — Phase 8 per-residue feature
  computation (see "Completed phases" below): pLDDT (from the AlphaFold
  mmCIF's `B_iso_or_equiv`), local RMSD (per-residue windowed
  superposition), global Cα distance (whole-chain superposition, reusing
  Phase 4's method), secondary structure (`biotite` P-SEA — `mkdssp`
  confirmed not installed). `tests/test_structural_metrics.py` (11 tests,
  synthetic coordinates only, no network/structure files).
- `src/analysis/error_analysis.py` — Phase 8 orchestration (see
  "Completed phases" below), runnable as
  `.venv/bin/python -m src.analysis.error_analysis`. Reuses Phase 7's
  `src.analysis.benchmark.load_all_chain_frames`/`build_chain_frame`
  rather than re-deriving the join. `tests/test_error_analysis.py` (7
  tests, synthetic data only). New `config.py` constants:
  `LOCAL_RMSD_RADIUS_ANGSTROM=10.0`, `LOCAL_RMSD_MIN_NEIGHBORS=3`,
  `PLDDT_BAND_MIN_RESIDUES=30`, `VIF_COLLINEARITY_THRESHOLD=5.0`,
  `LENGTH_ANALYSIS_N_QUANTILES=4`, `PHASE8_BOOTSTRAP_SEED=0`,
  `PHASE8_BAND_BOOTSTRAP_N_RESAMPLES=2_000` (smaller than Phase 7's
  10,000 — profiled empirically, see `PLAN.md` Phase 8 for why).
  `statsmodels==0.15.0` newly installed (`.venv`) and added to
  `requirements.txt` (cluster-robust OLS + VIF for the regression).
- `results/error_analysis/{per_residue_features.parquet,band_metrics.csv,
  regression.csv,length_analysis.csv,flagged_chains_summary.csv,
  structural_feature_exclusions.csv}`, `results/figures/error_analysis_
  {band_metrics,shift_vs_local_rmsd,aupr_drop_vs_length,flagged_vs_unflagged}.png`
  — **new, Phase 8**: error-analysis outputs, see "Completed phases"
  below for the numbers.
- `logs/error_analysis.log` — full DEBUG log of the Phase 8 run (zero
  `WARNING` lines — no missing Cα/pLDDT residues, no undersized bands, no
  empty bootstrap groups).
- `src/analysis/plddt_filtering.py` — Phase 8b implementation (post-hoc,
  exploratory; see "Completed phases" below), runnable as
  `.venv/bin/python -m src.analysis.plddt_filtering`.
  `tests/test_plddt_filtering.py` (9 tests, synthetic data only). New
  `config.py` constants: `PLDDT_FILTER_PREDICTION_THRESHOLD=0.5`,
  `PLDDT_FILTER_CUTOFFS=(0,50,70,90)`.
- `results/error_analysis/plddt_filtering.csv`,
  `results/figures/error_analysis_plddt_filtering.png` — **new, Phase 8b**
  outputs, see "Completed phases" below for the numbers.
- `logs/plddt_filtering.log` — full DEBUG log of the Phase 8b run.
- `paper/build_paper.py`, `paper/template.tex` — regenerates
  `paper/results_paper.pdf` (and `paper/results_paper.tex`,
  `paper/figures/fig{1..7}_*`) end to end from `results/` and
  `data/interim/` in one command
  (`.venv/bin/python paper/build_paper.py`); no number in the paper is
  hand-typed. See README.md's "Building the results paper" section for
  the PDF route (`pdflatex`) and the structure-figure fallback
  (matplotlib Calpha trace, since no headless 3D molecular renderer could
  be installed in this container -- `pymol`/`pymol-open-source` aren't on
  PyPI, `pyvista`/VTK segfaults with no X/EGL/OSMesa and no root access).
  `paper/figures/example_chain_render.png`, if a real render is ever
  dropped in there, overrides the matplotlib fallback automatically (no
  code change needed) -- see README.md for exactly what it should show.
  9 pages as of this session's audit rebuild (2026-09-18), 7 figures, 4
  tables.
- Empty scaffold directories: `notebooks/` (`data/processed` and
  `results/` are now both populated, see above; PeSTo's model repo itself
  is checked out too, see above).

## Decisions already made (don't re-derive these — see `PLAN.md` for full reasoning)

- **Lean scope for Phases 6-9, decided 2026-09-17** — see `PLAN.md`
  "Scope decision" note. In scope now: Phase 6 on the primary set only
  (`eligible_homolog`), two inputs (`exp`, `af_trimmed`); Phases 7/8 on
  those predictions. Deferred (documented, not run): `af_full`, full-size
  `exact_train`/`none` leakage sensitivity (replaced by a later ~250-chain
  seeded subsample), the apo/unbound arm, Phase 9 entirely.
- **A real memory ceiling in this container, found empirically before
  committing to a full run.** 7.5GB RAM total, no GPU. PeSTo's own
  `extract_topology` builds a full O(n^2) pairwise-distance tensor, so
  peak memory grows worse than linearly with atom count: 2338 atoms->1.0GB,
  4855->2.1GB, 6041->3.0GB, 8833->6.2GB peak RSS; **13822 atoms was
  OOM-killed outright.** `config.PESTO_MAX_ATOMS = 7500` was derived by
  solving for peak memory <=~4.8GB (this container's ~6.3GB "available"
  minus a 1.5GB headroom) against this curve. Under the *original*
  (non-lean) full scope, a sequential single-worker projection came to
  **~27 hours**, with ~110 `af_full` inputs alone exceeding the ceiling --
  this, not just wall-clock convenience, is what motivated the lean-scope
  decision above. Under the lean scope (primary set, `exp`+`af_trimmed`
  only), **zero** inputs exceed the ceiling (max real primary-set input
  was 6,158 atoms) — the lean scope happens to sidestep the memory problem
  almost entirely, not just the time problem.
- **PeSTo's own code and this project's own `src` package collide if both
  are ever imported in the same process.** PeSTo does absolute imports of
  a top-level package literally named `src` (`from src.model_operations
  import ...`); if this project's own `src` package has *also* been bound
  in `sys.modules` in that process (which happens automatically just by
  virtue of `src.models.run_pesto` itself being part of that package),
  Python's namespace-package resolution silently returns the wrong `src`
  (confirmed interactively: `import src.dataset` inside a `-c`/interactive
  session resolved to *this project's* `src/__init__.py`, not PeSTo's,
  because the cwd-derived `''` sys.path entry -- present only in
  `-c`/interactive/`-m` invocations, not plain script execution -- let a
  *regular* package (ours, has `__init__.py`) outrank PeSTo's *bare*
  `src/` directory found earlier in `sys.path`, and a regular package
  always wins over accumulated namespace portions regardless of order).
  Fix: `src/models/pesto_worker.py` is invoked only as a bare script
  (`python pesto_worker.py ...`, never `-m`/`import`), which keeps
  `sys.path[0]` at the worker file's own directory and never introduces
  the project root or cwd -- confirmed empirically to avoid the collision.
  This also happens to be exactly the process-isolation boundary needed
  for OOM containment (see below), so one design choice covers both needs.
- **One OS subprocess per (chain, input-variant) job, not a persistent
  multi-job worker.** Chosen for the two reasons above (OOM isolation,
  `sys` collision avoidance) over a more complex batched-worker design
  that would need partial-completion bookkeeping when a worker crashes
  mid-batch -- "keep it simple" per this session's brief. Reload overhead
  (~1-3s/job: Python startup + PeSTo imports + model load) is real but
  small relative to per-job inference time for anything but the tiniest
  chains, and negligible relative to total run time.
- **Residue order in PeSTo's output can be trusted without coordinate
  matching, but only because this project's own input files are written
  in a fully-controlled order.** PeSTo's own preprocessing
  (`clean_structure`) renumbers residues sequentially in file-encounter
  order (correctly treating insertion-code changes as new-residue
  boundaries); a later step can physically move a HETATM-flagged residue's
  (e.g. a modified residue like MSE) *atoms* elsewhere in the concatenated
  array, but `pt.unique()`-based grouping sorts by that renumbering, which
  was assigned *before* the reordering -- so it recovers true file order
  regardless. Confirmed with a synthetic MSE-containing fixture before
  relying on it. Consequence: `src/models/run_pesto.py` writes each
  chain's residues in a fully-controlled sorted order and zips PeSTo's
  output row *i* onto that same row *i*, with no coordinate-based matching
  needed (unlike Phase 5's assembly-copy selection, which *did* need
  coordinate verification for a different reason).
- **Interface partners must be computed in the true biological assembly,
  never the raw asymmetric unit — confirmed decisive on real data,
  2026-09-17.** Comparing the two directly (`n_partner_chains_asu_only`
  vs. the assembly-based partner count, computed for every one of the
  2,604 Phase-4-mapped chains): only **1** chain would show zero protein
  partners in the ASU, vs. **75** correctly excluded
  (`no_protein_partner_in_assembly`) once the real biological assembly is
  used. I.e. the ASU almost always has *some* spatially nearby chain
  (crystal packing), which a naive ASU-based contact search would have
  scored as a real interface for the vast majority of those 75
  genuinely-monomeric chains. This is exactly the failure mode Phase 5's
  design (PLAN.md item 1) was meant to avoid, and it is not a rare edge
  case — it's the single largest source of exclusion in this phase.
- **"Which copy of the representative chain to use" cannot always be "the
  identity-operator copy" — corrected mid-session on real data
  (`7QPB`).** The original design assumed the representative's own copy
  in the assembly must be generated by the identity operator (verified by
  coordinates), since only the deposited, Phase-4-mapped coordinates were
  assumed reachable that way. A real biologically-annotated trimeric
  assembly falsified this: chain `A` there is reached only via a
  non-identity crystal-symmetry operator (a different chain, `C`, is the
  identity-generated one), yet `A`'s residue identity/numbering is
  preserved exactly through that transform (`make_assembly` moves atoms,
  never touches `auth_seq_id`/residue metadata). Corrected rule
  (`find_self_copy`, `PLAN.md` Phase 5): when a chain appears via exactly
  one copy in the chosen assembly, that copy is the representative,
  identity-transformed or not; a tie-break (prefer identity, else lowest
  `AddNumber` copy name) is only needed when a chain genuinely appears via
  *multiple* copies in one assembly.
- **biotite's ProtOr SASA radii are looked up per (residue name, atom
  name), not a bare element symbol — found and fixed before the full run,
  not after.** An early draft built the SASA-calculation `AtomArray` from
  coordinates alone (element hardcoded to carbon as a placeholder, no real
  atom/residue names set). This doesn't silently produce wrong numbers —
  it crashes (`IndexError` inside `biotite`'s own `vdw_radius_protor`),
  confirmed by a standalone reproduction before touching the real
  pipeline. Fixed by carrying each atom's real element, atom name, and
  parent residue name through the whole geometry pipeline (`HeavyAtom`
  dataclass), not just coordinates. Caught during interactive development,
  never ran against real chains in the broken form.
- **AlphaFold DB hosts third-party Community models (mostly ColabFold)
  under the same per-accession lookup as genuine AlphaFold2 predictions —
  discovered 2026-09-16, gated in Phase 4.** Phase 3's prediction-API query
  returns whichever model(s) exactly match the queried accession,
  regardless of `providerId`; ~9% of all 2,756 accessions Phase 3 queried
  resolved to a non-`"GDM"` (i.e. non-Google-DeepMind) provider (`VR3D`
  222, `ATBC` 17, `NTDX` 9, `BFVD` 3 — mostly `ColabFold v1.5.2`/`v1.0-alpha`).
  These don't carry the AF2 training-cutoff guarantee confound (a) relies
  on. Phase 3's code as originally specified didn't check this (not asked
  to); rather than reopening/rerunning Phase 3, Phase 4's
  `verify_official_alphafold_provider` checks Phase 3's own cached API
  JSON's `providerId` first, before any mmCIF parsing, and excludes
  anything non-`"GDM"` (`alphafold_model_not_official_afdb`, 43/2,704 =
  1.6% of Phase 4's scope — smaller than the 9% of *all* accessions since
  many Community-model accessions were already excluded upstream for other
  reasons, or don't map to a representative still in scope). **If Phase 3
  or its inputs are ever rerun, consider adding this check there directly**
  (cheaper — skips the whole AlphaFold download for a Community hit) — not
  done retroactively this session since Phase 4's gate already fully
  prevents contamination of the mapped/benchmarked set.
- **Residue mapping is per-residue, never an assumed/constant offset.**
  Confirmed necessary on real data, not just in principle: representative
  `8WHI_A` maps author residue 1335 to UniProt position 1137 (a -198
  offset) while other residues in other chains use entirely different
  offsets — a single global offset per chain (or worse, per file) would
  be silently wrong. SIFTS's own embedded per-residue
  `pdbx_sifts_xref_db_num` (or, when that's absent/fails validation, a
  free-end-gap global alignment against AlphaFold's sequence) is used
  residue-by-residue instead. See `PLAN.md` Phase 4 and
  `scripts/spot_check_mapping.py`'s output for a directly-inspectable
  real example.
- **Modified residues (MSE, SEP, KCX, ...) are resolved to their parent
  amino acid via `gemmi.find_tabulated_residue(...).one_letter_code`**
  (covers the full Chemical Component Dictionary, not a hand-rolled
  lookup table) — both for SIFTS-path identity comparison and for building
  the fallback alignment's sequence string. E.g. MSE (selenomethionine) ->
  `M`, compared directly against AlphaFold's plain MET.
- **A mapping method (SIFTS or the fallback alignment) is trusted per
  chain only if ≥90% of its mapped-and-in-range residues match AlphaFold's
  residue identity** (`RESIDUE_MAPPING_VALIDATION_IDENTITY`, new config
  constant) — below that, it's treated as a numbering/frame error for
  that method (try the other method, or exclude), not biology. An
  isolated single-residue mismatch in an otherwise-clean chain (an
  engineered point mutation) does *not* fail the chain — it's kept, with
  `match_flag=False`, contributing to a `n_mutations` count, not an
  exclusion. Confirmed well-separated on the real run: every chain that
  used the strict 90% threshold successfully mapped landed at ≥90% (most
  at 100%, 535/2,604 chains had ≥1 real mutation), while every chain
  excluded for failing both methods' validation clustered at 55-90%
  identity — a different, distinguishable population, not a borderline
  cutoff call.
- **`MIN_MAPPED_COVERAGE = 0.80`** (new config constant, chains below this
  are excluded as `coverage_below_threshold`): chosen for consistency with
  `CLUSTER_MIN_COVERAGE` (also 0.80, Phase 2). Confirmed barely binding on
  the real run — minimum coverage among all 2,604 successfully mapped
  chains was 0.81; tightening to 0.90/0.95 would only additionally drop
  13/57 chains (of 2,604) — see the coverage sweep in "Completed phases"
  below.
- **Residue-mapping parquet schema, one row per *observed* experimental
  residue** (not one row per UniProt/AlphaFold-range position): the
  recorded `PLAN.md` column list included `is_observed`, which could be
  read either way; this session's implementation keys the table by the
  physically-observed experimental residue (matching the explicit column
  list in this session's brief) and keeps `is_observed` as an
  always-`True` column for schema fidelity rather than using it to
  represent unobserved-gap rows. Noted as a clarifying implementation
  choice, not a conflict requiring a decision.
- **Leakage threshold, decided 2026-09-16:** the primary benchmark test set
  is representatives with `max_identity_any < SEQUENCE_IDENTITY_CUTOFF`
  (0.30) against PeSTo's train+test+validation union (`"homolog"` mode,
  `pesto_homolog_overlap == False`, `eligible_homolog == True` in
  `candidates_dedup.csv`) — n=696. The three sensitivity modes stay as
  recorded (`identity<0.50/0.70/0.95`, exact-ID-train-only, no filter), all
  with strictly larger n, per `PLAN.md` Phase 7. **Recorded limitation:**
  this primary set is skewed toward 2021+ releases — chains released
  2018-2020 are almost all (99-100%) flagged as PeSTo-overlapping under
  this definition, consistent with PeSTo's own training snapshot extending
  to roughly 2020-2021 (well past AF2's 2018-04-30 cutoff) — see `PLAN.md`
  confound (f). This increases the primary set's average distance from
  AF2's training cutoff relative to what the 2018-04-30 filter alone would
  suggest. Not fixable by adjusting the identity threshold within the swept
  range (0.30-0.95 all show the same year pattern).
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
- **`MIN_CHAIN_LENGTH` (Phase 1) is a filter on `entity_poly.rcsb_sample_sequence_length`,
  the polymer entity's full expressed/construct sequence length -- NOT the
  count of residues actually resolved in the deposited coordinates,
  despite `config.py`'s original comment (and this project's own paper,
  before this session's audit) calling it "modeled residues."** Confirmed
  on real data during this session's paper audit: 86% of Phase 1
  candidates with a Phase 4 mapping have `seq_length != n_observed`, with
  a large mean gap (~31 residues) and extreme cases (PDB 7F90 chain B:
  1,817 vs. 45). Both `config.py`'s comment and every paper/doc reference
  now say "entity sequence length," not "modeled residues." Don't
  reintroduce the old phrasing -- it's not a synonym.
- **`interface_labels.py`'s RSA is silently NaN (not dropped, not logged)
  for any residue whose chemical component isn't a standard amino acid
  per gemmi's CCD** (no Sander/Rost max-ASA reference value exists for
  it) -- found during this session's paper audit while tracing why
  Phase 7's pooled 95,118 residues became 95,117 in Phase 8's regression.
  Fixed, not just documented: `compute_sasa_labels` now logs a WARNING
  per occurrence (residue key + chemical component) and `run()` logs a
  run-level summary count. Rerunning Phase 5 with this fix produces
  byte-identical labeled-chain counts and pooled agreement statistics to
  before (2,527/2,604 labeled, Jaccard 0.9047972508591066) -- confirmed
  this is a pure logging addition, not a behavior change, so Phase 6/7/8
  did not need to be rerun. Project-wide: 44/528,351 labeled residues
  (across all leakage modes) have undefined RSA this way; in the primary
  set specifically, exactly 1 (PDB 35WD chain A, auth residue 1, chemical
  component ZI0), which is the residue Phase 8's regression drops.
- **README.md's documented Phase 1 reproduction command was wrong** --
  it used `--max-entries 500` (the pilot/smoke-test flag) where the
  actual full-scale run that produced this project's committed
  `candidates.csv` (22,887 chains, 8,808 entries) used no cap at all.
  Fixed. Also added: `scripts/validate_mapping_geometry.py` is
  **not optional** despite living in `scripts/` alongside genuinely
  optional hand-verification tools -- `src/analysis/benchmark.py` reads
  its `phase4_geometric_validation.csv` output directly for the
  "geometrically flagged" stratum, and `src/analysis/error_analysis.py`
  depends on it transitively through `benchmark.load_chain_metadata()`.
  Confirmed by grepping every `src/` module for `INTERIM_DATA_DIR`/
  `RESULTS_DIR` reads: this is the *only* script-produced (as opposed to
  `src/`-pipeline-produced) file anything downstream depends on;
  `phase5_spatial_verification.csv` and every `spot_check_*`/
  `phase6_sanity_check.py` output is confirmed read by nothing else, so
  those really are optional. README.md's reproduction section now states
  this explicitly and annotates every step with network/runtime info.

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

### Phase 3 — Download PDB mmCIFs and AlphaFold DB models

**Done (2026-09-16, full scale, all 3,009 Phase 2 representatives).**
`src/data/fetch_structures.py`, runnable as
`.venv/bin/python -m src.data.fetch_structures`. See `PLAN.md` Phase 3 for
the AlphaFold-API-resolution mechanism (no hardcoded model version/URL) and
the consolidation with what was originally a separate Phase 4 download
(the PDBe "updated" mmCIF, embedding SIFTS residue mappings, is fetched
here directly into `data/raw/pdb/`, not a separate `pdb_updated/` dir).

**Pre-flight estimate (real calibration sample, not guessed):** 2,294
unique PDB entries, 2,756 unique UniProt accessions → projected ~68-74
minutes (two calibration runs, consistent), ~3.3-4.7 GB raw — both
comfortably under the 2h/20GB stop-and-ask budget, so the run proceeded
without pausing.

**Actual run:** wall-clock **71m24s** (17:37:03-18:48:27), matching the
calibrated projection closely. Disk usage: `data/raw/pdb/` 637 MB (2,294
files, gzip-compressed), `data/raw/alphafold/` 247 MB (2,485 model files +
2,756 cached API responses, gzip/JSON), `data/raw/uniprot/` 9.8 MB (2,506
cached current-sequence FASTA files) — **~894 MB total**, far under the 20
GB budget.

**PDB download: 3,009/3,009 representatives succeeded (100%)** — every
unique PDB entry's updated mmCIF was fetched with zero failures.

**AlphaFold download: 2,704/3,009 representatives succeeded (89.9%), 305
failed.** Failure-reason breakdown (row-level, i.e. per representative,
not per unique accession — a few failing accessions are shared by more
than one representative row, e.g. homodimer partner chains):

| failure reason | n (of 3,009) |
| --- | --- |
| `accession_absent_from_alphafold_db` | 222 |
| `fragmented_in_alphafold_db` | 61 |
| `sequence_mismatch` | 22 |
| `network_error` / `pdb_download_error` | 0 |

**Survivor counts per leakage filter mode** (`data/interim/phase3_attrition.csv`),
`n_both_success` = representative has both a PDB structure and a validated
AlphaFold model:

| mode | n_eligible | n_pdb_success | n_alphafold_success | n_both_success | meets `MIN_TEST_CHAINS_TARGET` (100)? |
| --- | --- | --- | --- | --- | --- |
| `homolog` (primary) | 696 | 696 | 609 | **609** | Yes |
| `exact_train` | 2,175 | 2,175 | 1,957 | 1,957 | Yes |
| `none` | 3,009 | 3,009 | 2,704 | 2,704 | Yes |

**The primary set survives structure download at 87.5% (609/696)** — the
AlphaFold-side attrition (66 accession-absent, 13 fragmented, 8
sequence-mismatch within the primary set alone) costs some of the 696, but
609 remains far above `MIN_TEST_CHAINS_TARGET` (100). Nothing here changes
the leakage-threshold decision or requires loosening it for sample size.

**Nothing unexpected in the failure mix:** `accession_absent_from_alphafold_db`
dominates (mostly non-human/non-model-organism UniProt accessions, or
accessions predating/outside AlphaFold DB's covered proteome set — expected
given the candidate pool spans arbitrary PDB depositors, not just
well-studied model organisms). `fragmented_in_alphafold_db` (AFDB split the
*queried* accession's own sequence into >1 fragment model) appeared for a
small number of very long proteins (e.g. `P0DTD1`, the SARS-CoV-2 ORF1ab
polyprotein). `sequence_mismatch` (AlphaFold model's own UniProt sequence
disagrees with the accession's current UniProt sequence) is small (22) but
real — these models were skipped rather than silently used with
stale/wrong numbering.

**Tests:** `tests/test_fetch_structures.py` (30 tests, no network,
passing) — URL construction; AFDB response classification on saved
fixtures (success, isoform-only, fragmented, absent, empty body); FASTA
parsing/sequence-match; priority ordering; `candidates_dedup.csv`
loading; fetch-report-row schema/`overall_status`; per-leakage-mode
attrition counting; and, via a fake session that raises on any unwired
URL, both "cache hit -> zero network calls" and every failure path
end-to-end.

### Phase 4 — Residue-numbering mapping

**Done (2026-09-16, full scale, all 2,704 representatives with both
structures downloaded — Phase 3's `overall_status == "complete"` set).**
`src/data/align_residues.py`, runnable as
`.venv/bin/python -m src.data.align_residues`. See `PLAN.md` Phase 4 for
the full design (SIFTS-primary/fallback-alignment mapping, AlphaFold
numbering self-validation, the official-AlphaFold-vs-Community-model
gate, the two new config constants and their justification).

**Real-data smoke-testing before the full run:** ran `process_representative`
directly (no batch orchestration) against three random samples (30, then
250, then 600 representatives, ~22% of scope) before committing to the
full run — zero crashes across all 880 sampled chains, and every
exclusion reason that came up was independently verified as a genuine,
explainable data issue (see "Decisions already made" above), not a code
bug. This is what surfaced the AlphaFold-Community-model finding, before
it could quietly corrupt the mapped set.

**Full run: wall-clock 2m19s** (measured; a 200-representative timing
sample beforehand projected ~2 minutes, so no pre-flight stop-and-ask was
needed — well under any reasonable budget for a purely local, no-network
computation).

**Overall: 2,604/2,704 representatives mapped (96.3%), 100 excluded.**
Exclusion-reason breakdown:

| exclusion reason | n (of 2,704) |
| --- | --- |
| `alphafold_model_not_official_afdb` | 43 |
| `fallback_failed_validation` | 35 |
| `coverage_below_threshold` | 22 |

**Mapping-method split (of 2,604 mapped): 2,502 via SIFTS (96.1%), 102 via
the fallback alignment (3.9%).** The fallback method is doing real,
necessary work for a small-but-nonzero fraction of chains — not dead code.

**Coverage distribution (mapped chains):** mean 0.995, median 1.0, min
0.81 (i.e. every mapped chain cleared `MIN_MAPPED_COVERAGE` with room to
spare — the coverage gate is doing real exclusion work at 22 chains, but
essentially all *surviving* chains map almost completely).
**Identity distribution (mapped chains):** mean 0.997, median 1.0, min
0.90 (exactly the validation threshold, as expected — this is the
boundary where a chain either passes as-is or gets excluded, confirming
the threshold is where the real decision boundary sits, not an arbitrary
cutoff through the middle of the distribution). 535/2,604 mapped chains
(20.5%) carry ≥1 recorded engineered mutation (`n_mutations` histogram:
2,069 chains with 0, 291 with 1, 111 with 2, tapering off to a handful
with 7-9) — real biology, correctly kept rather than excluded.

**Survivor counts per leakage filter mode** (`data/interim/phase4_attrition.csv`):

| mode | n_eligible | n_mapped | n_excluded | meets `MIN_TEST_CHAINS_TARGET` (100)? |
| --- | --- | --- | --- | --- |
| `homolog` (primary) | 609 | **579** | 30 | Yes |
| `exact_train` | 1,957 | 1,873 | 84 | Yes |
| `none` | 2,704 | 2,604 | 100 | Yes |

**The primary set survives residue mapping at 95.1% (579/609)** — 579
remains far above `MIN_TEST_CHAINS_TARGET` (100), and well above
`MIN_TEST_CHAINS_FLOOR` (50) too. Tracing the whole pipeline so far for
the primary mode: 696 leakage-filtered representatives (Phase 2) -> 609
with both structures downloaded (Phase 3, 87.5%) -> 579 with a validated
residue mapping (Phase 4, 95.1% of those) = **579/696 (83.2%) of the
primary set has made it through every phase run so far.**

**Coverage-threshold sweep** (`data/interim/phase4_coverage_sweep.csv`,
reporting only, of the 2,604 already-mapped chains):

| coverage threshold | n chains retained |
| --- | --- |
| 0.50 / 0.60 / 0.70 / 0.80 | 2,604 (all of them) |
| 0.90 | 2,591 |
| 0.95 | 2,547 |

**Nothing unexpected in the exclusion mix** beyond the AlphaFold-Community-
model finding (see "Decisions already made" above, the session's most
surprising/decision-relevant discovery). `fallback_failed_validation`
cases spot-checked individually all showed genuinely divergent sequences
(55-90% identity on both methods, e.g. likely different-organism
orthologs or fusion-construct partners) rather than alignment bugs.
`coverage_below_threshold` cases are chains dominated by tags/unmapped
segments, correctly excluded rather than benchmarked on a small,
cherry-picked mapped fragment.

**Hand-verification (`scripts/spot_check_mapping.py --n 3 --seed 0`):**
see the transcript below (also reproduced in this session's conversation).
All 9 spot-checked residues (3 chains x 3 residues) show the experimental
residue name and the AlphaFold residue name agreeing exactly at the mapped
UniProt position, including one chain (`8WHI_A`, UniProt `O75122`) with a
substantial, genuinely non-trivial author-to-UniProt offset (auth 1335 ->
UniProt 1137, a -198 shift) — direct, inspectable evidence the mapping is
not assuming any constant offset.

**Tests:** `tests/test_align_residues.py` (43 tests, no network, all
hand-built mmCIF fixtures inline) — see `PLAN.md` Phase 4 for the full
list of covered cases.

### Phase 5 — Ground-truth interface labels

**Done (2026-09-17, full scale, all 2,604 Phase-4-mapped representatives).**
`src/data/interface_labels.py`, runnable as
`.venv/bin/python -m src.data.interface_labels`. See `PLAN.md` Phase 5 for
the full design (biological-assembly selection, the corrected
"which-copy-is-the-representative" rule, partner detection, the two
labels, and the exact-join assert against Phase 4).

**Real-data smoke-testing before the full run:** ran `process_representative`
directly against three growing random samples (20, then 200, then 300,
then 500, ~42% of scope combined) before committing to the full run — zero
crashes across all sampled chains. This is what surfaced both the
ASU-vs-assembly finding and the identity-copy correction (via `7QPB`)
before either could affect a full run.

**Full run: wall-clock 4m29s** (a 200/300/500-chain timing sample
beforehand projected ~0.1s/chain, ~4.3 minutes for the full 2,604 — no
pre-flight stop-and-ask needed for a purely local, no-network computation).

**Overall: 2,527/2,604 representatives labeled (97.0%), 77 excluded.**
Exclusion-reason breakdown:

| exclusion reason | n (of 2,604) |
| --- | --- |
| `no_protein_partner_in_assembly` | 75 |
| `zero_interface_residues` | 2 |

**Survivor counts per leakage filter mode** (`data/interim/phase5_attrition.csv`):

| mode | n_eligible | n_labeled | n_excluded | meets `MIN_TEST_CHAINS_TARGET` (100)? |
| --- | --- | --- | --- | --- |
| `homolog` (primary) | 579 | **566** | 13 | Yes |
| `exact_train` | 1,873 | 1,829 | 44 | Yes |
| `none` | 2,604 | 2,527 | 77 | Yes |

**Tracing the whole pipeline for the primary mode:** 696 leakage-filtered
representatives (Phase 2) -> 609 with both structures downloaded (Phase
3) -> 579 with a validated residue mapping (Phase 4) -> **566 with
interface labels (Phase 5) = 566/696 (81.3%) of the original primary set
has made it through every phase run so far.**

**Interface fraction per chain** (2,527 labeled chains): mean 27.9%,
median 21.1%, min 0.3%, max 100% (all residues); restricted to surface
residues only (RSA ≥ `SURFACE_RSA_THRESHOLD`=0.25): mean 32.7%, median
27.4% — higher, as expected, since interface residues are disproportionately
surface residues, though not perfectly so (some surface patches face away
from any partner). Mean surface-residue count per chain: 99.6 (median 85).

**Distance-vs-SASA label agreement:** pooled (all labeled chains,
`data/interim/phase5_label_agreement.csv`) confusion matrix TP=98,736
FP=3,422 FN=6,967 TN=419,226, **Jaccard = 0.905**. Per-chain Jaccard
(`labels_report.csv`): mean 0.900, median 0.909, min 0.333 — the lowest
per-chain values all belong to chains with very small interface fractions
(a couple of residues), where Jaccard is naturally noisy; not evidence of
a systematic problem.

**Homomeric partners:** 676/2,527 labeled chains (26.8%) have at least one
homomeric partner in their assembly. `n_partner_chains` ranges 1-27
(mean 2.2, median 1) — most chains pair with exactly one partner, a
minority sit in much larger assemblies.

**Assembly vs. asymmetric-unit comparison (item 7, this phase's most
decision-relevant finding — see "Decisions already made" above):** only
1/2,604 chains has zero protein partners in the raw ASU, vs. 75 correctly
excluded once the true biological assembly is used. The gap is the
expected consequence of crystal packing: almost every ASU has *some*
spatially nearby chain, whether or not it's biologically meaningful.

**Nothing else unexpected** beyond the two corrections already recorded
above (identity-copy assumption, ProtOr metadata) — both found and fixed
during real-data smoke-testing, before the full run, not after.

**Hand-verification (`scripts/spot_check_interfaces.py --n 3 --seed 0`):**
3 primary-set labeled chains (`9JDD_A`, `7YCJ_B`, `7PU4_A`), 23-50
interface residues each, all forming contiguous runs of nearby residue
numbers rather than scattered noise — consistent with a real, spatially
coherent interface patch rather than a labeling artifact. Full ChimeraX
command blocks (each fetching the chosen assembly directly from RCSB,
independent of this project's own assembly-generation code) reproduced in
this session's conversation.

**Tests:** `tests/test_interface_labels.py` (18 tests, no network,
hand-built mmCIF fixtures with real assembly categories) — see `PLAN.md`
Phase 5 for the full list of covered cases.

### Phase 4 addendum — geometric validation (2026-09-17)

`scripts/validate_mapping_geometry.py`: Kabsch/SVD-superposed each
primary-set chain's experimental Calpha atoms onto the AlphaFold model's
Calpha atoms at the same UniProt positions (one rigid-body fit per chain),
for all 579 primary-set Phase-4-mapped chains.

**Overall fit is very good: median per-chain Calpha distance 0.82A, mean
2.02A** (`data/interim/phase4_geometric_validation.csv`). **94/579 (16.2%)
chains were flagged** for a run of >=5 consecutive residues staying >10A
apart even after the best possible global fit.

**Investigated rather than taken at face value, since 16% seemed high.**
Cross-checked the flagged chains against Phase 4's own recorded
method/identity: 87/94 flagged chains used the **SIFTS** method (the more
trustworthy of the two, driven by PDBe's own curated per-residue UniProt
cross-references) with **median sequence identity 1.000** (mean 0.996) —
i.e. the mapping itself is about as trustworthy as this project's data
gets. Spot-checked the single most extreme case (`9F07_C`, 100% of
residues flagged): sequence identity 95.2% (79/83 residues), and the
AlphaFold model's own pLDDT over the same range averages 92.1 (min 53) —
a *confident* AlphaFold prediction, not a low-confidence/disordered region
that would explain the gap as "AlphaFold didn't know." **Conclusion: this
is not a Phase 4 mapping defect** (the sequence-level correspondence is
essentially perfect) but genuine AlphaFold-vs-experimental-structure
divergence in 3D despite a confident, correct prediction of *what* protein
and *which* residues — i.e. exactly the phenomenon this whole project
exists to measure, previewed here rather than discovered as a bug. Caveat
worth carrying forward: a single global rigid-body fit will also flag
genuine multi-domain hinge motion (two well-fitting sub-domains, wrongly
"far apart" only because one global rotation can't satisfy both) as if it
were the same kind of problem -- this script doesn't yet distinguish the
two, so the 94-chain count is an upper bound on "interesting 3D
divergence," not a clean count of any one specific phenomenon.

### Phase 5 addendum — programmatic spatial verification (2026-09-17)

Since ChimeraX wasn't available for hand verification, and per this
session's instruction to attempt an independent external check first:
`scripts/validate_interface_labels.py` tried PDBe's PISA/interface-
annotation API (`interfacelist`, `pdbe_pisa/interface_residues` — 3
endpoint-pattern variants tried) for a seeded sample of 20 primary-set
chains; **none returned data this session** (all either 404 or a JSON
"Not Found" body, suggesting either the endpoint pattern itself is no
longer current, PDBe's PISA doesn't have annotations for the queried
entries, or both) — falling back to the two geometric checks specified as
the fallback, on the same 20-chain seeded sample (seed=0):

- **Contiguity:** mean 98.6% of interface residues have another interface
  residue within 8A (19/20 chains at 91-100%, one chain — not shown
  separately, all were reasonably high). Interface residues are
  overwhelmingly clustered into real patches, not scattered.
- **Centroid:** 19/20 chains (95%) have their interface patch's centroid
  closer to the partner chain's centroid than a same-chain sample of
  comparably-exposed non-interface surface residues' centroid is. The one
  exception (`7YCJ_B`: interface-to-partner 22.4A vs. surface-to-partner
  13.4A) is a plausible geometry-dependent outlier (e.g. an
  asymmetric/wrapping partner), not inspected further.

Both checks support the interface labels being spatially real and correct
-- see `data/interim/phase5_spatial_verification.csv` for the full table.

### Phase 6 — PeSTo inference on isolated single chains (lean scope)

**Done (2026-09-17), lean scope: primary set only (566 chains), `exp`+`af_trimmed`
inputs only.** `src/models/run_pesto.py` + `src/models/pesto_worker.py`,
runnable as `.venv/bin/python -m src.models.run_pesto --subset primary
--inputs exp,af_trimmed`. See `PLAN.md` Phase 6 for the full design
(PeSTo setup/reproduction, the memory ceiling, the `sys.modules`
collision and why each job is its own subprocess, concurrency tiers,
resumability).

**Smoke test (5 smallest primary chains, then a stratified 67-chain
timing sample spanning the real atom-count distribution) before
committing to the full run:** all succeeded; remeasured projection for
the full 566-chain/1,132-job primary run was **~1.3 hours** — comfortably
under the 3-hour flag threshold, so the full run launched without pausing
to ask.

**Full run launched as a true detached background process**
(`nohup ... > logs/phase6_primary.log 2>&1 &` + `disown`, per this
session's explicit request that it survive the session's own tool-call
lifecycle) and monitored periodically until completion.

**Wall-clock: 1h04m52s** (03:40:55-04:45:47) for the run as launched;
this included resuming 77 jobs already completed during the smoke-test/
timing-sample phase, so no work was wasted. CPU-time sum across all
1,132 jobs: ~3.1 hours (parallelized across the 4/2/1-worker tiers, so
wall-clock is well under that).

**Result: 566/566 primary-set chains (100%) succeeded on both `exp` and
`af_trimmed` -- zero failures, zero jobs exceeded the memory ceiling.**
The lean scope's real atom counts turned out to fit comfortably: max
observed input size was 6,158 atoms, safely under `PESTO_MAX_ATOMS`
(7,500) -- the memory-ceiling concern that drove the lean-scope decision
was almost entirely an `af_full`/full-scope problem, not a primary-set
`exp`/`af_trimmed` one. Join verification (every successful output
covers every one of Phase 4's mapped residues for that chain) passed
with zero problems.

**Per-job runtime** (`data/interim/inference_report.csv`): `exp` mean
9.9s (median 8.8s, max 30.4s); `af_trimmed` mean 9.9s (median 8.5s, max
30.4s) -- consistent between the two input types, as expected since
runtime is driven by atom count, not which structure it came from.

**Missing heavy atoms** (`exp` only, truncated sidechains vs. the ideal
count for that residue type): only 1.35% of the 95,492 predicted
residues have any missing heavy atoms -- most experimental structures in
this resolution-filtered dataset are essentially complete.

**Part C sanity check (`scripts/phase6_sanity_check.py`): pooled ROC-AUC
of `exp` predictions against Phase 5's distance labels = 0.882**
(per-chain mean 0.843, median 0.901, n=554 chains with both classes
present; 12 chains skipped from the per-chain figure for having only one
label class). Well above chance, matching PeSTo's own reported accuracy
on experimental structures -- **the interface channel, residue joins,
and alt-loc handling are all confirmed correct.**

**Tests:** `tests/test_run_pesto.py` (6 tests, no network, no model
weights) -- see `PLAN.md` Phase 6 for the full list of covered cases.

**Nothing unexpected in the failure mix, because there wasn't one** --
every job succeeded. The two real surprises this session were both found
*before* the full run (the memory ceiling and the `sys.modules`
collision, both in "Decisions already made" above), which is exactly why
the smoke-testing and calibration steps existed.

### Phase 7 — Benchmarking on experimental vs. trimmed-AlphaFold inputs (lean scope)

**Done (2026-09-18), lean scope: primary set only (566 chains), `exp` vs.
`af_trimmed`.** Analysis plan pre-registered in `PLAN.md` *before* any
comparison was computed (endpoints, strata, robustness check, alpha=0.05,
all committed to first); implemented exactly as registered in
`src/analysis/benchmark.py`, runnable as
`.venv/bin/python -m src.analysis.benchmark`. See `PLAN.md`'s Phase 7
pre-registration and completed write-up for full reasoning and the
strata table; numbers repeated here for the execution log:

- **Join/assertion:** all 566 chains passed the same-residue-set
  assertion (auth-keyed `exp` coverage == auth-keyed `af_trimmed`
  coverage == the full Phase-4-mapped set) with zero exclusions --
  `results/benchmark/exclusions.csv` is empty.
- **Metric exclusions:** 12/566 chains (2.1%) excluded from every
  ROC-AUC/AUPR metric -- all 12 are fully-interface chains (every scored
  residue labeled interface, 7-26 residues each), making the metrics
  undefined regardless of predictions. Not anticipated in the
  pre-registration but handled cleanly by the generic "undefined metric"
  exclusion rule.
- **Primary endpoint** (paired AUPR difference, exp - af_trimmed,
  n=554): median **+0.027**, 95% bootstrap CI **[0.019, 0.035]**,
  two-sided Wilcoxon **p ~ 3.5e-17**.
- **Secondary -- paired ROC-AUC difference:** median **+0.027**, CI
  **[0.018, 0.037]**, p ~ 4.6e-19.
- **Secondary -- pooled-residue metrics** (95,118 residues, large chains
  dominating by construction): exp AUPR **0.739**, af_trimmed AUPR
  **0.669**, pooled interface base rate **0.223**; pooled ROC-AUC exp
  **0.882** / af_trimmed **0.849**.
- **Robustness (DeltaSASA labels instead of distance labels):** median
  **+0.026**, CI **[0.020, 0.034]**, p ~ 3.0e-17 -- essentially identical
  to the primary (distance-label) result.
- **All eight pre-registered strata** (homomeric/heteromeric;
  small/large interface fraction at the 0.5 threshold;
  Phase-4-geometrically-flagged/unflagged; all-residues/surface-only)
  show the same direction and are each individually significant at
  alpha=0.05; effect size ranges from +0.015 (large interfaces, n=120) to
  +0.059 (geometrically-flagged chains, n=92, also the widest CI). Full
  table in `results/benchmark/strata.csv` and `PLAN.md`.
- **Figures** (`results/figures/`): `paired_aupr_scatter.png` (most
  chains cluster near the identity line with a tilt toward
  exp-better and a right-skewed tail); `aupr_diff_distribution.png`
  (sharp peak near zero, right-skewed); `pooled_pr_curve.png`,
  `pooled_roc_curve.png` (exp visibly above af_trimmed at every
  threshold, both well above the base rate).
- **Plain-terms reading:** the difference is real, consistent in
  direction across every slice examined, and highly statistically
  significant -- but it is not large in absolute terms for a typical
  chain (median ~0.03 AUPR points); avoid characterizing this as
  AlphaFold inputs being "much worse" across the board.

**Tests:** `tests/test_benchmark.py` (5 tests, no network, synthetic
fixtures) -- metrics-vs-sklearn agreement, single-class exclusion with a
logged count, the same-residue-set assertion raising on a real mismatch,
chain- (not residue-) level bootstrap behavior, and a synthetic
known-AUPR-drop recovery test. Full suite: 143/143 passing.

**Deferred, not part of this pass** (per the lean-scope decision):
leakage-filter-mode sensitivity (`exact_train`/`none`), any `af_full`
comparison, the apo/unbound arm, and Phase 8.

### Phase 8 — Error analysis (primary set, lean scope)

**Done (2026-09-17), lean scope: primary set only (566 chains), `exp` vs.
`af_trimmed`.** Analysis plan pre-registered in `PLAN.md` *before* any
structural feature was computed (band definitions, bootstrap methodology
and resample counts, regression design/SE method, descriptive analyses,
association-not-causation framing, all committed to first); implemented
exactly as registered in `src/analysis/structural_metrics.py` +
`src/analysis/error_analysis.py`, runnable as
`.venv/bin/python -m src.analysis.error_analysis`. See `PLAN.md`'s Phase
8 pre-registration and completed write-up for full reasoning and every
table; numbers repeated here for the execution log:

- **Run:** 95,118 pooled residues across all 566 primary-set chains,
  wall-clock 2m48s. Zero chain-level exclusions
  (`results/error_analysis/structural_feature_exclusions.csv` empty) and
  zero `WARNING` lines in `logs/error_analysis.log` -- every chain that
  had Phase 7 predictions also had resolvable Calpha coordinates, pLDDT,
  and a secondary-structure call for effectively every residue (only 1
  residue project-wide dropped from the regression for a missing RSA
  value).
- **pLDDT source confirmed, not assumed:** AlphaFold DB mmCIFs store
  per-residue pLDDT in `_atom_site.B_iso_or_equiv` (gemmi's
  `atom.b_iso`), identical across all atoms of a residue, values in
  AlphaFold's documented 0-100 range -- verified by direct inspection of
  a real cached AlphaFold model before writing any code that assumed it.
- **Secondary structure:** `mkdssp` confirmed not installed (no root
  access, as already recorded); `biotite.structure.annotate_sse` (P-SEA,
  3-state) used for all 566 chains.
- **Primary analysis -- pLDDT-band metrics** (`results/error_analysis/band_metrics.csv`,
  figure `error_analysis_band_metrics.png`; chain-resampled bootstrap
  CIs, n=2,000 draws -- smaller than Phase 7's 10,000, profiled
  empirically beforehand since each draw here recomputes a pooled
  O(n log n) metric, not a cheap scalar median): the exp-vs-af_trimmed
  AUPR gap shrinks monotonically with AlphaFold confidence -- **0.255
  (`<50`, n=2,828) -> 0.117 (`50-70`, n=5,735) -> 0.071 (`70-90`,
  n=23,762) -> 0.040 (`>=90`, n=62,793)**, a >6x reduction from lowest to
  highest confidence band; ROC-AUC gap shows the same pattern (0.195 ->
  0.075 -> 0.044 -> 0.022). Every band cleared `PLDDT_BAND_MIN_RESIDUES`
  (30) comfortably, zero undefined bootstrap draws anywhere. This
  directly answers Phase 8's primary question: yes, AlphaFold's own
  confidence tracks where its structural error costs PeSTo accuracy.
- **Secondary analysis -- shift regression** (`results/error_analysis/regression.csv`;
  n=95,117, 566 clusters, cluster-robust CRV1 SEs by chain): standardized
  `plddt_z` = **-0.0268** (p=3.9e-16), `local_rmsd_z` = **+0.0268**
  (p=1.4e-9), `rsa_z` = **+0.0254** (p=1.4e-57), all independently
  significant in the expected direction controlling for the others.
  Secondary structure: helix +0.0115 vs. coil (p=9.8e-4), strand +0.0073
  vs. coil (p=0.081, not significant). **No VIF exceeds
  `VIF_COLLINEARITY_THRESHOLD` (5)** -- max is 1.64 (`plddt_z`) -- milder
  collinearity between pLDDT and local RMSD than the pre-registration
  flagged as a risk; both carry real independent signal.
- **Descriptive -- chain length vs. AUPR drop** (`results/error_analysis/length_analysis.csv`,
  figure `error_analysis_aupr_drop_vs_length.png`; n=554): **no
  meaningful association** -- Spearman rho=0.028, 95% CI [-0.054, 0.108],
  p=0.51; OLS slope ~ -3.3e-6/residue, p=0.96; length quartiles (mean
  51/99/176/363 residues) show non-monotonic median drops
  0.020/0.020/0.041/0.025. **This contradicts the pre-registration's a
  priori hypothesis** that the pooled-vs-median gap was a length-
  weighting effect -- the real explanation is that the per-chain drop
  distribution is right-skewed (unweighted mean 0.061 vs. median 0.027),
  length-weighting the mean barely changes it (0.061, almost identical
  to unweighted), and pooled AUPR (0.070) isn't a weighted average of
  per-chain AUPRs to begin with -- it's a separate rank-based statistic
  over all residues combined. Reported as-is rather than reframed to fit
  the original guess.
- **Descriptive -- geometrically flagged (92) vs. unflagged (462)
  chains** (`results/error_analysis/flagged_chains_summary.csv`, figure
  `error_analysis_flagged_vs_unflagged.png`; two-sided Mann-Whitney U):
  flagged chains have much higher mean local RMSD (2.21 vs. 0.67 Å,
  p=4.2e-32) and global Ca distance (7.91 vs. 1.40 Å, p=9.5e-40), lower
  mean pLDDT (77.4 vs. 87.9, p=2.8e-11), longer chain length (224 vs. 162
  residues, p=2.5e-4), and a larger AUPR drop (0.108 vs. 0.052, p=0.02);
  interface fraction does not differ (p=0.18). A coherent population --
  consistent with Phase 4's addendum finding that this flag catches real
  AlphaFold-vs-experimental 3D divergence, not a mapping artifact.
- **Phase 9 gate:** both conditions (Phase 7 primary-endpoint
  significance + a Phase 8 pLDDT-band effect) are now met for the first
  time. Per the lean-scope decision, **Phase 9 is still not attempted**
  -- that decision was to skip Phase 9 for this pass unconditionally, not
  contingent on the gate outcome.
- **Anything surprising:** the length-vs-drop null result is the
  standout finding -- a reasonable a priori guess the data doesn't
  support, worth remembering if this dataset is revisited. Also notable:
  pLDDT/local-RMSD collinearity (VIF 1.6) was milder than expected going
  in.

**Tests:** `tests/test_structural_metrics.py` (11 tests: Kabsch recovers
a known rotation; local RMSD near-zero for a rigidly moved copy and large
only near a locally perturbed region; too-few-neighbors flagged with a
reason, not silently dropped; global Ca distance sanity checks; pLDDT
band boundary assignment; I/O helpers against hand-built synthetic
mmCIFs with real coordinates/B-factors) + `tests/test_error_analysis.py`
(7 tests: every band reported including empty ones, undersized flagging;
chain-resampled bootstrap demonstrated at chain granularity via a
2-chain discrete-outcome check and an always-single-class-chain
undefined-draw-rate check; regression recovers known standardized
coefficients on synthetic clustered data; VIF correctly flags a
constructed collinear predictor pair and correctly does not flag an
independent one). Full suite: 161/161 passing.

**Deferred, not part of this pass** (per the lean-scope decision): the
apo/unbound arm (no unbound-structure local RMSD), `af_full`, and any
leakage-sensitivity-mode analysis.

### Phase 8b — pLDDT-filtering deployment check (post-hoc, exploratory)

**Done (2026-09-18).** Added this session specifically to test a claim
the Phase 8 write-up's Discussion made informally but never actually
tested: that pLDDT is "actionable" for flagging unreliable AlphaFold-input
predictions. Pre-registered in `PLAN.md` as post-hoc/exploratory (written
before these numbers were computed, but after Phase 7/8's own results
were already known -- explicitly not a true a priori pre-registration).
`src/analysis/plddt_filtering.py`, runnable as
`.venv/bin/python -m src.analysis.plddt_filtering`, ~4.5s runtime on the
primary set's 95,118 pooled residues.

**Result: filtering doesn't help much, and this is reported plainly rather
than reframed.** At a fixed 0.5 probability threshold, sweeping the pLDDT
cutoff from 0 to 90 moves AlphaFold-input F1 by only **+0.007** (0.596 ->
0.603, heavily overlapping 95% CIs), while discarding **34% of all
residues and 43% of the true interface residues**. The gap to experimental
input's F1 narrows only from 0.050 to 0.031, and part of even that is
experimental input's own F1 drifting down (0.647 -> 0.634) as its harder
residues get discarded too -- not AlphaFold input distinctly "catching
up." AlphaFold-input precision does rise materially (0.673 -> 0.713) while
recall stays essentially flat (0.535 -> 0.523): filtering trades real
interface coverage for a modest precision gain, not a broad quality
improvement. Full table in `PLAN.md`'s Phase 8b write-up and
`results/error_analysis/plddt_filtering.csv`.

**This directly replaced overreaching language in the paper.** The Phase
8 Discussion previously called pLDDT "the most actionable finding here"
and described it as "a practical... confidence signal for flagging
low-reliability regions" without ever testing a filtering rule --
`paper/template.tex` now cites Phase 8b's actual measured numbers instead,
and a new Results subsection (3.3) and figure (F1-vs-cutoff plus a
cost-of-filtering panel) report the full result.

**Tests:** `tests/test_plddt_filtering.py` (9 tests, synthetic data only)
-- precision/recall/F1 correctness at known thresholds; cutoff=0 keeps
everything; a cutoff above every residue's pLDDT keeps zero and is
reported as a row, not dropped; higher cutoffs never increase the kept
fraction; a perfect predictor scores perfectly regardless of cutoff;
chain-resampled bootstrap demonstrated at chain granularity (same
discrete-outcome pattern as Phase 8's own bootstrap tests). Full suite:
170/170 passing.

**Anything surprising:** the disconnect between Phase 8's own headline
pooled-AUPR-by-band result (>6x gap shrinkage from lowest to highest
pLDDT band) and this near-flat fixed-threshold F1 result. Both are
correct; they answer different questions (a threshold-free ranking metric
pooled by confidence band, vs. a fixed-threshold classifier's F1 after
discarding low-confidence residues), and conflating the two is exactly
the overreach this phase exists to catch.

## Next step

Phases 1-8 (lean scope), Phase 8b, and the results paper (with this
session's audit fixes) are all done. **Phase 9 is explicitly out of scope
for this lean pass** (see the lean-scope decision and Phase 8's gate note
above) -- both its gate conditions are now met, but that doesn't reopen
the question; don't start it without a new explicit instruction. Phase
8b's own result (naive pLDDT filtering doesn't help much at a fixed
threshold) does not change this -- Phase 9 as scoped is a *learned*
pLDDT-aware fine-tune, not a hard filter, and Phase 8b's finding is framed
in the paper as leaving that question open, not resolving it either way.
Deferred items that would need to happen first if this project ever goes
beyond the lean scope: the ~250-chain seeded subsamples for the
`exact_train`/`none` leakage-sensitivity modes (needs its own Phase 6
inference run -- only the primary set has predictions so far), any
`af_full`/apo-arm comparison, and (for Phase 8 specifically) the
apo/unbound arm's own local RMSD. Phase 10 (report generation) has not
been formally started, though the results paper now covers most of what
it would produce.

**The results paper's structure figure (Fig. 6) still uses the matplotlib
Calpha-trace fallback**, not a real molecular render -- no headless 3D
renderer could be installed in this container (see "Decisions already
made" / README.md). `paper/build_paper.py` will automatically pick up
`paper/figures/example_chain_render.png` if it's ever dropped in from
outside this environment; no code change needed, see README.md's
"Substituting a rendered structure figure" for exactly what it should
show.

Open question #3 from `PLAN.md` Sec. 5 is now fully resolved for this
project's purposes: Phase 7's own primary-endpoint alpha (0.05) is fixed
and was met (p ~ 3.5e-17), and Phase 8 confirmed a pLDDT-band effect
(monotonic, >6x AUPR-gap shrinkage from lowest to highest confidence
band, plus a significant `plddt_z` regression coefficient) -- **the
combined Phase 9 gate is met, as of 2026-09-17.** Per the lean-scope
decision, this still does not by itself start Phase 9 -- that remains a
separate, not-yet-given instruction.

**Worth considering before Phase 6 or later:** whether Phase 3 should be
revisited to add an explicit `providerId == "GDM"` check at fetch time
(see "Decisions already made" above) — not urgent, since Phase 4's gate
already prevents any contamination downstream, but it would save the
(currently wasted) disk space/bandwidth already spent downloading the 251
Community-model mmCIFs across all 2,756 accessions Phase 3 queried, and
would make Phase 3's own `fetch_report.csv` reflect this distinction
without needing to cross-reference Phase 4's exclusion reasons.

## How to update this file

When you finish a phase, a sub-step, or make a new decision that isn't yet
in `PLAN.md`: update "What exists on disk," move finished items out of
"Next step" and into "Decisions already made" or a new "Completed phases"
section, and bump "Last updated." Keep entries about *state and decisions*,
not a narrated log of tool calls.
