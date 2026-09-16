# Implementation Plan

Source: `docs/Project_Proposal.pdf` (extracted text: `docs/proposal.txt`).
All thresholds referenced below live in `src/config.py`; no phase should
hardcode a number that already exists there.

## 1. Project question

PeSTo is a geometric transformer trained on experimental PDB structures to
label protein-protein interface residues from atomic coordinates alone. This
project measures how much of that accuracy survives when PeSTo is instead
given AlphaFold-predicted structures of the same proteins, using a paired
design (same ground-truth interface, same protein, only the input structure
changes) so any performance drop can be attributed to structural input
quality rather than protein difficulty. It further asks *why* the drop
happens — testing whether pLDDT and local structural deviation from the
experimental structure explain prediction errors — and, if they do, whether
adding pLDDT as a model input (fine-tuned on Colab, since this container has
no GPU) recovers some of the lost accuracy.

## 2. Phases

### Phase 1 — Candidate complex selection
- **Inputs:** RCSB Search API (REST, live query); SIFTS PDB↔UniProt segment
  mapping (bulk file, e.g. `pdb_chain_uniprot.tsv.gz` from EBI).
- **Config used:** `RESOLUTION_CUTOFF_ANGSTROM` (2.5), `PDB_RELEASE_DATE_CUTOFF`
  (2018-04-30) — see confound (a); candidates are filtered to
  `initial_release_date > PDB_RELEASE_DATE_CUTOFF` (RCSB
  `rcsb_accession_info.initial_release_date`, **not** deposit date — AF2's
  training cutoff is a release-date cutoff) to reduce AF2 memorization risk.
- **Outputs:**
  - `data/raw/rcsb_search/query_<hash>.json` (cached raw API response)
  - `data/raw/sifts/pdb_chain_uniprot.tsv.gz` (cached bulk file, chain-level
    lookups only — see Phase 4 for residue-level mapping)
  - `data/interim/candidate_complexes.parquet` — columns: `pdb_id`, `chain_id`,
    `entity_id`, `uniprot_acc`, `resolution`, `initial_release_date`,
    `n_protein_entities`, `has_nucleic_acid`
- **src modules:** `src/data/rcsb_search.py`, `src/data/sifts.py`,
  `src/pipeline/select_candidates.py`
- **Tests:** `tests/test_rcsb_search.py` (query JSON embeds config cutoffs,
  not literals; mocked response parsed correctly), `tests/test_sifts.py`
  (mapping parser on a small fixture), `tests/test_select_candidates.py`
  (entries with a nucleic-acid entity or <2 distinct protein entities are
  excluded and the drop is logged with a reason)
- **Success criterion:** query executes, cache file is non-empty,
  100% of retained rows satisfy resolution ≤ 2.5 Å, release date >
  2018-04-30, ≥2 UniProt-mapped protein chains, 0 nucleic-acid entities.
  Actual retained count is logged; if it's under ~200 complexes, flag before
  continuing to Phase 2 (see confound (e)).

### Phase 2 — Redundancy reduction + PeSTo-overlap flagging (sequence-level)

**Resolving the train/validation/test naming ambiguity.** Downloaded all
three PeSTo split files directly (now cached at `data/raw/pesto_splits/`)
and cross-checked three independent sources of evidence:

| file | line count | matches paper's stated count/role? | matches training-code usage? |
|---|---|---|---|
| `subunits_train_set.txt` | 376,216 | yes — paper: "training set (376216 chains)" | yes — `model/config.py`'s `train_selection_filepath`, fed to `dataloader_train` for gradient updates in `model/main.py` |
| `subunits_test_set.txt` | 97,424 | paper says "testing set (97424 chains)," used only for final quality assessment | **no** — it's `model/config.py`'s `test_selection_filepath`, and `model/main.py` evaluates `dataloader_test` (built from this exact file) at every `eval_step` and overwrites the single released `model.pt` checkpoint whenever its loss on this set improves ("save model and update min loss") — i.e. in the actual public training code, this file is used for **checkpoint/model selection**, not held out |
| `subunits_validation_set.txt` | 101,700 | paper says "validation set (101700 chains)," used for hyperparameter selection during training | **no** — this filename is **never referenced** anywhere in `model/config.py` or `model/main.py`; per `data/datasets/README.md` it's used only afterward "to evaluate the trained model and perform comparisons/benchmarks" |

The **training set identification is unambiguous**: `subunits_train_set.txt`
is by far the largest file, its count matches the paper exactly, and the
public training script loads it directly for gradient updates.

The **test-vs-validation identification is not confident**: the paper's
Methods text and the actual public training code disagree about which
named file was used for in-training model selection. Per the fallback
agreed in advance, **both `subunits_test_set.txt` and `subunits_validation_set.txt`
are treated as excluded** (unioned with the training set) rather than
assuming either one is a clean, untouched holdout — no chain is treated as
"only overlapping PeSTo's real test set and therefore still eligible,"
because we cannot say with confidence which file that is.

- **Inputs:** `data/interim/candidate_complexes.parquet`; chain sequences
  from the RCSB polymer-entity API; sequences for every chain in
  `data/raw/pesto_splits/subunits_{train,test,validation}_set.txt` (all
  three, per the resolution above), fetched from RCSB and cached.
- **Config used:** `SEQUENCE_IDENTITY_CUTOFF` (0.30).
- **Outputs:**
  - `data/interim/sequences.fasta` (candidate chains)
  - `data/interim/redundancy_clusters.tsv` (candidate chain → cluster ID,
    for the 30%-identity redundancy reduction among candidates themselves)
  - `data/raw/pesto_splits/subunits_{train,test,validation}_set.txt`
    (already cached — verified present, see table above)
  - `data/raw/pesto_splits/pesto_all_splits.fasta` (sequences for the union
    of all three PeSTo split files)
  - `data/interim/pesto_homology_search.tsv` (MMseqs2 search hits: candidate
    chain → matched PeSTo chain(s), % identity)
  - `data/processed/nonredundant_complexes.parquet` — adds `cluster_id`,
    `cluster_representative` (bool), `pesto_homolog_overlap` (bool — hit at
    ≥`SEQUENCE_IDENTITY_CUTOFF` against the train+test+validation union),
    `pesto_exact_train_overlap` (bool — exact `PDBID_CHAINID` match against
    `subunits_train_set.txt` only, secondary/narrower flag)
- **src modules:** `src/data/sequences.py`, `src/pipeline/redundancy_reduction.py`,
  `src/data/pesto_overlap.py` (runs MMseqs2 `search` of candidate sequences
  against `pesto_all_splits.fasta` at `SEQUENCE_IDENTITY_CUTOFF`; also does
  the exact-ID lookup against `subunits_train_set.txt` alone)
- **Tests:** `tests/test_sequences.py`, `tests/test_redundancy_reduction.py`
  (clustering on synthetic sequences with known identity, verified against
  `SEQUENCE_IDENTITY_CUTOFF`), `tests/test_pesto_overlap.py` (a synthetic
  candidate sequence built as a point mutant of a fixture PeSTo-split
  sequence, at identity above/below the cutoff, is flagged/not-flagged
  correctly for `pesto_homolog_overlap`; exact-ID match logic is tested
  separately for `pesto_exact_train_overlap`)
- **Success criterion:** no two retained candidate chains have pairwise
  identity above the config cutoff; `pesto_homolog_overlap` and
  `pesto_exact_train_overlap` are both non-null for 100% of rows; MMseqs2
  search log records how many candidate chains hit the PeSTo split union
  and at what identity; dropped-duplicate count logged with reason.

### Phase 3 — Download PDB mmCIFs and AlphaFold DB models
- **Inputs:** `data/processed/nonredundant_complexes.parquet`.
- **Outputs:**
  - `data/raw/pdb/{pdb_id}.cif.gz`
  - `data/raw/alphafold/{uniprot_acc}.cif.gz` (AlphaFold DB v4)
  - `data/interim/download_manifest.parquet` — `id`, `source`, `url`,
    `local_path`, `status`, `reason`
- **src modules:** `src/data/download_pdb.py`, `src/data/download_alphafold.py`,
  `src/pipeline/download_structures.py`
- **Tests:** `tests/test_download_pdb.py`, `tests/test_download_alphafold.py`
  (mocked HTTP; cached files are not re-fetched), `tests/test_download_structures.py`
  (manifest correctly tallies success/failure)
- **Success criterion:** ≥95% of nonredundant entries have both files cached;
  a rerun with files already present makes 0 HTTP calls (asserted via mock
  call count); every failure logged with ID + HTTP status/reason.

### Phase 4 — Residue-numbering mapping (most correctness-critical)
- **Inputs:** PDBe's "updated" mmCIF per entry — `https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id_lower}_updated.cif`
  — used **instead of** the plain RCSB mmCIF from Phase 3 for this step, since
  it embeds per-residue SIFTS UniProt cross-references directly in extra
  `_atom_site` fields (`pdbx_sifts_xref_db_name`, `pdbx_sifts_xref_db_acc`,
  `pdbx_sifts_xref_db_num`, `pdbx_sifts_xref_db_res`) alongside
  `_pdbx_poly_seq_scheme` — one download gets both structure and mapping, no
  per-entry mappings-API call needed. The bulk `pdb_chain_uniprot.tsv.gz`
  from Phase 1 still covers chain-level lookups; this phase is what resolves
  residue-level numbering.
- **Outputs:**
  - `data/raw/pdb_updated/{pdb_id}_updated.cif.gz` (cached)
  - `data/interim/residue_mappings/{pdb_id}_{chain}.parquet` —
    `auth_seq_id`, `auth_ins_code`, `label_seq_id`, `uniprot_acc`,
    `uniprot_resnum`, `alphafold_resnum` (= `uniprot_resnum`, AF model numbering
    is UniProt-native), `residue_name`, `is_observed`, `unmapped_reason`
- **src modules:** `src/data/download_pdb_updated.py`, `src/mapping/pdb_numbering.py`,
  `src/mapping/sifts_residue_map.py` (parses the embedded `_atom_site`
  `pdbx_sifts_xref_db_*` fields, not a separate API response),
  `src/pipeline/build_residue_mappings.py`
- **Tests:** `tests/test_pdb_numbering.py` (fixture mmCIF with an insertion
  code and a gap, parsed exactly), `tests/test_sifts_residue_map.py` (parses
  the embedded `pdbx_sifts_xref_db_*` columns from a fixture updated-mmCIF),
  `tests/test_build_residue_mappings.py` (round-trip against a hand-verified
  table for ≥5 real entries chosen to include insertion codes, engineered
  tags, and chain breaks; every unmapped residue carries a reason string,
  never a silent drop)
- **Success criterion:** 100% match against the hand-verified fixture table;
  every observed `auth_seq_id` in every mapped chain has either a UniProt
  residue number or an explicit `unmapped_reason`.

### Phase 5 — Ground-truth interface labels
- **Inputs:** `data/raw/pdb/{pdb_id}.cif.gz` (full assembly), residue mappings
  from Phase 4.
- **Config used:** `INTERFACE_DISTANCE_CUTOFF_ANGSTROM` (5.0).
- **Outputs:** `data/processed/interface_labels/{pdb_id}.parquet` —
  `chain_id`, `auth_seq_id`, `uniprot_resnum`, `is_interface_contact`,
  `delta_sasa`, `is_interface_sasa`
- **src modules:** `src/labels/heavy_atom_contacts.py` (cross-chain heavy-atom
  neighbor search via biotite/gemmi), `src/labels/sasa_labels.py`
  (`biotite.structure.sasa`, isolated vs. complexed), `src/pipeline/build_interface_labels.py`
- **Tests:** `tests/test_heavy_atom_contacts.py` (synthetic two-residue system,
  boundary case at exactly 5.0 Å), `tests/test_sasa_labels.py` (synthetic
  dimer/monomer sanity check), `tests/test_build_interface_labels.py` (label
  agreement rate is computed and logged, not assumed)
- **Success criterion:** labels for 100% of successfully-parsed complexes;
  ≥85% residue-level agreement (Cohen's κ reported) between the 5 Å-contact
  and ΔSASA labels, disagreement rate logged.

### Phase 6 — PeSTo inference on isolated single chains
- **Inputs:** `external/pesto/` (pretrained weights, git checkout — gitignored,
  cached so reruns work offline); experimental structures, AlphaFold models,
  and unbound/apo structures (matched via a second RCSB query for
  single-entity entries sharing a UniProt accession with a candidate).
  The unbound/apo arm is **best-effort**: it does not block the main
  experimental-vs-AlphaFold benchmark, and its sample size is reported
  separately in Phase 7 rather than folded into the primary comparison.
- **Outputs:** `results/predictions/pesto/{experimental,alphafold,unbound}/{id}_{chain}.parquet`
  — per-residue predicted interface probability;
  `results/predictions/unbound_coverage.json` — count/fraction of candidate
  proteins for which an unbound structure was found, logged regardless of
  how small.
- **src modules:** `src/inference/isolate_chain.py`, `src/inference/pesto_runner.py`,
  `src/pipeline/run_pesto_inference.py`
- **Tests:** `tests/test_isolate_chain.py` (isolation removes all other
  chains' atoms, preserves numbering), `tests/test_pesto_runner.py` (smoke
  test on a tiny fixture: one score per residue, all scores in [0,1])
- **Success criterion:** predictions for 100% of chains with a valid isolated
  structure (or a logged failure reason) across the experimental and
  AlphaFold arms — these two are required; predicted-residue count equals
  the observed-residue count of the isolated input for every chain, no
  silent drop between mapping and inference. The unbound arm has no minimum
  yield requirement — `unbound_coverage.json` just needs to exist and be
  accurate.

### Phase 7 — Benchmarking
- **Inputs:** `results/predictions/pesto/*`, `data/processed/interface_labels/*`.
- **Outputs:** `results/metrics/per_protein_metrics.parquet` (`id`, `source`,
  `roc_auc`, `aupr`, `base_rate`, `n_residues`), `results/metrics/pooled_metrics.json`,
  `results/metrics/paired_tests.json` (Wilcoxon signed-rank, experimental vs.
  AlphaFold and experimental vs. unbound), `results/metrics/bootstrap_cis.json`
- **src modules:** `src/eval/metrics.py`, `src/eval/paired_tests.py`,
  `src/eval/bootstrap.py`, `src/pipeline/run_benchmark.py`
- **Tests:** `tests/test_metrics.py` (fixture with known ROC-AUC/AUPR),
  `tests/test_paired_tests.py` (Wilcoxon direction on synthetic paired
  deltas), `tests/test_bootstrap.py` (fixed-seed reproducibility)
- **Leakage sensitivity analysis.** The headline experimental-vs-AlphaFold
  benchmark is computed three ways, controlled by `LEAKAGE_FILTER_MODES` in
  `src/config.py` (`"homolog"`, `"exact_train"`, `"none"`), each written to
  its own row set in `results/metrics/per_protein_metrics.parquet` (a
  `leakage_filter_mode` column) and reported with its own n:
  - **Primary — `"homolog"`:** exclude chains with `pesto_homolog_overlap`
    (sequence-identity hit against the union of PeSTo's train + test +
    validation splits — see Phase 2/confound (b) for why the exclusion set
    is all three files, not just train+validation).
  - **Sensitivity A — `"exact_train"`:** exclude only `pesto_exact_train_overlap`
    (literal `PDBID_CHAINID` match to the training file alone) — a weaker,
    ID-only filter, kept to show how much the sequence-level search actually
    catches beyond exact-ID matching.
  - **Sensitivity B — `"none"`:** no leakage filter at all — upper-bound
    sample size, for reference only, not the reported headline number.
- **Success criterion:** metrics computed for every paired chain in every
  mode with zero unexplained NaNs; experimental-vs-AlphaFold AUPR delta
  reported with a Wilcoxon p-value and bootstrap CI per mode, fixed random
  seed recorded. Target ≥100 chains for the primary (`"homolog"`) mode; 50
  is a hard floor. Below 50, do not loosen the leakage filter to hit the
  target — instead report the primary mode's result stratified by
  `pesto_homolog_overlap` status (overlap vs. non-overlap side by side) so
  the leakage-affected comparison stays visible rather than disappearing,
  and lean on sensitivity modes A/B (which will have larger n) to show
  whether the result is filter-sensitive. The unbound arm's metrics (from
  `unbound_coverage.json`) are reported with their own (likely much smaller)
  n, separately, and are not required to clear any size threshold.

### Phase 8 — Error analysis
- **Inputs:** AlphaFold predictions/labels, AlphaFold mmCIF B-factor column
  (pLDDT), residue mappings, experimental (and unbound, where available)
  structures for local RMSD.
- **Config used:** `PLDDT_BANDS` (50, 70, 90).
- **Outputs:** `results/error_analysis/plddt_stratified_metrics.parquet`,
  `results/error_analysis/local_rmsd.parquet`, `results/error_analysis/feature_table.parquet`
  (`plddt`, `local_rmsd`, `rsa`, `secondary_structure`, `prediction_error`),
  `results/error_analysis/logistic_regression_summary.json` (coefficients + VIF)
- **src modules:** `src/analysis/plddt_stratify.py`, `src/analysis/local_rmsd.py`
  (windowed superposition), `src/analysis/rsa.py` (`biotite.structure.sasa` +
  Sander/Rost max-ASA table — see risk on `freesasa`), `src/analysis/secondary_structure.py`
  (`biotite.structure.annotate_sse` — see risk on DSSP), `src/analysis/logistic_model.py`
  (logistic regression + `statsmodels` VIF — **not yet installed**, add when
  this phase starts)
- **Tests:** `tests/test_plddt_stratify.py` (band edges read from config),
  `tests/test_local_rmsd.py` (identical structure → RMSD ≈ 0; known
  translation → RMSD matches offset), `tests/test_rsa.py`,
  `tests/test_secondary_structure.py`, `tests/test_logistic_model.py`
  (synthetic collinear feature pair is flagged)
- **Success criterion:** metrics reported for every non-empty pLDDT band
  (n≥30 residues) — empty bands logged, not skipped silently; logistic
  regression converges and any feature with VIF > 5 is flagged in the output,
  not silently included.

### Phase 9 — pLDDT-augmented fine-tuning (gated, environment-agnostic)

**Gate: this phase only runs if both hold:**
1. Phase 7's primary (`"homolog"`-mode) benchmark shows a statistically
   significant experimental-vs-AlphaFold AUPR gap (Wilcoxon p-value below a
   preregistered threshold, e.g. 0.05, on the bootstrap CI excluding 0).
2. Phase 8 shows that gap is associated with pLDDT (e.g. AUPR/ROC-AUC
   materially lower in the pLDDT<70 band than pLDDT>90, and/or a
   significant pLDDT coefficient in the Phase 8 logistic model after the
   VIF check).

If either condition fails, Phase 9 is skipped and the project's answer is
"there is no measurable/pLDDT-linked gap to correct" rather than fine-tuning
anyway. This gate is recorded here explicitly so it isn't silently
skipped or silently run regardless of Phase 7/8 outcomes.

- **Inputs:** experimental structures (pLDDT fixed at 100) + AlphaFold
  structures (real pLDDT) as a mixed training set; PeSTo pretrained weights.
- **Environment:** the training script is environment-agnostic (a plain
  CLI script taking `--device`), runnable on Colab GPU or an HPC cluster —
  which one is a later decision, not fixed by this plan.
- **Outputs:** `src/models/pesto_plddt.py` (widened input-embedding layer —
  see deviation (b)), `src/training/train_pesto_plddt.py` (CLI-driven, no
  local-GPU assumptions), `notebooks/finetune_pesto_plddt.ipynb` (thin Colab
  wrapper), `results/checkpoints/pesto_plddt_v1.pt` (produced on Colab, not
  in this container), `results/metrics/finetuned_benchmark.parquet`
- **src modules:** `src/models/pesto_plddt.py`, `src/training/dataset.py`,
  `src/training/train_pesto_plddt.py`
- **Tests:** `tests/test_pesto_plddt_model.py` (CPU forward pass on a tiny
  synthetic batch; embedding width = original width + 1), `tests/test_dataset.py`
  (experimental nodes get pLDDT=100, AlphaFold nodes get real values, checked
  on a fixture)
- **Success criterion:** the CPU smoke test loads PeSTo's actual released
  checkpoint (not random init) into the widened model — original weights
  copied into the corresponding slice of the widened embedding layer, new
  pLDDT channel's weights zero- or randomly-initialized — and confirms that
  feeding an experimental structure with pLDDT hardcoded to 100 reproduces
  the original pretrained PeSTo's output (within float tolerance) on a tiny
  fixture; this verifies the widened layer is a correctness-preserving
  superset of the original before any fine-tuning happens on Colab.
  `train_pesto_plddt.py --device cuda` should run unmodified on whichever
  GPU environment is chosen (Colab or HPC) on a small subset before scaling
  up there; which environment, and when, is a later decision gated on
  Phase 7/8 results as described above.

### Phase 10 — Reproducible report generation
- **Inputs:** everything under `results/metrics/` and `results/error_analysis/`.
- **Outputs:** `results/figures/*.png`, `results/report/report.html`
- **src modules:** `src/report/figures.py`, `src/report/build_report.py`
- **Tests:** `tests/test_figures.py` (each figure function runs on synthetic
  input and writes a non-empty file)
- **Success criterion:** `.venv/bin/python -m src.report.build_report`
  regenerates the report and all figures from only `data/processed/` and
  `results/metrics/`, with no manual steps, and is idempotent at a fixed seed.

## 3. Known confounds and mitigations

**(a) AlphaFold2 training-cutoff memorization.** AF2's original training set
used PDB entries *released* before 2018-04-30 (release date, not deposit
date — deposit can precede release by months); testing on complexes AF2 has
memorized isn't a real generalization test. *Mitigation:* Phase 1 filters
candidates to `initial_release_date > PDB_RELEASE_DATE_CUTOFF`. AlphaFold DB
v4 models come from AlphaFold2 (not AF3/Multimer retraining), so this cutoff
is treated as authoritative rather than a proxy — resolves former open
question #4.

**(b) PeSTo training-set overlap.** PeSTo was trained on "all biological
assemblies from the PDB," and post-2018 PDB entries (our whole candidate
pool, per confound (a)) are exactly the ones likely to include close
homologs of PeSTo's training chains deposited under a *different* PDB ID —
so **exact `PDBID_CHAINID` matching alone would systematically miss
leakage**, especially for our candidate set. *Mitigation:* Phase 2 does a
sequence-level search (MMseqs2 at `SEQUENCE_IDENTITY_CUTOFF`) of every
candidate chain against PeSTo's published split files (verified present at
`github.com/LBM-EPFL/PeSTo/data/datasets/`), flagging any hit as
`pesto_homolog_overlap`. A narrower `pesto_exact_train_overlap` (literal ID
match against the training file only) is kept as a secondary column.
Resolving *which* published file is the true training set vs. the true
held-out set required cross-checking the paper's Methods text against the
actual public training code (`model/main.py`, `model/config.py`) — see
Phase 2 for the full evidence table. Result: the training file is
unambiguous (`subunits_train_set.txt`, 376,216 chains, matches the paper's
count and is what the code loads for gradient updates), but the other two
files' roles are **not** confidently resolved — the paper's Methods text and
the training code disagree about which of `subunits_test_set.txt` (used by
the actual training loop for checkpoint selection) and
`subunits_validation_set.txt` (never referenced in the training code, used
only for post-hoc benchmarking per the repo's README) was the untouched
holdout. Per the agreed fallback, **both are treated as excluded** — the
homolog search is run against the union of all three files, not just
train+validation, and no chain is treated as "leakage-free" solely for
matching one of the two ambiguous files. Rather than silently dropping
overlapping entries, Phase 7 reports metrics separately for overlap-flagged
vs. non-overlap subsets (ties to confound (e)) and adds two sensitivity
variants (see Phase 7).

**(c) Bound-vs-unbound conformational change mistaken for AF error.** A
holo (bound) crystal structure can differ from AlphaFold's prediction simply
because the protein moves upon binding, not because AlphaFold is wrong.
*Mitigation:* Phase 6's unbound/apo arm and Phase 8's local RMSD are computed
against both bound and (where available) unbound experimental structures;
residues with large apo→holo motion are flagged so AF-vs-bound deviation at
those positions isn't attributed to AF error by default. This arm is
best-effort (see Phase 6) — its coverage is likely partial, and its
findings are reported with their own sample size rather than treated as
equally powered to the main benchmark.

**(d) Residue coverage mismatch.** Full-length AlphaFold models cover the
whole UniProt sequence; PDB chains are often partially observed (missing
loops, tags, chain breaks). *Mitigation:* Phase 7 metrics are restricted to
the residue set observed in the experimental structure AND mapped in
AlphaFold (built in Phase 4/5); the count of AlphaFold-only residues excluded
from scoring is logged per protein, not silently dropped from labels.

**(e) Test set too small after leakage filtering.** Excluding all
`pesto_homolog_overlap` entries (against the train+test+validation union,
per confound (b)) could leave too few complexes for a reliable paired test.
*Mitigation:* target ≥100 non-overlap chains for the primary
(`"homolog"`-mode) analysis; 50 is a hard floor. Below 50, do not loosen the
filter to hit the target — Phase 7 instead (1) reports the primary mode's
result stratified by `pesto_homolog_overlap` status side by side, and (2)
reports the two sensitivity modes (`"exact_train"`, `"none"`), which have
strictly larger n, so the reader can see whether the finding survives a
looser filter rather than losing the comparison entirely to a sample-size
floor.

## 4. Deviations from the proposal

1. **PISA → direct coordinate computation.** The proposal cites PISA for
   interface annotations; we compute interfaces directly from coordinates
   (5 Å heavy-atom contact, primary; ΔSASA, robustness check) instead.
   *Why:* PISA has no reliable batch/programmatic API for thousands of
   automated, cacheable, offline-rerunnable queries, and computing directly
   gives us a single explicit, config-driven distance cutoff instead of
   PISA's internal criteria.
2. **"No architectural changes" is inaccurate.** Adding pLDDT as a node
   feature changes the input dimensionality of PeSTo's first embedding
   layer, so that layer must be widened (and its pretrained weights either
   zero-padded/extended or reinitialized for the new channel). *Why:* the
   proposal likely means "no change to the overall graph-transformer
   architecture," but the literal claim is wrong at the input layer; PLAN.md
   Phase 9 treats this as a minimal, localized architectural change.
3. **Yuan et al.'s numbers are not a direct benchmark target.** Yuan et al.
   benchmarked PeSTo/ScanNet on ESMFold structures, not AlphaFold. *Why:* our
   AlphaFold-based AUPR drop is not directly comparable to their
   0.797→0.691 ESMFold figure; we'll report our own numbers and only compare
   direction (both expected to drop), not magnitude, against Yuan et al.
4. **"First systematic baseline" needs softening.** Yuan et al. already ran
   a systematic experimental-vs-predicted-structure benchmark of PeSTo (on
   ESMFold), and the proposal itself notes PeSTo's authors have begun their
   own AlphaFold analysis. *Why:* the accurate framing is that this project
   is the first systematic **AlphaFold-specific** benchmark with a matched
   bound/unbound/AlphaFold design plus pLDDT-stratified error analysis and a
   fine-tuning attempt — not the first to study transferability at all.

## 5. Open questions

Resolved since the first draft: PeSTo's split lists are public, and the
train/test/validation naming ambiguity is resolved as far as evidence
allows — training file is unambiguous, but test-vs-validation roles
conflict between the paper's Methods text and the actual public training
code, so both are conservatively excluded (→ confound (b), Phase 2); the
homolog search is sequence-level via MMseqs2, not exact-ID (→ Phase 2);
non-overlap target is ≥100 chains with a 50 floor, with two additional
sensitivity modes (→ confound (e), Phase 7); the unbound arm is best-effort
and non-blocking (→ confound (c)); the AF cutoff is 2018-04-30 by *release*
date, treated as authoritative for AF2/AFDB v4 (→ confound (a)); residue
mapping uses PDBe's updated mmCIF instead of per-entry API calls (→ Phase
4); Phase 9 now only runs if Phase 7 shows a significant gap and Phase 8
ties it to pLDDT, and its environment (Colab vs. HPC) is deferred until
that gate is reached. Remaining:

1. **Root cause of the paper/code discrepancy.** Is the published
   `model/main.py` actually the script that produced the distributed
   `model.pt` checkpoint (e.g. under `model/save/i_v4_1_.../`), or could an
   earlier/different training run (matching the paper's stated roles) have
   produced it? We took the code as ground truth for what happened, but
   haven't confirmed the released checkpoint's provenance against it
   directly — worth a quick check if the leakage numbers turn out to matter
   a lot for the headline result.
2. Should we also attempt to identify PeSTo's exact `bc-30.out` clustering
   snapshot date (not stated in the repo), to sanity-check whether any of
   our post-2018-04-30 candidates could even have existed at PeSTo's
   training time? Currently we rely purely on the sequence-level homology
   search (Phase 2), which doesn't need this date, so this is optional
   extra corroboration rather than a blocker.
3. What preregistered significance threshold should the Phase 9 gate use
   (this plan assumes Wilcoxon p<0.05 on the primary mode plus a
   pLDDT-band effect in Phase 8) — confirm before Phase 7 runs, so the gate
   isn't chosen after seeing the results.

## 6. Risks

- **MMseqs2** — not installed, no `apt`/root access in this container.
  Fallback: download a static MMseqs2 binary into `external/mmseqs/` (no
  root required); if that's blocked too, fall back to a slower pairwise
  Biopython-alignment clustering (only viable if the candidate set stays in
  the low hundreds after Phase 1).
- **DSSP (`mkdssp`)** — not installed, no root access. Fallback:
  `biotite.structure.annotate_sse` (confirmed available; 3-state P-SEA-style
  assignment instead of DSSP's 8-state) — document this as a resolution
  limitation in the Phase 10 report.
- **freesasa** — not installed as a system/pip package. Fallback:
  `biotite.structure.sasa` (Shrake-Rupley, confirmed available in the
  installed biotite) for both Phase 5 ΔSASA and Phase 8 RSA; note the
  algorithm choice in the report since absolute SASA values can differ
  slightly from freesasa's.
- **statsmodels** — needed for Phase 8's VIF collinearity check, not yet in
  `requirements.txt`; install into `.venv` when Phase 8 starts and re-freeze.
- **External API/repo throughput or downtime** (RCSB, PDBe, AlphaFold DB,
  and the PeSTo GitHub repo for Phase 2's split files) — mitigated by
  caching every response under `data/raw/` (CLAUDE.md convention) with
  retry/backoff and resumable, idempotent manifests, so reruns don't
  re-fetch and a mid-run outage doesn't lose progress. Using PDBe's updated
  mmCIF for Phase 4 (one download instead of structure + per-residue API
  call) reduces request volume specifically for the highest-count phase.
- **PeSTo repo could reorganize/rename `data/datasets/`** — low risk (files
  verified present at the time of this plan), but `pesto_overlap.py` should
  fail loudly (not silently skip flagging) if the expected files are
  missing after a fresh clone, since confound (b)'s mitigation depends on
  them.
- **No local GPU** — Phase 9 training cannot run locally by design; only a
  CPU correctness smoke test is feasible here, real training must happen on
  Colab.
