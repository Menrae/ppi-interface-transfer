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

### Phase 1 — Candidate complex selection — **DONE** (2026-09-16, full-scale rerun same day)

Implemented as a single module rather than the three sketched below —
`src/data/select_complexes.py` covers search, entry/entity fetch, SIFTS
cross-check, and filtering together, since none of the pieces were reused
elsewhere. Output paths and cache layout also ended up simpler than
originally sketched; both are corrected below to match what actually
exists on disk. First run used `--max-entries 500` as a pilot/smoke test
(preserved at `data/interim/pilot_500/`, not overwritten); the full,
uncapped run (all 10,774 search hits) is now the authoritative
`data/interim/candidates.csv`.

- **Inputs:** RCSB Search API v2 (REST, live query, paginated); RCSB
  GraphQL API (`data.rcsb.org/graphql`, batched entry+polymer-entity fetch
  — see below) with a REST single-item fallback; SIFTS PDB↔UniProt bulk
  file (`pdb_chain_uniprot.tsv.gz` from the EBI SIFTS FTP site).
- **Config used:** `RESOLUTION_CUTOFF_ANGSTROM` (2.5), `PDB_RELEASE_DATE_CUTOFF`
  (2018-04-30) — see confound (a); candidates are filtered to
  `initial_release_date > PDB_RELEASE_DATE_CUTOFF` (RCSB
  `rcsb_accession_info.initial_release_date`, **not** deposit date — AF2's
  training cutoff is a release-date cutoff) to reduce AF2 memorization risk;
  `MAX_PROTEIN_ENTITIES` (10) bounds entry complexity; `MIN_CHAIN_LENGTH`
  (40) drops peptide fragments (both added to `src/config.py` this phase).
- **Batched metadata fetch (added for the full-scale rerun):** entry +
  nested polymer-entity metadata is fetched via the RCSB GraphQL API in
  batches of `GRAPHQL_BATCH_SIZE` (200) IDs per request instead of one REST
  call per entry plus one per entity — empirically ~125 entries/sec, vs.
  tens of thousands of sequential REST calls under the old approach. The
  GraphQL response is split (`split_graphql_entry`) into the exact same
  per-entry/per-entity JSON shape the old REST fetchers cached, so
  `parse_entry_summary`/`parse_polymer_entity` and every downstream
  consumer (including Phase 2) are unchanged. Any ID GraphQL doesn't return
  (rare) falls back to the original single-item REST fetch rather than
  being silently dropped.
- **Resumability:** every entry/entity is still cached one-file-per-item
  (`entries/{id}.json`, `polymer_entities/{id}_{entity_id}.json}`) — an
  already-cached ID is never re-requested regardless of GraphQL batch
  boundaries, so an interrupted run resumes exactly where it left off, and
  a full rerun after any partial completion is cheap (cache hits only).
  Verified live: two transient connection resets during the full-scale
  search-pagination pass were retried transparently by `build_session()`'s
  Retry adapter with no data loss.
- **Deposition-group collapse (new attrition stage, entry-level):**
  PanDDA-style fragment-screening campaigns deposit many near-identical
  entries (same complex, different soaked fragment) under one RCSB
  "deposit group" — a real, documented field
  (`rcsb_entry_group_membership`/`pdbx_deposit_group`, filtered to
  `aggregation_method == "matching_deposit_group_id"` specifically, not
  RCSB's unrelated sequence-similarity browsing groups). Before the
  expensive per-entry fetch, one best-resolution entry is kept per group
  (`collapse_deposition_groups`, deterministic: lowest resolution, then
  ascending PDB ID); this both cuts fetch volume and keeps entry-level
  counts from misrepresenting biological diversity. Implemented as its own
  logged attrition stage (`deposition_group_collapse`, entry-level units)
  rather than left to Phase 2 clustering alone, once this authoritative
  RCSB field was found — Phase 2's chain-level MMseqs2 clustering still
  independently collapses any remaining near-duplicate sequences (e.g.
  campaigns not tagged with this field, or true independent redeposits of
  the same protein pair), so the two mechanisms are complementary, not
  redundant.
- **Outputs (actual):**
  - `data/raw/rcsb/search/<query_hash>/page_start<N>.json` (cached search
    pages, one per 100-row page)
  - `data/raw/rcsb/graphql_group_batches/`, `data/raw/rcsb/graphql_entry_batches/`
    (cached raw GraphQL batch responses)
  - `data/raw/rcsb/entries/{pdb_id}.json`, `data/raw/rcsb/polymer_entities/{pdb_id}_{entity_id}.json`
    (cached per-entry/per-entity metadata, written by either the GraphQL or
    REST-fallback path — same shape either way)
  - `data/raw/sifts/pdb_chain_uniprot.tsv.gz` (cached bulk file)
  - `data/interim/candidates.csv` — one row per surviving protein chain:
    `pdb_id`, `release_date`, `resolution`, `n_protein_entities`,
    `entity_id`, `chain_id`, `seq_length`, `polymer_type`, `uniprot_ids`
    (`;`-joined), `n_uniprot_ids`, `sifts_uniprot_ids`, `sifts_agrees`
  - `data/interim/phase1_attrition.csv` — entry-level stages
    (`entries_from_search`, `deposition_group_collapse`) followed by
    chain-level filter stages (`stage`, `n_before`, `n_dropped`,
    `n_remaining`, `note`; the `note` column says which unit a stage uses)
  - `data/interim/pilot_500/` — the preserved 500-entry pilot run's
    `candidates.csv`/`phase1_attrition.csv` (and its Phase 2 outputs), kept
    for before/after comparison, not part of the active pipeline
  - `logs/select_complexes.log` — full DEBUG-level log incl. every dropped
    chain's reason, every SIFTS disagreement, and every deposition-group
    duplicate drop
- **src modules:** `src/data/select_complexes.py` (search query builder,
  paginated search, GraphQL batched fetch + REST fallback, deposition-group
  collapse, SIFTS load+cross-check, chain filters, CLI).
- **Tests:** `tests/test_select_complexes.py` (21 tests, all passing,
  no network) — query builder reads config not literals; entry/entity
  parsing on fixture JSON in `tests/fixtures/`; SIFTS mapping load +
  agreement/disagreement cross-check; each chain filter individually
  (including a chimera case and a short-chain case) plus the combined
  pipeline's attrition bookkeeping; deposition-group collapse (best
  resolution, ID tiebreak, missing-resolution handling, singleton
  passthrough); GraphQL batch response splitting matches the REST fixture
  shape end-to-end through the real parsers; group-membership batch
  parsing filters out non-deposit aggregation methods; cache-warming
  resumability (a fully-cached ID triggers zero network calls, asserted via
  a fake session that raises on any unexpected call) and REST fallback for
  IDs GraphQL omits.
- **Success criterion:** met on both runs. Pilot (`--max-entries 500`,
  2026-09-16): 500/500 entries fetched, 1275 candidate chains, 160 unique
  UniProt accessions (preserved at `data/interim/pilot_500/`). Full run (no
  cap, same day): 9988/9988 post-collapse entries fetched (0 fetch errors,
  0 local sanity-check failures), 22887 candidate chains across 8808
  entries, 4029 unique UniProt accessions, total wall-clock ~100 seconds.
  Full attrition and distribution numbers are in `PROGRESS.md`.

### Phase 2 — Redundancy reduction + PeSTo-overlap flagging (sequence-level) — **DONE** (2026-09-16)

Implemented as a single module, `src/data/cluster_and_split.py`, rather than
the three sketched below (`src/data/sequences.py`,
`src/pipeline/redundancy_reduction.py`, `src/data/pesto_overlap.py`) —
following the same precedent as Phase 1: none of the pieces are reused
elsewhere, so one module covers sequence loading, MMseqs2 clustering,
PeSTo-split sequence lookup, and the homology search together. Output
filenames and the on-disk format (CSV/TSV throughout, not parquet) are also
corrected below to match what actually exists. The train/validation/test
evidence table immediately below is unchanged and still authoritative.

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

*(Minor footnote, confirmed 2026-09-16: `subunits_test_set.txt` and
`subunits_validation_set.txt` are each missing a trailing newline, so
`wc -l` undercounts their true chain count by exactly one --
`src/data/cluster_and_split.py` correctly loads 97,425 and 101,701 chains
respectively. This still matches the paper's stated 97,424/101,700, which
are themselves `wc -l`-style counts.)*

- **Inputs (actual):** `data/interim/candidates.csv` (Phase 1 output); chain
  sequences reused from Phase 1's cached RCSB polymer-entity JSON
  (`entity_poly.pdbx_seq_one_letter_code_can`) — no re-fetch needed; PeSTo
  split-file sequences resolved via the bulk `pdb_seqres.txt.gz` file from
  `files.wwpdb.org` (one ~67 MB download covering every PDB chain, cached),
  rather than per-chain RCSB calls, since the three split files together
  list ~575k chains and per-chain API calls at that volume were infeasible.
- **Config used:** `SEQUENCE_IDENTITY_CUTOFF` (0.30), `CLUSTER_MIN_COVERAGE`
  (0.80, new this phase), `MIN_TEST_CHAINS_TARGET` (100, new — formalizes
  the number already used in prose), `MIN_TEST_CHAINS_FLOOR` (50, new,
  likewise), `LEAKAGE_SWEEP_IDENTITY_THRESHOLDS` (0.30/0.50/0.70/0.95, new).
- **MMseqs2:** installed as the official static AVX2 Linux build at
  `external/mmseqs/bin/mmseqs` (gitignored, not conda) — commit
  `d401e78c2d18a822cdb1527d7464a043f6035a15`.
- **Outputs (actual):**
  - `data/interim/sequences.fasta` — candidate chain sequences
  - `external/mmseqs/` — cached MMseqs2 static binary (see above)
  - `data/raw/pdb_seqres/pdb_seqres.txt.gz` — cached bulk sequence file
  - `data/raw/pesto_splits/pesto_all_splits.fasta` — sequences for the
    union of all three PeSTo split files (chains with no entry in
    `pdb_seqres.txt.gz`, e.g. obsolete/superseded PDB IDs, are logged and
    excluded from the search target set, not silently merged in)
  - `data/raw/mmseqs_work/` — MMseqs2 intermediate cluster/search DBs and
    tmp dirs (cache, not a final output)
  - `data/interim/clusters.tsv` — `cluster_id, pdb_id, chain_id,
    is_representative`, one row per candidate chain (`cluster_id` is the
    deterministically-chosen representative's own `PDBID_CHAINID`, not
    MMseqs2's internal representative choice — see representative-selection
    rule below)
  - `data/interim/pesto_homology_search.tsv` — one row per cluster
    representative: max identity + best-hit ID against each of PeSTo's
    train/test/validation splits individually, plus `max_identity_any`
  - `data/interim/candidates_dedup.csv` — one row per cluster
    representative: all `candidates.csv` columns plus `cluster_id`,
    `cluster_size`, the per-split identity columns above,
    `pesto_homolog_overlap` (bool — `max_identity_any` ≥
    `SEQUENCE_IDENTITY_CUTOFF` against the train+test+validation union, the
    recorded leakage-control decision from the evidence table above),
    `pesto_exact_train_overlap` (bool — exact `PDBID_CHAINID` membership in
    `subunits_train_set.txt` only, secondary/narrower flag). This supersedes
    the originally-sketched `nonredundant_complexes.parquet`.
  - `data/interim/leakage_threshold_sweep.csv` — for each swept identity
    threshold plus an exact-ID-only mode: `n_survive`, `meets_target`
    (≥`MIN_TEST_CHAINS_TARGET`), `meets_floor` (≥`MIN_TEST_CHAINS_FLOOR`).
    This is a reporting table only — it does not change
    `SEQUENCE_IDENTITY_CUTOFF` or which chains are flagged in
    `candidates_dedup.csv`; see PROGRESS.md for the first run's numbers and
    what they imply for confound (e) / Phase 7.
  - `data/interim/phase2_attrition.csv` — chain-level attrition
    (`initial_candidate_chains` → `missing_cached_sequence` →
    `redundancy_clustering`)
- **Representative selection (deterministic, not MMseqs2's internal pick):**
  best (lowest) resolution, then longest chain (`seq_length`), then
  ascending `PDBID_CHAINID` string as a final tiebreak.
- **src modules:** `src/data/cluster_and_split.py` — sequence loading
  (reused from Phase 1's cache), `mmseqs easy-cluster` invocation +
  deterministic representative selection, bulk `pdb_seqres.txt.gz`
  streaming lookup for PeSTo split sequences, `mmseqs easy-search` of
  representatives against the PeSTo split union, leakage flag computation,
  and the threshold sweep.
- **Tests:** `tests/test_cluster_and_split.py` (14 tests, all passing) —
  deterministic representative selection (including tie-breaks and
  order-independence) on a hand-built fixture; parsing of saved MMseqs2
  `easy-cluster`/`easy-search` output fixtures in `tests/fixtures/`; leakage
  flag computation (best hit + identity per split, `pesto_homolog_overlap`,
  `pesto_exact_train_overlap`) on synthetic hits; threshold-sweep survivor
  counts against `MIN_TEST_CHAINS_TARGET`/`FLOOR` on a synthetic table; bulk
  `pdb_seqres.txt.gz` lookup on a synthetic gzip fixture; one test that
  invokes the real MMseqs2 binary end-to-end on tiny synthetic sequences,
  marked to skip if the binary is absent.
- **Success criterion:** met on both runs. Every retained representative has
  non-null `pesto_homolog_overlap`/`pesto_exact_train_overlap`; the MMseqs2
  command lines and hit counts are logged (`logs/cluster_and_split.log`);
  the dropped-duplicate count is logged with reason at the
  `redundancy_clustering` attrition stage. Pilot run (2026-09-16, against
  Phase 1's 500-entry pilot slice, preserved at `data/interim/pilot_500/`):
  the primary (`"homolog"`-mode) survivor count (15) came in well under
  `MIN_TEST_CHAINS_FLOOR` on that small, campaign-heavy candidate pool.
  Full-scale rerun (same day, unchanged code, run against the full
  22,887-chain candidate pool from Phase 1's uncapped run): the primary
  mode now clears both `MIN_TEST_CHAINS_TARGET` and `MIN_TEST_CHAINS_FLOOR`
  at every swept threshold, confirming the pilot's shortfall was a small-
  sample artifact rather than a structural problem. See PROGRESS.md for the
  full cluster count, leakage rates, threshold-sweep table, and the two
  additional breakdowns (train-only vs. train+test+validation union;
  survivor rate by release year) produced to inform the leakage-threshold
  decision. **Decision recorded 2026-09-16** (see confound (e) and the new
  confound (f) below): the primary test set is `identity < 0.30` against
  PeSTo's train+test+validation union — i.e. `eligible_homolog` /
  `pesto_homolog_overlap == False` in `candidates_dedup.csv` — n=696. The
  sensitivity analyses (looser thresholds up to 0.95, exact-ID-only, no
  filter) remain as recorded, all with strictly larger n. Every candidate
  representative's `candidates_dedup.csv` row now carries one boolean
  column per `LEAKAGE_FILTER_MODES` entry (`eligible_homolog`,
  `eligible_exact_train`, `eligible_none`), computed from the same
  in-memory leakage flags rather than recomputed from CSV strings, so later
  phases can select any subset by column without rerunning the homology
  search.

### Phase 3 — Download PDB mmCIFs and AlphaFold DB models — **DONE** (2026-09-16)

Implemented as a single module, `src/data/fetch_structures.py`, following
the same one-module-per-phase precedent as Phases 1/2. Two corrections to
the original sketch below, made deliberately at implementation time (this
session, after the leakage-threshold decision was recorded — see
PROGRESS.md):

1. **Scope is all 3,009 Phase 2 representatives, not a
   `nonredundant_complexes.parquet` that doesn't exist** (superseded by
   `candidates_dedup.csv`, see Phase 2). Every representative is
   downloaded — including the sensitivity-mode-only ones — so the
   sensitivity analyses (Phase 7) have inputs too, not just the primary
   696. Representatives are processed in priority order (primary
   `eligible_homolog` set first, then `eligible_exact_train`, then the
   rest) so the primary set is fully usable even if the run is
   interrupted partway through.
2. **The experimental structure downloaded here is PDBe's "updated" mmCIF**
   (`https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id_lower}_updated.cif`),
   not a plain RCSB mmCIF — i.e. this phase fetches directly what Phase 4
   was originally going to fetch separately into its own
   `data/raw/pdb_updated/` directory (see the corrected Phase 4 below).
   One download now serves both this phase's raw-structure need and Phase
   4's residue-mapping need (the file embeds SIFTS UniProt cross-references
   per residue), instead of two separate downloads of overlapping content.
- **Inputs:** `data/interim/candidates_dedup.csv` (Phase 2 output, 3,009
  representative rows, including the `eligible_homolog` /
  `eligible_exact_train` / `eligible_none` columns used for both download
  priority and the attrition breakdown below).
- **AlphaFold model resolution (no hardcoded version/URL pattern):** the
  AlphaFold DB prediction API (`https://alphafold.ebi.ac.uk/api/prediction/{accession}`)
  is queried per UniProt accession; the response can include *isoform*
  entries (e.g. querying `Q9UPN9` can also return `Q9UPN9-2`), so only the
  entry whose own `uniprotAccession` exactly matches the queried accession
  counts as a hit. Model version (`latestVersion`), the model's own UniProt
  sequence (`uniprotSequence`), and its residue range (`uniprotStart`/
  `uniprotEnd`) are all read from this response, never assumed. The mmCIF
  itself is downloaded from the response's own `cifUrl` — never a
  constructed/guessed URL. No PAE file is downloaded (not called for
  anywhere in this plan).
- **Outputs:**
  - `data/raw/pdb/{PDB_ID}_updated.cif.gz` (gzip-compressed; ~4x smaller
    than the plain-text download)
  - `data/raw/alphafold/{uniprot_acc}.cif.gz` (gzip-compressed)
  - `data/raw/alphafold/api/{uniprot_acc}.json` (cached raw prediction-API
    response — both successful and permanently-failed (e.g.
    accession-absent) outcomes are cached, since "absent from AlphaFold
    DB" is a common, stable outcome here worth memoizing, not a transient
    error to keep retrying)
  - `data/raw/uniprot/{uniprot_acc}.fasta` (cached current UniProt sequence,
    fetched to detect `sequence_mismatch` against the AlphaFold model's own
    sequence)
  - `data/interim/fetch_report.csv` — one row per representative:
    `cluster_id`, `pdb_id`, `chain_id`, `uniprot_acc`, `release_date`, the
    three `eligible_*` leakage-mode columns (carried through from Phase 2
    so downstream phases can filter without rejoining), both sources' local
    paths/status/failure-reason, and `overall_status`
    (`complete`/`partial`/`failed`)
  - `data/interim/phase3_attrition.csv` — one row per `LEAKAGE_FILTER_MODES`
    entry: eligible-representative count, PDB/AlphaFold/both-success counts
    for that mode's subset, and whether `n_both_success` clears
    `MIN_TEST_CHAINS_TARGET`
- **Pre-flight estimate (this session's addition, not in the original
  sketch):** before the bulk download, a small *real* (not synthetic)
  calibration sample of not-yet-cached ids is fetched to empirically
  project total wall-clock time and disk usage (`calibrate_and_estimate`
  in `src/data/fetch_structures.py`); if the projection exceeds
  `config.PHASE3_TIME_BUDGET_SECONDS` (2h) or `config.PHASE3_DISK_BUDGET_BYTES`
  (20GB), the run stops before downloading the rest. See PROGRESS.md for
  this run's actual calibrated projection and outcome.
- **Failure reasons, logged distinctly (never a silent drop):**
  `accession_absent_from_alphafold_db`, `fragmented_in_alphafold_db` (AFDB
  has split the *queried* accession's own sequence into >1 fragment model
  — fragments are not stitched here), `sequence_mismatch` (the AlphaFold
  model's own UniProt sequence disagrees with the accession's *current*
  UniProt sequence — numbering from this model would be unsafe for Phase
  4, so its mmCIF is not downloaded), `network_error` (transient —
  uncached, retried on the next run), `pdb_download_error:<HTTP status>`
  (PDBe returned a non-200 status for the updated mmCIF).
- **src modules:** `src/data/fetch_structures.py` — pure URL-builder/
  response-classifier functions (unit-tested directly on fixtures, no
  network) plus the cached fetch/download functions and orchestration.
- **Tests:** `tests/test_fetch_structures.py` (30 tests, no network) — URL
  construction; AFDB response classification on saved fixtures (success,
  isoform-only, fragmented, absent, empty body); FASTA parsing and
  sequence-match comparison; priority ordering; `candidates_dedup.csv`
  loading/boolean-parsing; fetch-report-row schema and
  `overall_status` derivation; per-leakage-mode attrition counting; and,
  via a fake `requests.Session` that raises on any unwired URL, both "cache
  hit makes zero network calls" and every failure-classification path
  (network error, HTTP error, accession absent, fragmented, sequence
  mismatch) end-to-end through `fetch_pdb_structure`/`fetch_alphafold_model`.
- **Success criterion:** every representative gets a `fetch_report.csv` row
  with a definite status and, on failure, a specific reason (never blank);
  a rerun with every file already cached makes zero HTTP calls; the
  per-leakage-mode attrition table reports whether the primary mode's
  `n_both_success` clears `MIN_TEST_CHAINS_TARGET`. See PROGRESS.md for this
  run's actual counts, disk usage, and failure-reason breakdown.

### Phase 4 — Residue-numbering mapping (most correctness-critical) — **DONE** (2026-09-16)

Implemented as a single module, `src/data/align_residues.py`, per this
session's explicit instruction (superseding the three-module sketch
below, following the same one-module-per-phase precedent as Phases 1-3).
**Deviation from the recorded output path, flagged and resolved in favor
of the recorded version per standing instruction:** this session's brief
asked for `data/processed/residue_maps/{pdb_id}_{chain}.parquet`; the
recorded plan (immediately above, prior revision) says
`data/interim/residue_mappings/{pdb_id}_{chain}.parquet`. Proceeded with
the **recorded** `data/interim/residue_mappings/` path — the actual
directory and columns produced are documented below.

- **Inputs:** `data/interim/fetch_report.csv` rows with
  `overall_status == "complete"` (2,704 representatives with both
  structures downloaded, Phase 3) — sufficient on its own (already carries
  `pdb_cif_path`, `alphafold_cif_path`, `alphafold_uniprot_start/end`,
  `uniprot_acc`, and the `eligible_*` leakage columns), no rejoin against
  `candidates_dedup.csv` needed. PDBe's "updated" mmCIF is read directly
  from Phase 3's cache at `data/raw/pdb/{PDB_ID}_updated.cif.gz` — no
  separate download.
- **Primary mapping source:** the same embedded `_atom_site` SIFTS
  cross-reference columns recorded in the prior revision of this section
  (`pdbx_sifts_xref_db_name/_acc/_num`) — **discovery this session: the
  AlphaFold mmCIF (fetched from AlphaFold DB via the EBI/PDBe mirror) embeds
  the identical columns**, self-declaring its own UniProt-native numbering
  per residue. One parser (`load_raw_atom_site_rows`/`reduce_to_residues`)
  handles both files. Modified residues (e.g. MSE) are resolved to their
  parent amino acid via `gemmi.find_tabulated_residue(...).one_letter_code`
  (covers the full CCD, not a hand-rolled table) — this is how they're
  "treated through their parent residue."
- **AlphaFold numbering is verified per chain, not assumed:**
  `validate_alphafold_numbering` checks that the AlphaFold model's own
  `auth_seq_id` range exactly equals its recorded
  `[alphafold_uniprot_start, alphafold_uniprot_end]` *and* that every
  residue's own embedded SIFTS annotation agrees (`db_name == "UNP"`,
  `db_acc == uniprot_acc`, `db_num == auth_seq_id`). Any disagreement (or
  the columns being entirely absent from that particular file) excludes
  the chain (`alphafold_numbering_not_uniprot_native`/
  `alphafold_range_mismatch`) rather than trusting the Phase 3 API
  metadata blindly.
- **Surprising finding this session, now an explicit gate:** ~9% of the
  UniProt accessions Phase 3 resolved to an "AlphaFold DB" hit are actually
  third-party **Community submissions** (mostly `ColabFold v1.5.2`/`v1.0-alpha`,
  provider ids `VR3D`/`ATBC`/`NTDX`/`BFVD`), not genuine Google DeepMind
  AlphaFold2 predictions (provider id `GDM`) — AlphaFold DB's own prediction
  API (queried by Phase 3) returns these under the *same* accession lookup,
  and Phase 3's spec never checked provider identity. These do not carry
  the training-cutoff/pipeline guarantee this project's whole AF2
  comparison rests on (confound (a)), so `verify_official_alphafold_provider`
  now checks Phase 3's own cached API response's `providerId` and excludes
  (`alphafold_model_not_official_afdb`) anything that isn't `"GDM"`, checked
  first, before any mmCIF parsing. This is a real correctness gap in the
  *data*, not in Phase 3's code as originally specified — recorded here
  rather than reopening Phase 3, since Phase 4's gate fully prevents
  contamination of the mapped/benchmarked set either way. See PROGRESS.md
  for the exact count (43/2,704, ~1.6% of this phase's scope — smaller than
  9% of all accessions since most Community-model accessions had already
  been filtered out for other reasons upstream or don't recur across
  multiple representatives).
- **Fallback (PLAN.md item 4):** when SIFTS mapping is absent for a chain
  or fails validation (see below), a global alignment (Biopython
  `PairwiseAligner`, BLOSUM62, free end-gaps via `end_insertion_score`/
  `end_deletion_score = 0` so expression tags/linkers at either terminus
  align to nothing instead of forcing a bad internal alignment) of the
  observed experimental sequence against the AlphaFold model's own
  sequence is tried instead. If the fallback also fails validation, the
  chain is excluded (`fallback_failed_validation`, with both methods'
  compared-count/identity embedded in the reason string).
- **Validation gate (new config constant, justified below):**
  `RESIDUE_MAPPING_VALIDATION_IDENTITY = 0.90` — a mapping method (SIFTS or
  fallback) is trusted for a chain only if ≥90% of its UniProt-mapped,
  AlphaFold-range residues match AlphaFold's residue identity there.
  Deliberately generous to real engineered point mutations (which must
  *not* fail the chain, per this session's instruction — a single mismatch
  in an otherwise-clean chain is recorded as a mutation, not an exclusion)
  while still catching frame-shift/misalignment-scale errors, which in
  practice produce near-zero identity, not a borderline value — confirmed
  empirically on the full run (see PROGRESS.md: real `fallback_failed_validation`
  cases clustered at 55-90% identity, well below both real point-mutation
  rates and the threshold, while every successfully mapped chain landed at
  ≥90%, most at 100%).
- **`MIN_MAPPED_COVERAGE = 0.80`** (new config constant): a chain is
  excluded (`coverage_below_threshold`) if fewer than 80% of its *observed*
  residues end up kept (mapped to UniProt *and* within the AlphaFold
  model's range). Justified by consistency with `CLUSTER_MIN_COVERAGE`
  (also 0.80, Phase 2) — "a meaningful majority of the chain" — while still
  tolerating some loss to expression tags/linkers and residues outside the
  AlphaFold range. The coverage-threshold sweep this session ran
  (`data/interim/phase4_coverage_sweep.csv`) shows this threshold is barely
  binding in practice: the minimum coverage among all 2,604 successfully
  mapped chains was 0.81, and tightening to 0.90 or 0.95 would only drop
  13 or 57 chains respectively (of 2,604) — most chains that map at all
  map almost completely.
- **Outputs:**
  - `data/interim/residue_mappings/{pdb_id}_{chain}.parquet` — one row per
    experimental residue *observed* in the crystal structure (not one row
    per UniProt position; see below): `auth_seq_id`, `auth_ins_code`,
    `label_seq_id`, `residue_name` (observed 3-letter, e.g. `"MSE"`),
    `is_observed` (always `True` in this table — kept for schema
    compatibility with the column recorded above, not used to represent
    unobserved gaps, which are simply absent rows), `uniprot_acc`,
    `uniprot_resnum`, `alphafold_resnum` (`= uniprot_resnum` once
    verified — both kept for schema fidelity with the recorded columns),
    `alphafold_residue_name`, `match_flag` (`None` if not compared),
    `method` (`"sifts"`/`"fallback"`, chain-level, repeated per row),
    `unmapped_reason` (`""` iff kept, including kept-but-mismatched
    engineered mutations — never blank without an explicit reason
    otherwise).
  - `data/interim/mapping_report.csv` — one row per representative in
    scope: `cluster_id`, `pdb_id`, `chain_id`, `uniprot_acc`, the three
    `eligible_*` leakage columns (carried through, no rejoin needed),
    `method`, `n_observed`, `n_kept`, `coverage`, `n_compared`,
    `n_mutations`, `identity`, `status` (`"mapped"`/`"excluded"`),
    `exclusion_reason`.
  - `data/interim/phase4_attrition.csv` — per `LEAKAGE_FILTER_MODES` entry:
    eligible-representative count, mapped/excluded counts, whether
    `n_mapped` clears `MIN_TEST_CHAINS_TARGET`.
  - `data/interim/phase4_coverage_sweep.csv` — reporting only (PLAN.md item
    6): how many already-mapped chains would additionally survive at each
    of `(0.50, 0.60, 0.70, 0.80, 0.90, 0.95)` coverage, without changing
    `MIN_MAPPED_COVERAGE` or which chains Phase 4 actually excluded.
- **src modules:** `src/data/align_residues.py`.
- **Tests:** `tests/test_align_residues.py` (43 tests, no network) — hand-built
  mmCIF fixtures covering insertion codes, alternate locations (occupancy-
  and alt-id-tiebreak), non-sequential/negative author numbering
  (including one fixture asserting a naive constant-offset mapping would
  be wrong where ours isn't), internal gaps, expression tags with no SIFTS
  mapping, residues outside the AlphaFold range, a modified residue (MSE),
  an engineered point mutation (kept, not excluded), a chain requiring the
  fallback alignment, AlphaFold-numbering self-validation
  (success/range-mismatch/wrong-accession/missing-columns), the
  official-vs-Community-provider gate, coverage/attrition/sweep
  arithmetic, and end-to-end `process_representative` on fixture files for
  each terminal outcome (sifts success, fallback success, AF-numbering
  exclusion, coverage exclusion, mutation-chain success).
- **Hand-verification:** `scripts/spot_check_mapping.py` -- seeded random
  sample of N primary-set mapped chains x 3 kept residues each, printed as
  a table plus a ChimeraX command block (opens both structures, selects
  the corresponding residues in each) for visual confirmation. Not part of
  the pipeline itself; reads Phase 4's already-written outputs.
- **Success criterion:** every representative in scope gets a
  `mapping_report.csv` row with a definite status and, on exclusion, a
  specific reason; every parquet row is either kept (`unmapped_reason ==
  ""`) or carries an explicit reason, never a silent drop. Met on the full
  run — see PROGRESS.md for the actual counts, method split, and
  coverage/identity distributions.
- **Geometric validation, added 2026-09-17 (not in the original plan):**
  Phase 4 only ever validated its mapping by *sequence* identity (SIFTS or
  the fallback alignment, gated at 90%), which can't catch a mapping that
  is numerically self-consistent but doesn't correspond in 3D. `scripts/validate_mapping_geometry.py`
  closes this gap for the primary set: Kabsch/SVD-superposes the
  experimental chain's Calpha atoms onto the AlphaFold model's Calpha
  atoms at the same UniProt positions (one rigid-body fit per chain), then
  flags chains where a long contiguous run of residues stays far apart
  after the best-possible global fit. Writes
  `data/interim/phase4_geometric_validation.csv`. See PROGRESS.md for the
  full result and why the ~16% flagged rate turned out to reflect genuine
  AlphaFold-vs-experiment structural disagreement (this project's actual
  subject matter), not a mapping defect — every flagged chain's underlying
  sequence identity was ~100%, mostly via the trusted SIFTS method.

### Phase 5 — Ground-truth interface labels — **DONE** (2026-09-17)

Implemented as a single module, `src/data/interface_labels.py`, per this
session's explicit instruction (superseding the three-module sketch
below, following the same one-module-per-phase precedent as Phases 1-4).

**Deviation from the recorded output path, flagged and resolved in favor
of the session's brief, unlike Phase 4's analogous conflict:** this
session's brief asked for `data/interim/interface_labels/{pdb_id}_{chain}.parquet`
(per-chain, matching the Phase 3/4 convention); the line above (unchanged
since the very first planning commit, `ac34866`) says
`data/processed/interface_labels/{pdb_id}.parquet` (per-entry, with a
`chain_id` column). Unlike Phase 4's Phase-4-vs-this-session path conflict
(where the *recorded* path had already been deliberately corrected in a
previous implementation session), this Phase 5 line was never revised
after being implemented — it is still the original, provisional sketch,
exactly like every other phase's sketch was before its own implementation
session corrected it (see Phases 1-4 above). Treating it as an
authoritative "already-made decision" would be inconsistent with how
every other phase was actually handled. Implemented the session's
explicit per-chain path, matching Phase 3/4's own established convention;
documented below as the actual output.
- **Inputs:** `data/interim/mapping_report.csv` rows with `status ==
  "mapped"` (2,604 chains, Phase 4's full output) plus each chain's own
  `data/interim/residue_mappings/{pdb_id}_{chain}.parquet` (Phase 4) for
  its exact observed-residue key set and UniProt positions;
  `data/raw/pdb/{PDB_ID}_updated.cif.gz` (Phase 3 cache) for geometry.
- **Config used:** `INTERFACE_DISTANCE_CUTOFF_ANGSTROM` (5.0, existing).
  New: `SASA_BURIAL_CUTOFF_ANGSTROM2` (1.0, given directly by this
  session's brief) and `SURFACE_RSA_THRESHOLD` (0.25 — justified below).
- **Assembly choice (this session's item 1): the biological assembly, not
  the asymmetric unit.** For each representative chain, the assembly
  chosen is the (author-preferred, else first-listed) `pdbx_struct_assembly`
  whose generators actually cover that chain (`choose_assembly`); if a
  file defines no assemblies at all, the documented fallback is to treat
  the asymmetric unit itself as "the assembly" (a real, common case — the
  deposited coordinates already *are* the biological unit). The
  assembly's own `pdbx_struct_oper_list` operators are applied via
  `gemmi.make_assembly`, which the AlphaFold mmCIF discovery in Phase 4
  turned out not to need (only the experimental mmCIF has real
  crystallographic symmetry) — see PROGRESS.md for the empirical
  ASU-vs-assembly comparison this unlocks (item 7).
  - **"Which copy is the representative," corrected mid-session on real
    data:** the original design required the representative's copy in the
    assembly to be generated by the *identity* operator specifically
    (verified by coordinates), on the theory that this is the only copy
    guaranteed to correspond to the deposited, Phase-4-mapped coordinates.
    Real data (`7QPB`, chain `A`) falsified the premise: a biologically
    correct trimeric assembly can generate its OWN listed chain only via a
    *non-identity* crystal-symmetry operator (chain `C`'s copy is the
    identity-generated one; chain `A`'s is not) — `make_assembly` still
    preserves `A`'s own `auth_seq_id`/residue identity through that
    transform, so requiring identity specifically wrongly excluded this
    (and similar) chains. Corrected: when a chain appears via exactly one
    copy in the chosen assembly, that copy is unambiguously "the
    representative" regardless of which operator generated it. Only when
    a chain appears via *multiple* copies in one assembly (a true
    multi-copy homomer within that assembly) is a tie-break needed at
    all — documented as: prefer the identity-transformed copy if one
    exists, else the lowest `AddNumber`-assigned copy name, both fully
    deterministic (`find_self_copy`).
- **Partners (item 2):** any other chain in the chosen assembly with at
  least one residue whose chemical-component data
  (`gemmi.find_tabulated_residue(...).is_amino_acid()`) says it's an
  amino acid — checked per residue, not via chain- or entity-level
  metadata (which can be absent/unreliable for non-polymer residues, see
  Phase 4's precedent) — so a partner chain never needs to have passed any
  candidate-selection filter itself, and short peptides count. Ligands,
  ions, water, and nucleic acids are excluded by the same per-residue
  check. A partner is "homomeric" with the representative if any of its
  residues share the representative's own polymer `entity_id`.
- **Primary label:** heavy-atom (`element.is_hydrogen` excluded) distance
  ≤ `INTERFACE_DISTANCE_CUTOFF_ANGSTROM` to any partner heavy atom, via a
  `scipy.spatial.cKDTree` over the pooled partner atoms (one tree per
  chain, queried once per representative atom) — scales to the largest
  assemblies encountered (up to 27 partner chains on the real run).
  Alternate locations are resolved per atom name (highest occupancy,
  ascending altloc on ties) — same rule as Phase 4, generalized from a
  single representative row to the full atom set needed for geometry.
- **Robustness (ΔSASA) label and RSA:** `biotite.structure.sasa`
  (Shrake-Rupley, `vdw_radii="ProtOr"` — freesasa is not installable in
  this container, recorded risk, PLAN.md Sec 6). ΔSASA = isolated-chain
  SASA minus SASA in the assembly context (representative + protein
  partners only, no ligands/waters, for consistency with the distance
  label's partner set), summed per residue; `is_interface_sasa` =
  `ΔSASA > SASA_BURIAL_CUTOFF_ANGSTROM2`. **Correctness note (found before
  the full run, not after):** ProtOr radii are looked up per *(residue
  name, atom name)*, not a bare element symbol — an early draft built the
  SASA `AtomArray` from coordinates alone (element hardcoded to carbon,
  no residue/atom names), which biotite doesn't silently miscompute but
  actively crashes on (`IndexError` inside `vdw_radius_protor`); fixed by
  carrying each atom's real element, atom name, and parent residue name
  through the whole pipeline. RSA = isolated-chain residue SASA / a
  per-amino-acid theoretical maximum (Tien et al. 2013, "Theoretical"
  column — `MAX_SASA_ANGSTROM2`), via the same parent-residue resolution
  as Phase 4 for modified residues (e.g. MSE via its MET parent value);
  `NaN` for non-standard residues with no tabulated maximum.
  `SURFACE_RSA_THRESHOLD = 0.25`: a widely used structural-biology
  rule-of-thumb for "exposed" (residues below this are considered buried
  in the monomer, before any partner is even considered), used to report
  interface fraction among surface residues specifically, not just all
  residues (Phase 7 will need both).
- **Restricted to Phase 4's exact residue set (item 6):** `representative_residues`
  builds the geometry-derived key set (`auth_seq_id`, `ins_code`) directly
  from the chosen assembly copy's own residues (`label_seq is not None`,
  the same "is this part of the polymer sequence" criterion Phase 4 used)
  and asserts it matches Phase 4's parquet key set *exactly* — any
  mismatch (found in one but not the other) excludes the chain
  (`residue_key_mismatch_with_phase4_map`) rather than silently
  reconciling. Confirmed to hold for every successfully-labeled chain on
  the full run.
- **Outputs:**
  - `data/interim/interface_labels/{pdb_id}_{chain}.parquet` — one row per
    Phase-4-observed residue: `auth_seq_id`, `auth_ins_code`,
    `uniprot_resnum`, `is_interface_contact`, `min_partner_distance`,
    `delta_sasa`, `is_interface_sasa`, `rsa`.
  - `data/interim/labels_report.csv` — one row per Phase-4-mapped
    representative: `cluster_id`/`pdb_id`/`chain_id`/`uniprot_acc`, the
    three `eligible_*` leakage columns, `assembly_id`, `n_partner_chains`,
    `has_homomeric_partner`, `n_observed`, `interface_fraction_distance`,
    `interface_fraction_sasa`, `interface_fraction_distance_surface`,
    `n_surface_residues`, `jaccard`, `n_partner_chains_asu_only`,
    `would_exclude_if_asu` (item 7's sanity comparison, computed for every
    chain regardless of outcome), `status`, `exclusion_reason`.
  - `data/interim/phase5_attrition.csv` — per `LEAKAGE_FILTER_MODES` entry.
  - `data/interim/phase5_label_agreement.csv` — pooled (all labeled
    chains) distance-vs-SASA confusion matrix (`tp`/`fp`/`fn`/`tn`) and
    Jaccard, alongside the per-chain `jaccard` column in `labels_report.csv`
    (item 8).
- **Exclusion reasons, logged distinctly:** `no_protein_partner_in_assembly`,
  `zero_interface_residues`, `chain_not_in_any_assembly`,
  `assembly_self_copy_not_found`, `residue_key_mismatch_with_phase4_map`,
  `structure_read_error`, `phase4_mapping_not_available` (Phase 4 parquet
  missing for a nominally-"mapped" row — a data-integrity guard, not
  expected to fire).
- **src modules:** `src/data/interface_labels.py`.
- **Hand-verification:** `scripts/spot_check_interfaces.py` — seeded
  random sample of N primary-set labeled chains, printed as a table plus
  a ChimeraX command block that fetches the chosen assembly *directly from
  RCSB* (`open <pdb_id> assembly <id>`, independent of this project's own
  local assembly-generation code — a genuinely independent check, not a
  visualization of our own derived coordinates) and highlights the
  distance-labeled interface residues.
- **Tests:** `tests/test_interface_labels.py` (18 tests, no network,
  hand-built mmCIF fixtures with real `pdbx_struct_assembly`/
  `pdbx_struct_oper_list` categories) — known contacts just inside/outside
  the cutoff in one fixture; hydrogen atoms, water, and a ligand ignored
  even when geometrically closer than the true nearest protein partner; a
  partner that exists only after applying an assembly symmetry operator
  (plus the ASU-alone sanity check confirming zero partners there); an
  asymmetric-unit crystal contact from a chain excluded from the chosen
  assembly's `asym_id_list`, confirmed not counted even though closer than
  the real (in-assembly) partner; the exact-join assert on both a match
  and two mismatch directions (extra/missing residue); alt-loc occupancy
  tie-break; Jaccard/confusion-matrix arithmetic; `choose_assembly`'s
  no-assemblies-defined and chain-not-covered paths.
- **Success criterion:** every Phase-4-mapped representative gets a
  `labels_report.csv` row with a definite status and, on exclusion, a
  specific reason; every kept residue's labels are computed from exactly
  the Phase 4 residue set, asserted, never assumed. Met on the full run —
  see PROGRESS.md for the actual counts, label-agreement statistics, and
  the ASU-vs-assembly comparison (this phase's most decision-relevant
  finding).
- **Programmatic spatial verification, added 2026-09-17 (not in the
  original plan; requested as a ChimeraX-free alternative to hand
  verification):** `scripts/validate_interface_labels.py` first attempts
  PDBe's PISA/interface-annotation API for an independent external check
  on a seeded sample of primary-set chains; none of the endpoint patterns
  tried returned data this session (see PROGRESS.md), so it falls back to
  two purely-geometric checks on the same sample: (1) contiguity -- most
  labeled interface residues have another interface residue within 8A
  (real interfaces are patches, not scattered noise); (2) centroid -- the
  interface patch's centroid sits closer to the partner chain than a
  same-size sample of comparably-exposed non-interface surface residues'
  centroid does (confirms the patch actually faces the partner). Writes
  `data/interim/phase5_spatial_verification.csv`.

### Scope decision (2026-09-17): lean path for Phases 6-9

This project is taking a **lean path** from here through Phase 9, recorded
here as the current scope. Everything below this note in Phases 6-9
describes the *full* originally-planned design; the lean path is a
deliberate narrowing of scope, not a discovery that the full design was
wrong — the deferred items remain valid future work if there's ever a
reason to go back for them.

**In scope now:**
- Phase 6, run on the **primary set only** (`eligible_homolog`, 566
  Phase-5-labeled chains as of the last full run), with exactly **two**
  input variants per chain: the isolated experimental chain (`exp`) and
  the AlphaFold model trimmed to the experimental chain's mapped UniProt
  range (`af_trimmed`).
- Phases 7 and 8, computed only on those `exp`/`af_trimmed` predictions
  for the primary set.

**Deferred (documented, not run):**
- The `af_full` (untrimmed AlphaFold model) input variant. The code path
  stays working (`run_pesto.py --inputs exp,af_full,af_trimmed` still
  functions) but isn't part of the default/lean invocation.
- The `exact_train` and `none` leakage-filter modes at full size. When
  Phase 7 does report leakage-filter sensitivity, the plan is a **seeded
  random subsample of ~250 chains** from each non-primary mode, not the
  full mode (which would mean re-running Phase 6 inference on thousands
  more chains for a sensitivity check) — a decision that *replaces* the
  full-size sensitivity design recorded in Phase 7 below, not just an
  unimplemented sketch being filled in for the first time.
- The apo/unbound arm (Phase 6 Part D and everywhere it's referenced in
  Phase 7/8/confound (c)).
- Phase 9 (pLDDT-augmented fine-tuning) entirely — consistent with its
  existing gate (only runs if Phase 7 shows a significant gap AND Phase 8
  ties it to pLDDT), which this lean path doesn't attempt to satisfy or
  bypass; Phase 9 simply isn't attempted this pass regardless of what
  Phase 7/8 find.

### Phase 6 — PeSTo inference on isolated single chains — **DONE, lean scope** (2026-09-17)

Implemented as `src/models/run_pesto.py` (orchestrator) +
`src/models/pesto_worker.py` (standalone inference worker), run on the
primary set only with `exp`+`af_trimmed` inputs per the lean-scope
decision above. `af_full` and the apo arm are deferred (code path for
`af_full` still works, just isn't invoked by default).

- **PeSTo setup (Part A):** cloned `external/PeSTo` at commit
  `ba651aa29aaa839d0ee2c458dee90ca1f934d90a` (license: CC BY-NC-SA 4.0).
  Checkpoint: `model/save/i_v4_1_2021-09-07_11-21/model_ckpt.pt` (the one
  `apply_model.ipynb` actually loads -- `model.pt` in the same directory
  is a *different* file, not used), SHA-256
  `6a0f070a90c909f0eb4737ed89050e68c95e4cc59cebd6ae791d41ded163dffe`.
- **Code-reading findings, confirmed not assumed:**
  - Input atoms: all heavy atoms; hydrogens/deuterium and water stripped
    by PeSTo's own `clean_structure`; ligands/ions are *not* auto-removed
    by PeSTo itself (why this project's own input-prep strips them).
  - B-factor/occupancy: not read as model input at all (only element/
    resname/atom-name one-hots) -- confirmed directly in `encode_features`.
  - `i_v4_1` specifically uses *only* the element one-hot as input
    (`encode_features(structure)[0]` in `apply_model.ipynb` -- the
    resname/atom-name one-hots are computed but discarded), matching the
    README's "v4 uses only atom element" note.
  - Per-atom -> per-residue reduction happens *inside* the model via a
    learned attention-pooling layer (`StatePoolLayer`), not an external
    aggregation rule to replicate.
  - **Interface channel: index 0** (`z[:, 0]`), confirmed two ways: (a)
    `config.py`'s `r_types` list order (`protein, nucleic, ion, ligand,
    lipid`), and (b) reproduction below.
  - **Residue-ordering gotcha, verified empirically:** PeSTo's own
    preprocessing (`clean_structure`) renumbers residues sequentially in
    file-encounter order (correctly treating insertion-code changes as
    residue boundaries); a later step (`tag_hetatm_chains` +
    `concatenate_chains`) can physically move a HETATM-flagged residue's
    (e.g. a modified residue) *atoms* elsewhere in the array, but sorting
    by the renumbering recovers true file order regardless -- confirmed
    with a synthetic MSE-containing fixture. Consequence: this project's
    own input files can be written in a fully-controlled residue order
    (sorted ascending by the final key), and PeSTo's output row *i*
    corresponds to that same row *i*, without needing coordinate-based
    matching.
- **Reproduction (Part A.3):** ran a from-scratch wrapper against
  `examples/channel/1BL8.pdb` and compared to PeSTo's own precomputed
  `examples/channel/1BL8_i0.pdb`. Pearson correlation 0.996; max residue
  difference 0.005, exactly bounded by that reference file's 2-decimal-
  place storage precision (`%6.2f`) -- i.e. the difference is fully
  explained by their own file's rounding, not a discrepancy. Self-
  determinism (repeated runs, same environment): bit-identical.
- **A real memory ceiling, found empirically (Part A.4):** this container
  has 7.5GB RAM total, no GPU. PeSTo's `extract_topology` builds a full
  O(n^2) pairwise atom-distance tensor, so peak memory grows worse than
  linearly with atom count. Measured peak RSS: 2338 atoms->1.0GB,
  4855->2.1GB, 6041->3.0GB, 8833->6.2GB; **13822 atoms was OOM-killed
  outright.** `PESTO_MAX_ATOMS = 7500` (new config constant) was derived
  by solving for peak memory <= ~4.8GB (this container's ~6.3GB
  "available" minus a 1.5GB headroom) against this empirical curve. A
  sequential, single-worker projection across *all* Phase-5-labeled
  chains x all 3 input variants (the original, non-lean scope) came to
  **~27 hours** with ~110 `af_full` inputs alone exceeding the memory
  ceiling -- this is what motivated the lean-scope decision recorded
  above, not just wall-clock convenience.
- **One-subprocess-per-job, by design, for two independent reasons:**
  (1) OOM isolation -- an OS-level SIGKILL can only be contained (logged
  by the caller, without losing sibling jobs) if each job is its own
  process; (2) a genuine `sys.modules` collision -- PeSTo's own code does
  absolute imports of a top-level package literally named `src`
  (`from src.model_operations import ...`), which collides with this
  project's own top-level `src` package if both ever resolve in the same
  process. `pesto_worker.py` is deliberately *not* part of the normal
  `src` import graph (always invoked as a bare script, never `-m`/`import`)
  so `sys.path[0]` is always the worker file's own directory, never the
  project root -- confirmed empirically (see PROGRESS.md for the exact
  failure mode this avoids).
- **Concurrency (Part B.3):** `config.PESTO_CONCURRENCY_TIERS`, a fixed,
  conservative schedule -- jobs are sorted by atom count and processed in
  ascending tiers (small: 4 workers x 2 threads; medium: 2 workers x 3
  threads; large, up to `PESTO_MAX_ATOMS`: 1 worker x 4 threads), one tier
  fully before the next, so memory from different tiers is never
  concurrent and each tier's own worst case stays within the same ~4.8GB
  budget the ceiling was derived from.
- **Resumability (Part B.4):** a job is skipped (`is_valid_output`) if its
  final output parquet already exists, is a readable parquet, and has the
  expected column -- checked *before* input preparation, so a resumed run
  does no wasted work at all for already-done jobs.
- **Inputs (Part B.1):** written to `data/interim/pesto_inputs/{pdb_id}_{chain}_{variant}.pdb`
  + a companion `.keys.csv` (same row order as the PDB, header names the
  final join key(s)) so the worker stays completely generic (it never
  needs to know "auth" vs. "uniprot" key semantics -- it just zips its
  model output onto whatever key columns the orchestrator already wrote).
  `exp`: the same alt-loc rule as Phase 4/5 (highest occupancy, ascending
  tie-break), restricted to exactly Phase 4's observed-residue set
  (asserted to match exactly, same philosophy as Phase 5). `af_trimmed`:
  AlphaFold's own numbering (already UniProt-native, Phase 4) sliced to
  `[min, max]` of the experimental chain's *mapped* positions inclusive
  (internal gaps filled, since AlphaFold models are always dense).
- **Per-residue missing heavy atoms (Part B.3, `exp` only):** compared
  against a textbook ideal heavy-atom count per amino acid
  (`IDEAL_HEAVY_ATOM_COUNTS`), parent-resolved for modified residues (same
  convention as Phases 4/5). Stored as `n_missing_heavy_atoms` in the
  `exp` output parquet.
- **Outputs:**
  - `data/interim/pesto_inputs/*.pdb` / `*.keys.csv` — prepared inputs,
    inspectable.
  - `data/processed/predictions/{exp,af_trimmed}/{pdb_id}_{chain}.parquet`
    — `exp`: `auth_seq_id`, `auth_ins_code`, `pesto_interface_prob`,
    `n_missing_heavy_atoms`; `af_trimmed`: `uniprot_resnum`,
    `pesto_interface_prob`. (`af_full` writes the same schema as
    `af_trimmed` if ever invoked.)
  - `data/interim/inference_report.csv` — one row per chain: per-variant
    status/reason/atom-count/residue-count/runtime.
  - `data/interim/phase6_attrition.csv` — per `LEAKAGE_FILTER_MODES`
    entry; since only the primary subset was run this pass, the
    `exact_train`/`none` rows report `n_run`/`n_succeeded` restricted to
    that same primary-set overlap, noted explicitly rather than implying a
    full run.
- **Failure/skip reasons, logged distinctly:** `exceeds_memory_ceiling`,
  `worker_crashed_or_oom:returncode=N` (N negative = killed by signal, e.g.
  -9 for SIGKILL/OOM), `input_prep_error:<detail>`, `no_matching_concurrency_tier`
  (defensive; should not fire given `PESTO_MAX_ATOMS` is itself the last
  tier's bound).
- **Join verification (Part B.2's assert):** after the run, every
  successfully-produced output is checked to cover every one of Phase 4's
  *mapped* residues for that chain (`verify_joins`) -- logged loudly, not
  silently trusted, even though the ordering argument above already makes
  a mismatch unexpected.
- **src modules:** `src/models/run_pesto.py`, `src/models/pesto_worker.py`.
- **Tests:** `tests/test_run_pesto.py` (6 tests, no network, no model
  weights) — `--subset`/`--inputs` selection; an oversized input is
  skipped while the other input for the same chain still runs; an
  existing valid output is skipped on resume without re-preparing its
  input; a simulated worker crash (via a monkeypatched `subprocess.run`)
  is logged and a sibling job in the same batch still completes.
- **Success criterion:** met — see PROGRESS.md for the full run's actual
  counts, runtime, failure breakdown, and the Part C sanity-check AUC.

**Geometric/spatial verification added this session (Phase 4/5
addenda, since ChimeraX-based hand-verification wasn't available):** see
the new bullets under Phase 4 and Phase 5 above
(`scripts/validate_mapping_geometry.py`, `scripts/validate_interface_labels.py`).

### Phase 7 analysis plan (pre-registered 2026-09-18, before any AlphaFold comparison was computed)

Scope this pass: primary set only (566 chains), `exp` vs. `af_trimmed`
(per the lean-scope decision above). Committed to below *before* running
`src/analysis/benchmark.py` against real predictions.

- **Primary endpoint:** per-chain paired difference in AUPR, `exp` minus
  `af_trimmed`, computed only over chains where both are defined (see
  exclusions below). Tested with a two-sided Wilcoxon signed-rank test on
  the paired differences; a 95% bootstrap CI on the **median** paired
  difference (not the mean, since AUPR differences are not expected to be
  symmetric). Significance threshold: **alpha = 0.05**, two-sided — this
  also resolves `PLAN.md` Sec. 5 open question #3 for Phase 7's own use
  (Phase 9's gate is separately out of scope this pass regardless of what
  Phase 7 finds, per the lean-scope decision, so this does not by itself
  trigger Phase 9).
- **Secondary endpoints:**
  - The same paired-difference/Wilcoxon/bootstrap-CI construction, for
    ROC-AUC instead of AUPR.
  - Pooled-residue AUPR and pooled-residue ROC-AUC, computed separately
    for `exp` and for `af_trimmed` (not a paired/difference statistic —
    each input's own pooled number, reported side by side).
- **Bootstrap:** resamples **chains**, not residues (residues within a
  chain are not independent observations) — with replacement, `n =
  config.BOOTSTRAP_N_RESAMPLES` (10,000), seeded
  (`config.PHASE7_BOOTSTRAP_SEED`), 95% CI via the percentile method.
- **Pre-specified strata (descriptive only — reported side by side with
  n per stratum, not additional hypothesis tests, not multiple-comparison
  corrected):**
  1. Homomeric vs. heteromeric partners (`labels_report.csv`'s
     `has_homomeric_partner`).
  2. Interface fraction below vs. at/above a threshold —
     `config.INTERFACE_FRACTION_STRATUM_THRESHOLD` (0.5, new constant),
     applied to `labels_report.csv`'s `interface_fraction_distance`.
  3. Chains flagged by the Phase 4 geometric validation
     (`phase4_geometric_validation.csv`, `longest_far_run >= 5`, the same
     rule already used to report that addendum) vs. unflagged.
  4. All residues vs. surface-only residues (`rsa >=
     config.SURFACE_RSA_THRESHOLD`, existing constant, 0.25).
- **Robustness check:** the primary endpoint (paired AUPR difference,
  Wilcoxon, bootstrap CI) repeated with the ground-truth label swapped
  from the distance label (`is_interface_contact`) to the ΔSASA label
  (`is_interface_sasa`) — same chains, same predictions, only the
  "ground truth" changes.
- **Exclusions, logged with a count, never silent:** a chain is excluded
  from a given metric (ROC-AUC or AUPR, for a given input) if that
  metric is undefined for it (all-one-class labels over the scored
  residue set). A chain is excluded from a *paired* comparison (primary
  or secondary endpoint) if the metric is undefined for *either* input.
  Pooled-residue metrics use every scored residue regardless of any
  single chain's per-chain definedness.
- **Deferred, not part of this pass** (per the lean-scope decision):
  leakage-filter-mode sensitivity (the ~250-chain seeded subsamples for
  `exact_train`/`none`) and any `af_full` comparison.

### Phase 7 — Benchmarking on experimental vs. trimmed-AlphaFold inputs — **DONE, lean scope** (2026-09-18)

Implemented exactly the pre-registered plan above, no deviations. Code:
`src/analysis/benchmark.py` (run as `.venv/bin/python -m src.analysis.benchmark`),
tests in `tests/test_benchmark.py` (5 tests, all passing; full suite 143/143).

**Join / same-residue-set assertion.** For each of the 566 primary-set
chains: `data/interim/interface_labels/{pdb}_{chain}.parquet` (all
Phase-4-observed residues, `uniprot_resnum` null for unmapped ones) is
filtered to mapped residues, then inner-joined to `exp` predictions on
`(auth_seq_id, auth_ins_code)` and to `af_trimmed` predictions on
`uniprot_resnum`. `assert_same_residue_set` independently checks that the
auth-keyed residue set covered by `exp` and the auth-keyed set
corresponding (via the mapping) to `af_trimmed` both equal the full mapped
set, raising rather than silently trusting the join. All 566 chains passed
this assertion with zero exclusions (`results/benchmark/exclusions.csv` is
empty) — Phase 4/5/6's residue bookkeeping held up with no join-level
surprises.

**Results** (full numbers in `results/benchmark/summary.csv` and
`strata.csv`; figures in `results/figures/`):

- **Primary endpoint** — paired AUPR difference (exp − af_trimmed), n=554
  chains (12 excluded, see below): median difference **+0.027**, 95%
  bootstrap CI **[0.019, 0.035]**, two-sided Wilcoxon p ≈ 3.5e-17.
  Experimental inputs are better, consistently and significantly, but by a
  small-to-moderate margin on the AUPR scale, not a dramatic one — see
  "plain-terms" note below.
- **Secondary — paired ROC-AUC difference:** median **+0.027**, CI
  **[0.018, 0.037]**, p ≈ 4.6e-19. Same direction and similar magnitude to
  AUPR.
- **Secondary — pooled-residue metrics** (95,118 residues pooled across
  all 566 chains, large chains dominating by construction): exp AUPR
  **0.739**, af_trimmed AUPR **0.669**, both far above the pooled interface
  base rate of **0.223** (so both are informative, not just "beating a
  coin flip" — but the base rate matters for reading the AUPR numbers at
  face value). Pooled ROC-AUC: exp **0.882**, af_trimmed **0.849**.
- **Prediction-shift descriptives:** per-chain Spearman correlation
  between exp and af_trimmed probabilities and mean absolute probability
  shift are in `per_chain_metrics.csv` (`spearman_exp_vs_af_trimmed`,
  `mean_abs_prob_shift`) — not a hypothesis test, just a description of
  how much the predictions move between inputs.
- **Robustness (ΔSASA labels):** primary endpoint recomputed with
  `is_interface_sasa` instead of `is_interface_contact`: median
  **+0.026**, CI **[0.020, 0.034]**, p ≈ 3.0e-17 — indistinguishable in
  practice from the distance-label result. The primary finding is not an
  artifact of the distance-cutoff label definition.
- **Strata (descriptive, not additional hypothesis tests; all as paired
  AUPR-difference median/CI/p within the stratum)**:

  | stratum | value | n | median diff | 95% CI | p |
  |---|---|---|---|---|---|
  | partner_type | heteromeric | 390 | 0.034 | [0.024, 0.047] | 9.5e-12 |
  | partner_type | homomeric | 164 | 0.017 | [0.007, 0.025] | 3.8e-07 |
  | interface_size | small (<0.5) | 434 | 0.035 | [0.025, 0.046] | 6.9e-12 |
  | interface_size | large (≥0.5) | 120 | 0.015 | [0.009, 0.022] | 3.4e-10 |
  | geometric_validation | unflagged | 462 | 0.024 | [0.015, 0.033] | 6.0e-13 |
  | geometric_validation | flagged | 92 | 0.059 | [0.023, 0.108] | 7.0e-06 |
  | residue_scope | all_residues | 554 | 0.027 | [0.019, 0.035] | 3.5e-17 |
  | residue_scope | surface_only | 554 | 0.020 | [0.014, 0.030] | 8.2e-15 |

  All eight strata show the same direction (exp better) and all are
  individually significant at alpha=0.05 — the effect doesn't reverse or
  vanish in any pre-registered slice, though its size varies (largest for
  the small subset of Phase-4-geometrically-flagged chains, n=92, where
  the CI is also widest; smallest for large-interface and homomeric
  chains). These are descriptive, per the pre-registration, not additional
  confirmatory tests.

**Exclusions.** 12 of 566 chains (2.1%) excluded from every AUPR/ROC-AUC
metric (and so from the paired comparisons): all 12 have every scored
residue labeled interface (single-class ground truth — small chains
entirely at a binding interface, e.g. `n_positive == n` ranging 7-26
residues), making ROC-AUC/AUPR undefined regardless of predictions. Zero
chains excluded for join/residue-set problems. Pooled metrics use all 566
chains' residues regardless (pooling doesn't require per-chain
definedness).

**Plain-terms summary.** PeSTo is measurably and consistently more
accurate on real (experimental) structures than on trimmed AlphaFold
models of the same chains, and this holds up under a stricter label
definition (ΔSASA), across homomeric/heteromeric partners, across
small/large interfaces, across Phase-4-flagged/unflagged chains, and
whether all residues or only surface residues are scored. The effect is
statistically unambiguous (p-values are tiny, driven by a large,
mostly-consistent-direction paired sample, not a borderline result) but
the typical chain-level gap is a few AUPR points (median ≈0.03 on a 0-1
scale) — most chains cluster near the identity line in
`results/figures/paired_aupr_scatter.png`, with a right-skewed tail of
chains where AlphaFold input costs noticeably more. It would overstate
this to call AlphaFold inputs "much worse" across the board; it is more
accurate to say AlphaFold inputs are worse on average, reliably so, and
occasionally by a lot, while for many individual chains the two inputs
are close.

**Anything surprising:** no residue-set/join surprises (0 exclusions there
— Phase 4-6's bookkeeping is solid); the 12 metric-exclusion chains are
all fully-interface small chains, not an edge case anticipated in the
pre-registration but handled cleanly by the "undefined metric" exclusion
rule without any special-casing; the ΔSASA robustness result is
essentially identical to the primary result, which is reassuring (the
finding isn't sensitive to label definition) but also means ΔSASA didn't
add much independent information here; the Phase-4-geometrically-flagged
stratum (n=92) shows the *largest* effect and also the widest CI, worth
keeping in mind if those chains are revisited in Phase 8.

**Deferred, not part of this pass** (per the lean-scope decision):
leakage-filter-mode sensitivity (`exact_train`/`none`), any `af_full`
comparison, and the apo/unbound arm.

### Phase 7 — Benchmarking (original full-scope sketch)

**Superseded for this pass by the pre-registered analysis plan above and
the actual implementation below** — the sketch immediately following was
never revised after Phase 6 changed the actual output paths/schema
(`results/predictions/pesto/*` and `data/processed/interface_labels/*`
don't exist; the real paths are `data/processed/predictions/{exp,af_trimmed}/`
and `data/interim/interface_labels/`), so, per the same precedent as
Phases 3-6, it's treated as the original provisional design, not a
decision to preserve as-is. Kept below for the full (non-lean) design's
historical record.

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

### Phase 8 analysis plan (pre-registered 2026-09-17, before any structural
feature was computed)

Scope this pass: primary set only (the 566 Phase-7-benchmarked chains),
`exp` vs. `af_trimmed` (per the lean-scope decision above; no apo/unbound
arm, no `af_full`). Committed to below *before* running
`src/analysis/structural_metrics.py` or `src/analysis/error_analysis.py`
against real predictions.

- **Primary question:** is the per-residue change in PeSTo predictions
  between `exp` and `af_trimmed` associated with local structural
  divergence and AlphaFold confidence? (Association, not causation — see
  below.)
- **Primary analysis:** per-residue ROC-AUC and AUPR of `af_trimmed`
  predictions within pLDDT bands (`config.PLDDT_BANDS` = 50/70/90, giving
  four bands: `<50`, `50-70`, `70-90`, `>=90`), alongside `exp`
  predictions on the *same* residues (same band membership, defined by
  AlphaFold's own confidence, not by which input is being scored) — this
  isolates whether AlphaFold's own confidence signal predicts *where*
  degradation concentrates. Chain-resampled bootstrap CIs (same
  methodology as Phase 7: resample chains with replacement, not residues,
  percentile method), with a new fixed seed `config.PHASE8_BOOTSTRAP_SEED`.
  **Resample count for this specific CI is `config.PHASE8_BAND_BOOTSTRAP_N_RESAMPLES`
  (2,000), not Phase 7's 10,000** — profiled empirically before running on
  real data: each draw here recomputes a pooled O(n log n) ROC-AUC/AUPR
  over tens of thousands of residues (~0.02s/draw at 50,000 residues),
  unlike Phase 7's bootstrap of a cheap ~500-value per-chain-median array,
  so 10,000 draws across every band x input combination would cost
  20-30 minutes; 2,000 still gives a stable percentile CI at this n. Every
  other Phase 8 bootstrap (length-vs-drop correlation, flagged-chain
  comparisons) resamples cheap per-chain scalars and uses the full
  `config.BOOTSTRAP_N_RESAMPLES` (10,000), same as Phase 7. A band is
  reported only if it has
  `>= config.PLDDT_BAND_MIN_RESIDUES` (30) residues; empty/undersized
  bands are logged, not silently skipped.
- **Secondary analysis:** a regression of the per-residue absolute
  probability shift `|p_af_trimmed - p_exp|` on pLDDT, local RMSD, RSA,
  and secondary-structure class (categorical: helix/strand/coil, coil as
  reference), with **cluster-robust (CRV1) standard errors, clustered by
  chain** (`statsmodels` OLS, `cov_type="cluster"`) — chosen over a mixed
  model for simplicity and because the question here is "are these
  associations there at all, correctly accounting for within-chain
  correlation," not "how much does the effect vary by chain," which is
  what a random-effects/mixed model would additionally buy. Continuous
  predictors (pLDDT, local RMSD, RSA) are z-scored before fitting so their
  coefficients are directly standardized; secondary-structure dummies stay
  0/1 (not standardized in the same sense — noted, not a gap). Variance
  inflation factors (VIF) are reported for every predictor (pLDDT and
  local RMSD are expected to correlate — both reflect the same underlying
  AlphaFold-vs-experimental divergence — so collinearity here is expected,
  not a bug); any predictor with VIF > `config.VIF_COLLINEARITY_THRESHOLD`
  (5) is flagged in the output, not silently included as if uncorrelated.
  Rows missing any regression predictor are dropped with a logged count
  and reason, never silently.
- **Descriptive analyses (associations, side-by-side, not additional
  confirmatory hypothesis tests — same framing as Phase 7's strata):**
  1. Whether the per-chain AUPR drop (`exp - af_trimmed`, from Phase 7's
     `per_chain_metrics.csv`) grows with chain length (`n_observed` from
     `labels_report.csv`) — Spearman correlation with a chain-resampled
     bootstrap CI, a single-predictor OLS slope, and a length-quartile
     table (n, mean/median drop, mean length per quartile) — specifically
     to explain the gap between Phase 7's per-chain **median** drop
     (+0.027) and its **pooled**-residue drop (AUPR 0.739 → 0.669, i.e.
     0.070): pooled metrics implicitly weight by chain length, so if drop
     and length are positively associated, longer chains (which dominate
     the pooled residue count) pull the pooled number above the
     unweighted per-chain median.
  2. What distinguishes the 92 chains Phase 4's geometric validation
     flagged (`longest_far_run >= 5`, same rule as Phase 7's
     `geometric_validation` stratum) from the 462 unflagged chains in the
     same primary set: per-chain mean pLDDT, mean local RMSD, mean global
     Cα distance, chain length, interface fraction, and AUPR drop,
     reported side by side (n, mean, median per group) with a two-sided
     Mann-Whitney U p-value per metric as descriptive context (same
     precedent as Phase 7 reporting p-values on descriptive strata).
- **All results — primary, secondary, and descriptive — are described as
  associations, not as evidence that low pLDDT or high local RMSD
  *causes* worse predictions.** This matters concretely for pLDDT and
  local RMSD, which are expected to correlate (both track the same
  underlying AlphaFold-vs-experimental divergence) — the regression
  reports this via VIF rather than treating either coefficient as an
  isolated causal effect.
- **Local RMSD design (implemented in `structural_metrics.py`):** Cα RMSD
  over mapped residues within `config.LOCAL_RMSD_RADIUS_ANGSTROM` (10 Å,
  per this session's brief) of the residue in question, in the
  *experimental* structure (neighborhood membership is defined by real
  observed geometry, not AlphaFold's), with a **separate rigid-body
  (Kabsch) superposition per residue's own neighborhood**, not one global
  chain-wide superposition. *Why local, not global:* Phase 4's own
  geometric-validation addendum already found that a single global fit
  both (a) can be dragged off by multi-domain hinge motion (two
  individually well-fitting sub-domains that no single global rotation
  can satisfy simultaneously), which would misattribute hinge motion to
  "AlphaFold got this residue wrong," and (b) is exactly the mechanism the
  92-chain geometric-validation flag was built to detect in the first
  place — reusing it as Phase 8's *only* structural-divergence measure
  would conflate "this specific residue's local geometry is off" with
  "this chain has some global misfit somewhere." A per-residue local
  refit isolates genuine local divergence from that confound. Residues
  whose neighborhood has fewer than `config.LOCAL_RMSD_MIN_NEIGHBORS` (3
  — the minimum for a non-degenerate 3D rigid-body fit) mapped residues
  with resolvable Cα coordinates in both structures get `local_rmsd =
  NaN`, logged with a reason, not dropped from the output table.
  **The existing Phase 4 per-residue Cα distance (one global superposition
  per chain, post-fit per-residue distance) is reused as a separate,
  complementary *global* measure** — Phase 4's own script
  (`scripts/validate_mapping_geometry.py`) only persisted chain-level
  summary statistics (median/mean/max distance, `longest_far_run`), not
  the per-residue distances themselves, so `structural_metrics.py`
  recomputes them with the identical method (same Kabsch/SVD superposition
  logic) rather than importing from a one-off, non-package `scripts/`
  file — a deliberate small duplication, not a re-derivation of the
  method.
- **pLDDT source:** AlphaFold DB mmCIFs store per-residue pLDDT as
  `_atom_site.B_iso_or_equiv` (i.e. gemmi's `atom.b_iso`) — confirmed by
  direct inspection this session (not assumed): all atoms within one
  residue carry the identical value (spot-checked), and the values fall
  in AlphaFold's documented 0-100 pLDDT range, not a real crystallographic
  B-factor. `structural_metrics.py` reads the CA atom's `b_iso` per
  residue.
- **Secondary structure source:** `mkdssp`/DSSP is confirmed **not
  installed** in this container (no root access, per the recorded risk) —
  `structural_metrics.py` uses `biotite.structure.annotate_sse` (P-SEA,
  3-state: helix/strand/coil) on the experimental chain's Cα trace. If
  `mkdssp` is ever installed later, this is the one place to switch.
- **RSA source:** reused directly from Phase 5's already-computed
  per-residue `rsa` column (`data/interim/interface_labels/*.parquet`) —
  not recomputed.
- **Deferred, not part of this pass** (per the lean-scope decision):
  the apo/unbound arm (no unbound-structure local RMSD), `af_full`, and
  any leakage-sensitivity-mode analysis (primary/`"homolog"` set only,
  same as Phase 7).

### Phase 8 — Error analysis — **DONE, lean scope** (2026-09-17)

Implemented exactly the pre-registered plan above, no deviations. Code:
`src/analysis/structural_metrics.py` (per-residue feature computation) +
`src/analysis/error_analysis.py` (orchestration, run as
`.venv/bin/python -m src.analysis.error_analysis`). Tests:
`tests/test_structural_metrics.py` (11 tests, synthetic coordinates only)
+ `tests/test_error_analysis.py` (7 tests, synthetic data only); full
suite 161/161 passing. `statsmodels==0.15.0` installed and added to
`requirements.txt` (needed for cluster-robust OLS + VIF).

**Run: all 566 primary-set chains, 95,118 pooled residues, wall-clock
2m48s.** Zero chain-level exclusions (`results/error_analysis/structural_feature_exclusions.csv`
is empty — every chain that had `exp`+`af_trimmed` predictions also had
resolvable Calpha coordinates, pLDDT, and a secondary-structure assignment
for effectively every residue) and zero warnings logged
(`logs/error_analysis.log`) — Phase 4-7's residue bookkeeping held up
again with no new join-level surprises. Secondary structure used
`biotite.structure.annotate_sse` (P-SEA) for all 566 chains — `mkdssp`
confirmed not installed, as expected. pLDDT confirmed read from
`_atom_site.B_iso_or_equiv`.

**Primary analysis — pLDDT-band metrics** (`results/error_analysis/band_metrics.csv`,
figure `results/figures/error_analysis_band_metrics.png`; chain-resampled
95% bootstrap CIs, `n=2,000` draws):

| pLDDT band | n residues | n chains | exp AUPR | af_trimmed AUPR | AUPR gap | exp ROC-AUC | af_trimmed ROC-AUC | ROC-AUC gap |
|---|---|---|---|---|---|---|---|---|
| `<50` | 2,828 | 216 | 0.798 [0.734, 0.858] | 0.543 [0.406, 0.705] | **0.255** | 0.868 [0.836, 0.894] | 0.673 [0.602, 0.738] | **0.195** |
| `50-70` | 5,735 | 445 | 0.766 [0.717, 0.810] | 0.649 [0.587, 0.707] | **0.117** | 0.843 [0.816, 0.867] | 0.768 [0.737, 0.797] | **0.075** |
| `70-90` | 23,762 | 535 | 0.734 [0.699, 0.766] | 0.663 [0.623, 0.701] | **0.071** | 0.853 [0.836, 0.869] | 0.809 [0.790, 0.828] | **0.044** |
| `>=90` | 62,793 | 493 | 0.730 [0.697, 0.761] | 0.690 [0.655, 0.724] | **0.040** | 0.893 [0.882, 0.904] | 0.871 [0.858, 0.885] | **0.022** |

Every band cleared `PLDDT_BAND_MIN_RESIDUES` (30) with wide margin, zero
undefined bootstrap draws in any band. **The AlphaFold-vs-experimental gap
shrinks monotonically as AlphaFold's own confidence increases** — from a
0.26 AUPR gap in the lowest-confidence band down to 0.04 in the
highest-confidence band, a >6x reduction — directly answering the primary
question: yes, AlphaFold's own pLDDT tracks where its structural error
translates into worse PeSTo predictions. (Side observation, not
pre-registered: the interface base rate also drops monotonically with
pLDDT band — 0.40 → 0.33 → 0.26 → 0.19 — consistent with low-confidence
regions disproportionately being flexible loops/termini, which are also
disproportionately interface-labeled; noted, not investigated further.)

**Secondary analysis — shift regression** (`results/error_analysis/regression.csv`;
n=95,117 residues, 1 dropped for missing RSA, logged; 566 clusters;
cluster-robust CRV1 SEs clustered by chain):

| term | coef | cluster-robust SE | p | VIF |
|---|---|---|---|---|
| const | 0.102 | 0.0029 | <1e-270 | — |
| `plddt_z` | **-0.0268** | 0.0033 | 3.9e-16 | 1.64 |
| `local_rmsd_z` | **+0.0268** | 0.0044 | 1.4e-9 | 1.48 |
| `rsa_z` | **+0.0254** | 0.0016 | 1.4e-57 | 1.18 |
| `sse_a` (helix, vs. coil) | +0.0115 | 0.0035 | 9.8e-4 | 1.15 |
| `sse_b` (strand, vs. coil) | +0.0073 | 0.0042 | 0.081 (n.s.) | 1.18 |

All three continuous predictors are significantly associated with the
absolute probability shift in the expected direction (higher pLDDT ->
smaller shift; higher local RMSD and higher RSA -> larger shift), each
controlling for the others. **No VIF exceeds `VIF_COLLINEARITY_THRESHOLD`
(5)** — the largest is 1.64 (`plddt_z`) — so despite pLDDT and local RMSD
being expected to correlate (both track AlphaFold-vs-experimental
divergence), the actual collinearity here is mild; both retain
independent, significant explanatory power, not just one proxying for the
other. Helix residues show a small additional shift over coil; strand
does not reach significance. As pre-registered, these are associations,
not causal claims — pLDDT and local RMSD are two lenses on the same
underlying divergence, and the regression doesn't establish which one, if
either, is doing the "real" work.

**Descriptive — chain length vs. AUPR drop** (`results/error_analysis/length_analysis.csv`,
figure `results/figures/error_analysis_aupr_drop_vs_length.png`; n=554,
same chains as Phase 7's primary endpoint): **no meaningful length
association** — Spearman rho=0.028, 95% CI [-0.054, 0.108], p=0.51; OLS
slope of drop on raw length ≈ -3.3e-6 per residue, p=0.96; length
quartiles (mean lengths 51/99/176/363 residues) show median drops
0.020/0.020/0.041/0.025 — not monotonic. **This contradicts the
pre-registration's a priori hypothesis** (that pooled-vs-median gap was
mainly a length-weighting effect) — stated plainly rather than
reframed after the fact. What actually explains Phase 7's pooled AUPR
drop (0.070, from pooled AUPR 0.739 -> 0.669) being larger than its
per-chain **median** drop (0.027): the per-chain drop distribution is
strongly **right-skewed** (unweighted mean 0.061, roughly 2.2x the
median), and length-weighting the mean barely moves it (length-weighted
mean 0.061, almost identical to the unweighted mean 0.061) — i.e. long
chains are not systematically worse, so weighting by length isn't what's
happening. The remaining gap between even the (weighted or unweighted)
mean of ~0.061 and the pooled value of 0.070 is a reminder that pooled
AUPR is **not a weighted average of per-chain AUPRs at all** — it's a
single rank-based statistic over all residues combined, which can differ
from any linear combination of the per-chain numbers. Net picture: a
minority of chains with unusually large AlphaFold-vs-experimental
degradation pull the mean well above the median, and pooling residues
(rather than averaging per-chain metrics) departs from even that mean —
chain length has nothing to do with it.

**Descriptive — geometrically flagged vs. unflagged chains**
(`results/error_analysis/flagged_chains_summary.csv`, figure
`results/figures/error_analysis_flagged_vs_unflagged.png`; same 92
flagged / 462 unflagged split as Phase 7's `geometric_validation`
stratum, Mann-Whitney U, two-sided):

| metric | flagged mean (median) | unflagged mean (median) | p |
|---|---|---|---|
| mean local RMSD (Å) | 2.21 (1.39) | 0.67 (0.55) | 4.2e-32 |
| mean global Cα distance (Å) | 7.91 (7.37) | 1.40 (0.98) | 9.5e-40 |
| mean pLDDT | 77.4 (84.1) | 87.9 (90.6) | 2.8e-11 |
| chain length (residues) | 224 (202) | 162 (117) | 2.5e-4 |
| per-chain AUPR drop | 0.108 (0.059) | 0.052 (0.024) | 0.02 |
| interface fraction | 0.314 (0.215) | 0.327 (0.266) | 0.18 (n.s.) |

The 92 flagged chains are a coherent, distinguishable population: lower
AlphaFold confidence, much larger local *and* global structural
divergence, longer on average, and a roughly 2x larger PeSTo AUPR drop —
consistent with confound (c)/Phase 4's addendum finding that this flag
catches genuine AlphaFold-vs-experimental 3D divergence (not a mapping
artifact), and consistent with the band-metrics/regression results above.
Interface fraction is the one metric that doesn't differ — being flagged
isn't related to how much of the chain is interface.

**Anything surprising:** the length-vs-drop null result (directly
contradicting the pre-registration's stated hypothesis, see above) is the
most notable finding of this phase — worth remembering if this dataset is
revisited, since "longer chains degrade more" was a reasonable-sounding a
priori guess that the data does not support; the actual explanation
(right-skew + pooled-AUPR nonlinearity) is less intuitive but is what the
data shows. Also notable: pLDDT and local RMSD turned out to be only
mildly collinear (VIF 1.6, not the "expected" strong correlation the
pre-registration flagged as a risk) — both retain independent signal in
the regression, which wasn't guaranteed going in.

**Deferred, not part of this pass** (per the lean-scope decision): the
apo/unbound arm, `af_full`, and any leakage-sensitivity-mode analysis.
**Phase 9 gate:** Phase 7's primary endpoint was significant (p ≈ 3.5e-17)
and Phase 8's band metrics now do tie the gap to pLDDT (monotonic,
>6x shrinkage from lowest to highest confidence band) — **both Phase 9
gate conditions are met** for the first time this project. Per the
lean-scope decision, **Phase 9 is still not attempted this pass
regardless** — that decision was to skip Phase 9 entirely for this lean
pass, not conditional on the gate; if Phase 9 is ever picked back up,
this note is where to look to confirm the gate was already satisfied
2026-09-17, not something to re-derive.

### Phase 8b — pLDDT-filtering deployment check (post-hoc, exploratory — pre-registered 2026-09-18, *after* Phase 7/8's results were already known)

**This is not a pre-registration in the same sense as Phase 7/8's above.**
Phase 8's Discussion argued informally that pLDDT is an "actionable" signal
for flagging unreliable PeSTo predictions on AlphaFold-only input, without
actually testing a filtering rule. This section is written down *before
computing the filtering numbers themselves*, so the specific test isn't
tuned after seeing its own result — but the decision to run this analysis
at all, and its motivating hypothesis, came *after* seeing Phase 7 and
Phase 8's results, which is exactly what makes it exploratory rather than
a true a priori pre-registration. Labeled as such throughout PLAN.md,
PROGRESS.md, and the paper — never presented as if it were fixed before
Phase 7/8 ran.

- **Question:** if a user excludes low-pLDDT residues before trusting
  PeSTo's AlphaFold-input predictions, do precision/recall/F1 on the
  *remaining* residues actually improve, and by how much of the true
  interface does that cost?
- **Primary analysis:** for both inputs (`exp`, `af_trimmed`) and each
  pLDDT cutoff in `config.PLDDT_FILTER_CUTOFFS` (0, 50, 70, 90 — 0 is the
  no-filtering baseline), keep only primary-set residues with
  `plddt >= cutoff` (pLDDT is a property of the AlphaFold model, shared by
  both inputs' rows in `per_residue_features.parquet` — not recomputed per
  input), binarize each input's own continuous probability at a **fixed**
  threshold (`config.PLDDT_FILTER_PREDICTION_THRESHOLD = 0.5`, chosen as
  the natural midpoint, not tuned per input/cutoff to avoid picking each
  curve's own best case), and compute precision, recall, and F1 on the
  kept residues.
- **Reference curve:** the identical sweep run on `exp` input (same
  cutoffs, same residues kept per cutoff, same fixed threshold). Since
  pLDDT is a property of the chain/AlphaFold model, not of which input is
  being scored, running the same cutoff on `exp` isolates whatever
  precision/recall/F1 change is specifically attributable to
  `af_trimmed`'s pLDDT-linked error (Section 3.2) from the generic effect
  of discarding hard-to-predict low-confidence residues (which could
  improve either curve for reasons unrelated to AlphaFold's structural
  error, e.g. low-pLDDT regions being disordered loops that are
  intrinsically harder to classify regardless of input).
- **Cost side, reported alongside every cutoff, never omitted:** the
  fraction of all scored residues discarded, and the fraction of *true*
  interface residues discarded — discarding flexible loops has a real
  cost, since many genuine interfaces sit exactly there (PLAN.md
  Introduction/confound rationale).
- **Bootstrap:** chain-resampled 95% CIs for precision/recall/F1 at each
  (input, cutoff) pair, same method as Phase 7/8
  (`config.PHASE8_BAND_BOOTSTRAP_N_RESAMPLES`, `config.PHASE8_BOOTSTRAP_SEED`
  — reused rather than a third redundant bootstrap-constant pair, since
  this is the same "pooled metric over tens of thousands of residues"
  cost profile as Phase 8's band metrics, not Phase 7's cheap per-chain
  scalar bootstrap). A resample where a cutoff+input combination yields
  zero predicted positives (precision undefined) or zero true positives
  in the kept set (recall undefined) is logged and excluded from that
  metric's CI, not silently treated as 0 or 1.
- **Outputs:** `results/error_analysis/plddt_filtering.csv` (one row per
  input x cutoff: n_residues, frac_residues_kept, n_positives,
  frac_positives_kept, precision/recall/f1 with CIs, n_bootstrap_undefined
  per metric), figure `results/figures/error_analysis_plddt_filtering.png`
  (F1 vs. cutoff for both inputs; fraction of all residues vs. fraction of
  true-interface residues discarded vs. cutoff).
- **src module:** `src/analysis/plddt_filtering.py`, reading
  `results/error_analysis/per_residue_features.parquet` (already has
  `plddt`, `is_interface_contact`, `exp_prob`, `af_trimmed_prob`,
  `pdb_id`, `chain_id` — no new joins needed).
- **Tests:** `tests/test_plddt_filtering.py` (no network, synthetic data):
  precision/recall/F1 computed correctly at a known threshold; cutoff=0
  keeps every residue; a cutoff above every residue's pLDDT keeps zero and
  is reported as such, not silently dropped from the output table;
  bootstrap resamples chains (same discrete-outcome test pattern as
  Phase 8's own bootstrap tests); undefined-metric draws (all-negative or
  all-non-interface resample) are counted, not treated as 0.
- **How this changes the paper:** a new Results subsection (3.3) reports
  the measured numbers in place of the prior unqualified "actionable
  finding" language, and the Discussion is revised to cite them instead of
  asserting the implication. If filtering turns out not to help much at a
  fixed threshold, that is reported plainly, not reframed.

### Phase 8b — pLDDT-filtering deployment check — **DONE, post-hoc/exploratory** (2026-09-18)

Implemented exactly the post-hoc plan above. Code: `src/analysis/plddt_filtering.py`
(run as `.venv/bin/python -m src.analysis.plddt_filtering`), tests in
`tests/test_plddt_filtering.py` (9 tests, synthetic data only); full suite
170/170 passing. Runtime: ~4.5s on the primary set's 95,118 pooled residues.

**Result: naive pLDDT filtering at a fixed 0.5 threshold does not
meaningfully close the experimental-vs-AlphaFold gap, and the finding this
session's brief anticipated ("if filtering turns out not to help much, say
so plainly") is exactly what happened.**

| Input | Cutoff | % residues kept | % true interface kept | Precision (95% CI) | Recall (95% CI) | F1 (95% CI) |
|---|---|---|---|---|---|---|
| Experimental | 0 | 100.0 | 100.0 | 0.730 [0.706, 0.752] | 0.580 [0.549, 0.611] | 0.647 [0.623, 0.669] |
| AlphaFold trimmed | 0 | 100.0 | 100.0 | 0.673 [0.646, 0.700] | 0.535 [0.505, 0.567] | 0.596 [0.571, 0.621] |
| Experimental | 50 | 97.0 | 94.7 | 0.726 [0.703, 0.750] | 0.573 [0.543, 0.605] | 0.640 [0.616, 0.664] |
| AlphaFold trimmed | 50 | 97.0 | 94.7 | 0.687 [0.661, 0.711] | 0.529 [0.500, 0.563] | 0.598 [0.573, 0.624] |
| Experimental | 70 | 91.0 | 85.7 | 0.726 [0.701, 0.750] | 0.566 [0.535, 0.600] | 0.636 [0.610, 0.661] |
| AlphaFold trimmed | 70 | 91.0 | 85.7 | 0.693 [0.667, 0.718] | 0.525 [0.492, 0.559] | 0.598 [0.569, 0.625] |
| Experimental | 90 | 66.0 | 57.0 | 0.737 [0.707, 0.765] | 0.556 [0.517, 0.594] | 0.634 [0.603, 0.662] |
| AlphaFold trimmed | 90 | 66.0 | 57.0 | 0.713 [0.682, 0.743] | 0.523 [0.487, 0.560] | 0.603 [0.573, 0.633] |

Zero undefined bootstrap draws for any (input, cutoff, metric) combination.

**Reading it plainly:** going from no filtering to the strictest cutoff
(pLDDT >= 90) moves AlphaFold-input F1 by only **+0.007** (0.596 -> 0.603),
with heavily overlapping 95% CIs -- not distinguishable from no
improvement -- while discarding **34% of all residues and 43% of the true
interface residues**. The gap to experimental input's F1 narrows only from
0.050 to 0.031, and part of even that narrowing is experimental input's
own F1 drifting down (0.647 -> 0.634) as its harder low-pLDDT residues are
discarded too, not AlphaFold input distinctly "catching up." The
composition of AlphaFold input's small gain is also informative: precision
rises materially (0.673 -> 0.713) while recall is essentially flat
(0.535 -> 0.523) -- filtering trades real coverage of the true interface
for a modest, not-clearly-significant precision gain, not a broad quality
improvement.

**This directly corrects an overreaching claim in the Phase 8 write-up's
Discussion** (audited 2026-09-18, see "Decisions already made" below): the
prior text called pLDDT "the most actionable finding here" and "a
practical... confidence signal for flagging low-reliability regions,"
without ever testing a filtering rule. That language has been replaced
throughout `paper/template.tex` with the measured numbers above. pLDDT
remains a genuine, independent statistical correlate of prediction shift
(Phase 8's regression, Section 3.2) -- that is a real but more modest
finding than a validated deployment filter, and this project no longer
claims the latter.

**Anything surprising:** the size of the disconnect between the pooled
per-band AUPR gap (Phase 8's own headline number, a >6x shrinkage from
lowest to highest pLDDT band) and the near-flat F1-at-fixed-threshold
result here. AUPR/ROC-AUC are threshold-free ranking metrics and pLDDT
banding pools by AlphaFold's own confidence, so the >6x band-level shrink
is real -- but it does not imply that simply cutting off low-confidence
residues at a single fixed operating threshold recovers similarly large
practical gains, since discarding residues costs recall/coverage that a
ranking metric doesn't have to pay for. The two findings are not in
tension; they just answer different questions ("does confidence track
where AUPR/ROC-AUC are worse" vs. "does thresholding on confidence improve
a fixed-threshold classifier's F1"), and conflating them is exactly the
overreach this phase was added to catch.

**Deferred, not attempted here:** a per-input-optimized threshold (rather
than the fixed 0.5), and any non-hard-cutoff use of pLDDT (e.g. a learned
pLDDT-aware model per Phase 9) -- both noted in the paper's Future Work as
open questions this result does not resolve either way.

### Phase 8 — Error analysis (original sketch, superseded)

**Superseded for this pass by the pre-registered analysis plan above and
the actual implementation below** — same precedent as Phases 3-7: this
sketch predates Phase 6/7's real output paths/schemas and was never
revised to match them (e.g. `results/metrics/`-style paths that don't
exist; the real predictions live at
`data/processed/predictions/{exp,af_trimmed}/` and labels at
`data/interim/interface_labels/`), and it specified a logistic regression
of a binary "prediction error" rather than the more informative continuous
absolute-probability-shift regression this session's brief asks for.
Kept below for the record, not as a target to implement literally.

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
floor. **Resolved 2026-09-16 (decision, not just mitigation design):** at
full scale this floor/target concern turned out to be moot — the primary
mode alone clears both `MIN_TEST_CHAINS_TARGET` (100) and
`MIN_TEST_CHAINS_FLOOR` (50) by a wide margin (n=696), so the leakage
threshold was kept at `SEQUENCE_IDENTITY_CUTOFF` (0.30) against the full
train+test+validation union rather than loosened for sample size. See
confound (f) for the tradeoff this choice does *not* resolve (temporal
coverage).

**(f) Primary test set is skewed toward 2021+ releases (limitation, not
mitigated).** Stratifying the primary (`"homolog"`-mode, 0.30 identity)
leakage flag by candidate release year (`data/interim/leakage_by_release_year.csv`,
Phase 2 full-scale run) shows chains released 2018-2020 are almost all
(99-100%) flagged as PeSTo-overlapping — essentially none of them survive
the primary filter — while the flag rate drops to and plateaus around
59-63% from 2021 onward. Practically, almost all of the primary set's 696
survivors are drawn from 2021+ releases; there is essentially no primary-set
coverage of the 2018-2020 window immediately after the AF2 training cutoff
(`PDB_RELEASE_DATE_CUTOFF`, confound (a)). This pattern suggests PeSTo's own
training-data snapshot extends to roughly 2020-2021 (PeSTo was published in
2023, so this is plausible), well past AF2's 2018-04-30 cutoff — i.e. the
two models' effective "knowledge cutoffs" are not aligned, and the primary
test set's chains are systematically farther in time from AF2's cutoff than
from PeSTo's. This is inherent to the recorded leakage definition and is
**not** resolved by the identity-threshold sweep (0.30-0.95 all show the
same year pattern, `data/interim/leakage_threshold_sweep_by_source.csv`) —
report it as a limitation in Phase 10, not something later phases should
try to fix by loosening the filter.

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
   extra corroboration rather than a blocker. *Partial empirical evidence
   found 2026-09-16:* stratifying the full-scale Phase 2 leakage flags by
   candidate release year shows near-total (99-100%) `pesto_homolog_overlap`
   for 2018-2020-released chains, dropping to a ~60% plateau from 2021
   onward — consistent with (not proof of) PeSTo's training snapshot
   extending to roughly 2020-2021, well past the AF2 2018-04-30 cutoff. See
   `data/interim/leakage_by_release_year.csv` and PROGRESS.md.
3. **Partially resolved (2026-09-18):** Phase 7's own primary-endpoint
   significance threshold is now fixed at alpha=0.05, two-sided (see the
   Phase 7 analysis plan pre-registration and the completed Phase 7
   write-up above) — Phase 7 did find a significant primary-endpoint
   result (p ≈ 3.5e-17) under that threshold. Still open: the *combined*
   Phase 9 gate (this plan assumes Phase 7 significance plus a pLDDT-band
   effect in Phase 8) has not been separately confirmed, and per the
   lean-scope decision Phase 9 is out of scope regardless of what Phase 7
   found, so this Phase 7 result does not by itself trigger Phase 9 —
   confirm the combined gate before Phase 8 is run, if Phase 8/9 are ever
   picked back up.

## 6. Risks

- **MMseqs2** — **resolved (2026-09-16):** the official static AVX2 Linux
  build installed cleanly (no root required) at `external/mmseqs/bin/mmseqs`
  (commit `d401e78c2d18a822cdb1527d7464a043f6035a15`); the Biopython-fallback
  plan was not needed.
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
