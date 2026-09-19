# PPI Interface Transfer

**Do machine-learning models that find protein binding sites still work when
the input structure comes from AlphaFold instead of a real experiment?**

This repository benchmarks a protein-protein interface predictor called
PeSTo on AlphaFold-predicted protein structures, compares it against the
same predictor run on experimentally determined structures, and tries to
explain any gap it finds.

## 1. What this project is about

Proteins rarely act alone. Most biological processes — from a cell
responding to a hormone to your immune system recognizing a virus — depend
on two or more proteins physically touching each other in a specific way.
The patch of surface where two proteins touch is called a **binding
interface**, and the specific amino acids (the residues that make up a
protein's sequence) that sit on that patch are called **interface
residues**. Knowing which residues form an interface matters for
understanding disease mutations, designing drugs, and mapping how proteins
work together across a whole cell.

Figuring out a protein's interface experimentally — by solving its 3D
structure with a bound partner — is slow and expensive. In the last few
years, machine learning models have gotten good at *predicting* interfaces
directly from a protein's 3D shape, without needing a partner protein
present at all. One such model, called **PeSTo**, was trained entirely on
experimentally determined structures. This project asks a simple but
previously under-tested question: **if you feed PeSTo a computationally
predicted structure instead of a real one, does it still work?**

## 2. Background

**AlphaFold** is a deep-learning system that predicts a protein's 3D
structure from its amino-acid sequence alone, without needing a laboratory
experiment. It was a major breakthrough because experimental structure
determination (X-ray crystallography, cryo-electron microscopy) is slow: the
Protein Data Bank (**PDB**), the main repository of experimentally solved
structures, contains roughly 220,000 entries after decades of work. The
**AlphaFold Protein Structure Database** (AlphaFold DB), by contrast,
now hosts predicted structures for over 200 million protein sequences —
essentially every protein science knows the sequence of. That gap is why
predicted structures increasingly stand in for experimental ones in
downstream analysis: there just aren't enough experimental structures to go
around.

But a prediction isn't a fact. AlphaFold reports its own uncertainty for
every residue it models, called **pLDDT** (predicted Local Distance
Difference Test), a 0–100 confidence score. Scores above 90 mean "very high
confidence," 70–90 "confident," 50–70 "low confidence," and below 50 often
means the region has no fixed structure at all (it's intrinsically
disordered). Confidence isn't evenly distributed: flexible loops and the
ends of protein chains — exactly the regions that often form binding
interfaces — tend to get the lowest scores.

**PeSTo** (Protein Structure Transformer) is a geometric deep-learning model
that reads in a protein's atoms (just their element type and 3D position, no
hand-crafted chemistry features) and predicts, per-residue, whether that
residue sits on a binding interface. It was trained exclusively on
experimentally determined PDB structures. A **GNN** (graph neural network) —
the broader family PeSTo belongs to — is a model that represents its input
(here, a protein) as a graph of connected nodes (atoms or residues) and
learns patterns from how they're spatially arranged.

## 3. The central question

It's tempting to assume a model trained on experimental structures will
"just work" on predicted ones, since AlphaFold structures often look very
close to the real thing at a glance. But PeSTo was never shown a predicted
structure during training, and the residues most relevant to interface
prediction — flexible surface loops — are exactly where AlphaFold's
confidence tends to be lowest and its geometry least reliable. Whether
PeSTo's accuracy survives that shift is an empirical question, not a given,
and it hadn't been directly tested for AlphaFold structures before this
project (a related study benchmarked PeSTo on ESMFold, a different
structure predictor — see [Deviations](#deviations-from-the-original-proposal)).

Concretely, this project asks: **how much does PeSTo's interface-prediction
accuracy drop when it's given an AlphaFold structure instead of an
experimental one for the same protein, and can that drop be explained by
AlphaFold's own confidence scores?** The answer turned out to be yes, and
tied to low pLDDT (Phase 8) — a natural follow-up is whether giving PeSTo
access to pLDDT as an input, or simply filtering on it, can recover some of
the lost accuracy. This project tested the simpler filtering version
directly (Phase 8b) and found it doesn't help much; training PeSTo on
pLDDT as an input feature (Phase 9) remains future work, not yet attempted.

## 4. Approach

The full phase-by-phase specification, including every threshold and known
caveat, lives in [`PLAN.md`](PLAN.md); this section is a plain-language
summary of *why* each phase is designed the way it is.

1. **Candidate selection.** Query the PDB for recent, high-resolution X-ray
   structures of protein complexes. "Recent" specifically means *released
   after 2018-04-30* — the date AlphaFold's own training data was cut off —
   so that AlphaFold is genuinely predicting these structures rather than
   having memorized them during its own training.
2. **Redundancy reduction and leakage control.** Many PDB entries are
   near-duplicates of each other (e.g. the same protein pair crystallized
   dozens of times with minor tweaks), so candidates are clustered by
   sequence similarity and only one representative per cluster is kept.
   Separately, since PeSTo itself was trained on "all biological assemblies
   from the PDB," candidates are checked for sequence-level similarity
   (not just identical PDB IDs — a protein solved after 2018 can still be a
   near-duplicate of something PeSTo trained on) against PeSTo's own
   published training data, so the benchmark doesn't unintentionally test
   PeSTo on data it already knows the answer to.
3. **Structure download.** Fetch both the experimental (PDB) and predicted
   (AlphaFold DB) structure for every surviving candidate protein — the
   matched pair the whole comparison depends on.
4. **Residue-numbering mapping.** Experimental structures and AlphaFold
   models don't always number residues the same way (missing loops, tags,
   gaps), so residue numbers are mapped onto a common reference before
   anything is compared, or a comparison would silently misalign.
5. **Ground-truth interface labels.** For each experimental structure, the
   true interface residues are computed directly from 3D coordinates (which
   residues from different chains sit close together).
6. **Run PeSTo.** The pretrained PeSTo model is run on each protein in
   isolation (as if its partner weren't present, matching how it would be
   used in practice) — once on the experimental structure, once on the
   AlphaFold structure, and, where available, once on an independently
   solved "unbound" (partner-free) experimental structure. The unbound arm
   exists to help separate *AlphaFold's* structural error from ordinary
   *biology*: a protein's shape can genuinely change a little when it binds
   a partner, and that's not a modeling failure.
7. **Benchmarking.** Compare PeSTo's predictions against the true interface
   labels for both the experimental and AlphaFold inputs, on the same
   proteins, and measure the accuracy gap.
8. **Error analysis.** Break the accuracy gap down by AlphaFold's pLDDT
   confidence score and by how far the AlphaFold structure's local geometry
   deviates from the experimental one, to see whether low-confidence,
   geometrically-off regions really are where PeSTo fails.
9. **pLDDT-aware fine-tuning (conditional).** If — and only if — steps 7–8
   show a real, confidence-linked accuracy gap, PeSTo is given pLDDT as an
   extra input and lightly retrained to see whether that recovers some of
   the lost accuracy.
10. **Report.** Regenerate all figures and a final write-up directly from
    the saved results, so the findings are reproducible from data alone.

## 5. Current status

*Last updated: 2026-09-18. See [`PROGRESS.md`](PROGRESS.md) for the
detailed, actively-maintained log this section is drawn from.*

- ✅ **Phase 1 — Candidate selection**: done, run at full scale against the
  entire filtered search (no subset cap): 22,887 candidate protein chains
  across 8,808 PDB entries, spanning 4,029 unique UniProt accessions. (An
  earlier 500-entry pilot run is kept on disk for comparison but is no
  longer the active dataset.)
- ✅ **Phase 2 — Redundancy reduction + leakage flagging**: done, also run
  at full scale. Sequence clustering collapsed the 22,887 candidate chains
  down to 3,009 non-redundant representatives. Checking those
  representatives against PeSTo's own published training data found that
  roughly three-quarters are sequence-similar enough to count as
  "overlapping" under the project's leakage definition. **Leakage threshold
  decided (2026-09-16):** the primary benchmark test set is the 696
  representatives with <30% sequence identity to PeSTo's combined
  train+test+validation data; sensitivity analyses at looser thresholds
  remain available. One caveat worth flagging: under this definition,
  almost none of the primary set's chains were released in 2018-2020 (right
  after AlphaFold2's own training cutoff) — they're nearly all flagged as
  overlapping PeSTo's data — so the primary set skews toward 2021-and-later
  releases, which are farther in time from AlphaFold2's cutoff than the
  minimum required by this project's filter.
- ✅ **Phase 3 — Download structures**: done, full scale, all 3,009
  representatives. Every one of the 2,294 unique PDB entries downloaded
  successfully. AlphaFold models were obtained for 2,704/3,009
  representatives (90%) — the rest failed for one of three distinct,
  logged reasons: the UniProt accession isn't in AlphaFold DB (222), the
  accession's sequence was split across multiple AlphaFold fragment models
  rather than one (61), or the AlphaFold model's own sequence no longer
  matches the accession's current UniProt sequence, so it was skipped
  rather than risk mismatched residue numbering later (22). For the
  primary (non-overlapping) benchmark set specifically, 609 of 696
  representatives now have both structures ready to go — comfortably above
  the ≥100 target. Total download: ~900 MB, ~71 minutes.
- ✅ **Phase 4 — Residue-numbering mapping**: done, full scale, all 2,704
  representatives with both structures downloaded. This is the step that
  turns "we have two structure files" into "we know, residue by residue,
  which atom in the experimental structure corresponds to which position
  in the AlphaFold model" — done via each PDB entry's own embedded
  UniProt cross-references (not any assumed numbering offset, which real
  entries in this dataset would get wrong: one representative's mapping
  shifts author residue 1335 to UniProt position 1137). 2,604/2,704
  chains (96%) mapped successfully; 100 were excluded, each for one of
  three specific, logged reasons — most notably, this step caught that
  ~9% of the "AlphaFold" models Phase 3 downloaded were actually
  third-party community submissions built with a different tool
  (ColabFold), not genuine AlphaFold2, and excluded those rather than
  silently treating them as equivalent. For the primary benchmark set,
  579 of 609 representatives (95%) now have a validated mapping —
  579/696 (83%) of the original leakage-filtered primary set has survived
  every phase run so far, still comfortably above the ≥100 target.
- ✅ **Phase 5 — Ground-truth interface labels**: done, full scale, all
  2,604 Phase-4-mapped representatives. For each chain, this step finds
  its real partner protein chain(s) in the crystal and labels which of
  its residues sit at that interface — using the true biological
  assembly (the depositor's annotated functional unit) rather than the
  raw crystal contents, since crystals routinely pack unrelated copies
  of a protein right next to each other with no biological meaning. That
  distinction mattered enormously in practice: only 1 of 2,604 chains
  would appear to have zero partners if the raw crystal contents were
  used directly, but 75 chains are correctly recognized as having no real
  biological partner once the proper functional assembly is used instead
  — i.e. crystal packing would have manufactured a fake interface for
  nearly every genuinely single-chain protein in the dataset had this
  distinction not been made. 2,527/2,604 chains (97%) were labeled
  successfully; two independent labeling methods (a distance-based
  measure and a buried-surface-area-based measure) agree on 90% of
  residues pooled across all chains. For the primary benchmark set, 566 of
  579 representatives (98%) now have interface labels — 566/696 (81%) of
  the original leakage-filtered primary set has survived every phase run
  so far, still comfortably above the ≥100 target.
- ✅ **Phase 6 — PeSTo inference**: done, on a deliberately narrowed
  ("lean") scope decided this session: the primary (non-overlapping)
  benchmark set only, comparing PeSTo's predictions on the real
  experimental structure against its predictions on the AlphaFold model
  trimmed to the same covered region (a third variant, the untrimmed
  full-length AlphaFold model, and the apo/unbound arm, are deferred as
  documented future work, not run). All 566 of 566 primary-set chains
  succeeded on both inputs — no failures, and none exceeded this
  container's memory limits, which a smaller pilot run showed would
  otherwise be a real risk for the largest proteins. Before trusting any
  of this, a sanity check compared PeSTo's predictions on real structures
  against the measured ground-truth interfaces from Phase 5: they agree
  strongly (a pooled discrimination score of 0.88 out of a possible 1.0,
  where 0.5 would mean no better than random) — confirming the whole
  pipeline is wired correctly before any real benchmarking begins.
  Two additional automated checks (substituting for visual inspection)
  were also run this session: the AlphaFold-vs-experimental structural
  agreement for the primary set, and the spatial plausibility of Phase
  5's interface labels — both described further in `PROGRESS.md`.
- ✅ **Phase 7 — Benchmarking (experimental vs. trimmed AlphaFold)**: done,
  same lean scope as Phase 6 (566 primary-set chains, real structure vs.
  trimmed AlphaFold model). Before looking at any result, the exact
  statistical tests, significance threshold, and descriptive breakdowns
  were written down in `PLAN.md` and committed to in advance, so the
  comparison couldn't be tuned after the fact. **Headline finding:**
  PeSTo is more accurate on real structures than on AlphaFold models of
  the same chains — reliably so (this held in every one of eight
  different ways the data was sliced, and also under a stricter,
  alternative definition of "interface"), but the typical difference per
  chain is modest, not dramatic: on a 0-to-1 accuracy scale where higher
  is better, the median chain-level gap is about 0.03. In plain terms —
  AlphaFold inputs are worse on average and the pattern is very unlikely
  to be due to chance, but for many individual chains the two inputs give
  similar results, and it would overstate the finding to call AlphaFold
  inputs "much worse" across the board. 12 of the 566 chains (2%) couldn't
  be scored at all because every one of their measured residues was
  labeled part of the interface, leaving no non-interface residues to
  compare against (a technical requirement of the accuracy measures used,
  not a data problem). Full numbers, figures, and the strata breakdown are
  in `results/benchmark/` and `results/figures/`, and described further in
  `PROGRESS.md` and `PLAN.md`.
- ✅ **Phase 8 — Error analysis**: done, same lean scope (566 primary-set
  chains). This phase asks *why* AlphaFold inputs are less accurate: is it
  associated with AlphaFold's own confidence score (pLDDT) and with how
  much the AlphaFold model's 3D structure actually differs from the real
  one nearby each residue? As before, the exact analyses were written down
  and committed to in `PLAN.md` before any of them were run.
  **Headline finding:** yes — the accuracy gap between real and AlphaFold
  structures shrinks steadily as AlphaFold's own confidence increases,
  from a large gap in AlphaFold's least-confident regions to a small one
  in its most-confident regions (more than a 6x reduction end to end). A
  statistical model confirms that both AlphaFold's confidence score and
  local structural divergence are each independently linked to how much a
  prediction shifts between the two inputs, after accounting for the
  other and for how exposed the residue is. All results here describe
  associations, not proven causes. One pre-registered guess turned out to
  be wrong, and is reported as such rather than quietly dropped: whether
  longer protein chains show a bigger accuracy drop — they don't
  (essentially no relationship); the reason Phase 7's "typical per-chain"
  number (0.03) looked smaller than its "everything pooled together"
  number (0.07) turned out to be that a minority of chains have an
  unusually large drop and skew the average upward, not that long chains
  are systematically worse. A side comparison also confirmed that the 92
  chains flagged earlier (Phase 4) for unusual 3D disagreement with
  AlphaFold really are a distinct group: lower AlphaFold confidence, much
  larger structural divergence, and a bigger accuracy drop than the rest.
  Full numbers and figures are in `results/error_analysis/` and
  `results/figures/`, and described further in `PROGRESS.md` and
  `PLAN.md`.
- ✅ **Phase 8b — pLDDT-filtering deployment check (post-hoc, exploratory)**:
  done. Phase 8's write-up had informally called pLDDT an "actionable"
  signal for flagging unreliable AlphaFold-input predictions without ever
  testing that directly — this phase tests it. **Headline finding: it
  doesn't help much.** Discarding low-pLDDT residues at a fixed prediction
  threshold barely moves AlphaFold-input accuracy (a ~0.007-point gain in
  F1 from no filtering to the strictest cutoff, not clearly distinguishable
  from no improvement given the uncertainty), while throwing away 43% of
  the true interface residues at that same strictest cutoff. Reported
  plainly as a limited result, not reframed as a success — this replaced
  the earlier, untested "actionable" language in the paper. Explicitly
  marked as post-hoc/exploratory (the idea only came up after seeing Phase
  7/8's results), not a pre-registered confirmatory test. Full numbers in
  `results/error_analysis/plddt_filtering.csv`.
- ⬜ Phase 9 (pLDDT-augmented fine-tuning) — both of its pre-registered
  gate conditions are now satisfied (Phase 7 found a significant gap, and
  Phase 8 tied it to AlphaFold's confidence score), but this project is
  deliberately not attempting Phase 9 this pass regardless, per the lean
  scope decided at Phase 6 — it remains documented future work, not
  something currently planned. Phase 8b's result doesn't resolve whether a
  *learned* pLDDT-aware model (Phase 9) would fare better than the naive
  fixed-threshold filter it actually tested. Phase 10 (report generation)
  has not been formally started, though `paper/results_paper.pdf` (see
  below) now covers most of what it would have produced.
- ✅ **Results paper** (`paper/results_paper.pdf`, built by
  `paper/build_paper.py` — see "Building the results paper" below):
  written this session, then adversarially audited against the underlying
  data in the same session. The audit found and fixed a terminology
  inconsistency (the 696-representative leakage-filtered set and the
  566-chain final analyzed set were both being called "the primary
  benchmark set"), a real mischaracterization of a Phase 1 filter
  ("modeled residues" — it's actually the polymer entity's full sequence
  length, which differs from the number of residues actually resolved in
  the structure for 86% of candidates), an unexplained 1-residue
  discrepancy between two pooled-residue counts (traced to a single
  non-standard residue with an undefined solvent-accessibility reference
  value, now logged), and two bugs in the documented reproduction
  commands (see "Reproducing the environment" below). See `PROGRESS.md`
  for the full list.

## 6. Repository layout

```
src/            importable Python package (data pipeline, models, analysis)
data/raw/       downloaded data, cached exactly as fetched (not committed)
data/interim/   intermediate pipeline outputs (not committed)
data/processed/ final, analysis-ready datasets (not committed)
results/        figures, tables, metrics (not committed)
external/       third-party tools (e.g. the MMseqs2 binary, PeSTo checkout)
notebooks/      exploratory notebooks
tests/          pytest test suite
logs/           run logs
docs/           the original project proposal
```

Large or regenerable data, results, and logs are intentionally left out of
version control (see `.gitignore`); only the code that produces them is
committed, so everything is reproducible by rerunning the pipeline.

### Reproducing the environment

```bash
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest
```

All file paths and tunable thresholds (resolution cutoffs, date cutoffs,
sequence-identity cutoffs, etc.) are centralized in `src/config.py` — see
`CLAUDE.md` for the full set of repository conventions.

### Running the full pipeline from a clean checkout

The sequence below is the exact, complete set of commands that reproduces
every number in this repository and in `paper/results_paper.pdf`, run in
order, on a fresh checkout. Nothing here depends on a file that exists only
because of some other, undocumented one-off command — every input any step
reads was itself produced by an earlier step in this same list (verified
2026-09-18; see PROGRESS.md for the audit). Each entry notes whether it
needs network access and its approximate wall-clock time, measured on this
project's container (10 CPUs, no GPU, 7.5GB RAM). **Steps marked NETWORK
cache everything they fetch under `data/raw/`** (CLAUDE.md convention), so
re-running them after an interruption doesn't re-fetch anything already on
disk, and only the first run pays the full network cost.

```bash
# Phase 1: candidate complex selection -- NETWORK (RCSB Search/GraphQL APIs)
# ~100s for the full, uncapped search. Omit --max-entries entirely: that
# flag exists only for a fast (~10s) smoke test on a small, non-representative
# subset -- the real analysis (22,887 candidate chains, 8,808 entries as of
# the run that produced this project's numbers) used no cap. Note that PDB
# grows over time, so an uncapped rerun today will not byte-for-byte match
# the committed candidates.csv from a run on an earlier date.
# (writes data/interim/candidates.csv)
.venv/bin/python -m src.data.select_complexes

# Phase 2: redundancy reduction + PeSTo leakage flagging -- NETWORK the
# first time only (downloads PeSTo's published split files and the bulk
# PDB pdb_seqres.txt.gz; both cached after that), then local MMseqs2
# clustering + homology search. Not separately timed this session; local
# compute only once its inputs are cached, comparable in scale to Phase 4
# (minutes, not hours) -- see logs/cluster_and_split.log for a full trace.
# (writes data/interim/candidates_dedup.csv, clusters.tsv, leakage_threshold_sweep.csv)
.venv/bin/python -m src.data.cluster_and_split

# Phase 3: download PDB (updated mmCIF) + AlphaFold DB structures -- NETWORK,
# ~71 minutes for the full 3,009-representative download. Add
# --estimate-only first to preview the time/disk projection without
# downloading anything.
# (writes data/interim/fetch_report.csv, phase3_attrition.csv)
.venv/bin/python -m src.data.fetch_structures

# Phase 4: residue-numbering mapping -- LOCAL ONLY, ~2.5 minutes
# (writes data/interim/mapping_report.csv, residue_mappings/*.parquet)
.venv/bin/python -m src.data.align_residues

# Hand-verification (optional, not read by any later step): prints a table +
# ChimeraX commands for N random mapped primary-set chains
.venv/bin/python scripts/spot_check_mapping.py --n 3 --seed 0

# Phase 5: ground-truth interface labels -- LOCAL ONLY, ~4.5 minutes
# (writes data/interim/labels_report.csv, interface_labels/*.parquet)
.venv/bin/python -m src.data.interface_labels

# Hand-verification (optional, not read by any later step): prints a table +
# ChimeraX commands (fetching the chosen biological assembly directly from
# RCSB) for N random labeled primary-set chains
.venv/bin/python scripts/spot_check_interfaces.py --n 3 --seed 0

# Phase 6: PeSTo inference (primary set, exp + af_trimmed; lean scope) --
# LOCAL ONLY (one-time NETWORK git clone of external/PeSTo the first time
# only), ~65 minutes for the full 566-chain/1,132-job primary-set run
# (writes data/processed/predictions/{exp,af_trimmed}/*.parquet)
.venv/bin/python -m src.models.run_pesto --subset primary --inputs exp,af_trimmed

# Sanity check (optional, not read by any later step): pooled ROC-AUC of
# exp predictions vs. Phase 5 labels -- worth checking before trusting
# Phase 7/8, but nothing downstream reads its output
.venv/bin/python scripts/phase6_sanity_check.py

# Geometric verification -- LOCAL ONLY, a few seconds. NOT OPTIONAL despite
# living in scripts/: src/analysis/benchmark.py reads
# data/interim/phase4_geometric_validation.csv directly for the
# "geometrically flagged" stratum, and src/analysis/error_analysis.py
# depends on it transitively through benchmark.py's load_chain_metadata().
# Phase 7 will still run without it, but silently drops that stratum.
.venv/bin/python scripts/validate_mapping_geometry.py

# Spatial verification (optional, not read by any later step) -- best-effort
# PDBe PISA network calls with a local geometric fallback; safe to skip
.venv/bin/python scripts/validate_interface_labels.py --n 20 --seed 0

# Phase 7: benchmarking PeSTo on exp vs. af_trimmed -- LOCAL ONLY, ~10s
# (writes results/benchmark/{per_chain_metrics,summary,strata,exclusions}.csv,
# results/figures/*.png)
.venv/bin/python -m src.analysis.benchmark

# Phase 8: error analysis -- pLDDT bands, local RMSD, RSA, secondary
# structure vs. prediction shift -- LOCAL ONLY, ~3 minutes
# (writes results/error_analysis/*.csv + *.parquet, results/figures/error_analysis_*.png)
.venv/bin/python -m src.analysis.error_analysis

# Phase 8b: pLDDT-filtering deployment check (post-hoc, exploratory -- see
# PLAN.md; not pre-registered before Phase 7/8's own results were seen) --
# LOCAL ONLY, ~5 seconds
# (writes results/error_analysis/plddt_filtering.csv,
# results/figures/error_analysis_plddt_filtering.png)
.venv/bin/python -m src.analysis.plddt_filtering

# Results paper -- LOCAL ONLY, well under a minute (matplotlib figure
# generation + two pdflatex passes)
# (writes paper/results_paper.pdf, paper/figures/*)
.venv/bin/python paper/build_paper.py
```

### Rebuilding just the paper

Once `results/` and `data/interim/` are populated (i.e. after the sequence
above, or from an existing checkout that already has them), the
print summary (9 pages as of the 2026-09-18 audit revision) written for a general-science-background reader can
be rebuilt on its own, without rerunning any phase:

```bash
.venv/bin/python paper/build_paper.py
```

This computes every reported statistic directly from the pipeline's own
output files (never hand-typed), regenerates all 7 figures as print-ready
vector PDFs (colorblind-safe Okabe-Ito palette; the one exception is the
example-structure figure, a raster PNG -- see below) into `paper/figures/`,
fills `paper/template.tex`'s placeholders, and compiles the PDF.
**PDF route:** `pdflatex` (found on `PATH`; run twice for cross-references),
checked for before falling back to a pip-installable route, per the build
script's own preference order. **Structure figure:** Figure 6 is a real
molecular render (PyMOL, cartoon representation) produced by
`scripts/render_example_chain.py` -- see below for how it was made and how
to reproduce or repoint it. If `paper/figures/example_chain_render.png` is
ever absent (e.g. a fresh checkout before that script has been run),
`fig_example_structure` in `build_paper.py` falls back automatically to a
plain matplotlib 3D rendering of the Calpha trace, colored the same way but
not a rendered cartoon/surface -- no code changes needed either way.

#### The structure figure: how it was rendered, and how to reproduce it

**Tool: PyMOL (`pymol-open-source` 3.1, via conda-forge), in its own isolated
conda environment at `external/pymol_render_env` (gitignored, never `.venv`
or the `astro` conda env).** Two other options were
tried first and ruled out for this container specifically (not as a general
judgment on either tool):

- **UCSF ChimeraX** (free for academic/non-commercial use): the Ubuntu 24.04
  `.deb` unpacks fine with `dpkg-deb -x` (no root needed), but its `ChimeraX`
  binary requires glibc >=2.34/2.38, newer than this container's Ubuntu
  20.04 (glibc 2.31). The last ChimeraX release with an Ubuntu 20.04 build
  (1.9, glibc-compatible) runs, but its offscreen rendering needs
  `libGL.so.1`, which is not present system-wide and has no root-installable
  path here -- resolving Mesa's own further transitive dependencies
  (`libglx-mesa0`, a software (llvmpipe) DRI driver, a matching `libLLVM`,
  `libdrm2`, `libsensors5`, ...) by hand-downloading individual `.deb`s
  without `apt` was judged not worth the remaining time budget once PyMOL
  (below) was already confirmed working.
- **py3Dmol/NGL via headless Chromium (playwright)**: not attempted, since
  PyMOL succeeded first.

**Install** (isolated, ~1.6 GB, does not touch `.venv` or `requirements.txt`):

```bash
conda create -y -p external/pymol_render_env -c conda-forge pymol-open-source
```

**Produce the figure** (regenerates `paper/figures/example_chain_render.png`
for a given chain; two-stage internally -- data prep runs under `.venv`,
rendering runs under the isolated PyMOL env above, invoked as a subprocess
-- but it's one command either way):

```bash
.venv/bin/python scripts/render_example_chain.py \
    --pdb-id 7DZ9 --chain-id C --output paper/figures/example_chain_render.png
```

To point it at a different chain: any `{pdb_id}_{chain_id}` with Phase
4/5/6 outputs on disk (`data/interim/residue_mappings/`,
`data/interim/interface_labels/`, `data/processed/predictions/{exp,af_trimmed}/`)
works. The chain `build_paper.py` currently selects is seeded but
data-dependent -- look for "PDB \_\_\_\_ chain \_\_\_\_, UniProt \_\_\_\_" in
Section 3.2 of the compiled paper, or `EX_PDB_ID`/`EX_CHAIN_ID` in its
console output, before assuming a particular ID.

**What the figure shows, and how:** three panels, one shared camera
orientation. (a) the experimental structure (all Phase-5-mapped residues),
cartoon representation, true interface residues green vs. non-interface
grey. (b) the same experimental structure, colored by PeSTo's per-residue
interface probability from the experimental input (viridis, 0-1, matching
build_paper.py's own colorbar exactly -- the RGB values are computed with
the project's own matplotlib in `.venv` and handed to PyMOL as literal
colors, not approximated by a PyMOL palette). (c) the AlphaFold model,
trimmed to the experimental chain's mapped UniProt range (the same
definition `src/models/run_pesto.py` uses for the `af_trimmed` input),
colored the same way by its own PeSTo predictions. The AlphaFold model in
(c) is superposed onto the experimental frame with the project's existing
Kabsch method (`src.analysis.structural_metrics.kabsch_superpose`'s
algorithm, computed once under `.venv` over the matched CA pairs and
applied to every atom of the AlphaFold model via a single
`get_coords`/`load_coords` round-trip in PyMOL) -- **not** PyMOL's own
`pair_fit` command, which was tried first and found to silently stop
applying its transform past roughly 110 explicit atom-pair arguments in
this PyMOL build (an internal selection-string length limit; it still
reports a plausible-looking RMSD even when it hasn't moved anything past
that point, so this is worth knowing if reusing `pair_fit` elsewhere on a
chain longer than ~110 residues).

## 7. Deviations from the original proposal

The [original proposal](docs/proposal.txt) differs from the implemented
plan in a few ways, documented in full in `PLAN.md`:

- Interfaces are computed directly from 3D coordinates rather than via the
  PISA web service (PISA has no practical batch API for this scale).
- Adding pLDDT as a model input does require a small architectural change
  (widening PeSTo's first layer), contrary to the proposal's framing.
- This project's numbers are not directly comparable to Yuan et al. (2024),
  who benchmarked PeSTo/ScanNet on ESMFold-predicted structures rather than
  AlphaFold — only the *direction* of the effect (a real accuracy drop) is
  expected to be consistent between the two.
- Given that prior AlphaFold-specific benchmarking work exists, this project
  is best described as the first systematic **AlphaFold-specific** PeSTo
  benchmark with a matched bound/unbound/AlphaFold design, rather than the
  first study of experimental-to-predicted transferability in general.

## 8. Key references

- Krapp, L.F., Abriata, L.A., Cortés Rodriguez, F. et al. PeSTo:
  parameter-free geometric deep learning for accurate prediction of protein
  binding interfaces. *Nat Commun* **14**, 2175 (2023).
  https://doi.org/10.1038/s41467-023-37701-8
- Tubiana, J., Schneidman-Duhovny, D., Wolfson, H.J. ScanNet: an
  interpretable geometric deep learning model for structure-based protein
  binding site prediction. *Nat Methods* **19**, 730–739 (2022).
  https://doi.org/10.1038/s41592-022-01490-7
- Varadi, M. et al. AlphaFold Protein Structure Database in 2024:
  providing structure coverage for over 214 million protein sequences.
  *Nucleic Acids Res.* **52**, D368–D375 (2024).
  https://doi.org/10.1093/nar/gkad1011
- Yuan, Q. et al. Genome-scale annotation of protein binding sites via
  language model and geometric deep learning (GPSite). *eLife* **13**,
  e93695 (2024). https://doi.org/10.7554/eLife.93695
