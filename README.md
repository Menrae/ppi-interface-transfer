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
AlphaFold's own confidence scores?** If the answer is "yes, and it's tied to
low pLDDT," a natural follow-up is whether giving PeSTo access to pLDDT as
an input can recover some of the lost accuracy — this project also makes a
first attempt at that.

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

*Last updated: 2026-09-16. See [`PROGRESS.md`](PROGRESS.md) for the
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
- ⬜ Phases 5–10 — not yet started.

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

### Running the completed phases

```bash
# Phase 1: candidate complex selection (writes data/interim/candidates.csv)
.venv/bin/python -m src.data.select_complexes --max-entries 500

# Phase 2: redundancy reduction + PeSTo leakage flagging
# (writes data/interim/candidates_dedup.csv, clusters.tsv, leakage_threshold_sweep.csv)
.venv/bin/python -m src.data.cluster_and_split

# Phase 3: download PDB (updated mmCIF) + AlphaFold DB structures
# (writes data/interim/fetch_report.csv, phase3_attrition.csv)
.venv/bin/python -m src.data.fetch_structures
# add --estimate-only to see the pre-flight time/disk projection without downloading

# Phase 4: residue-numbering mapping (experimental <-> AlphaFold)
# (writes data/interim/mapping_report.csv, residue_mappings/*.parquet)
.venv/bin/python -m src.data.align_residues

# Hand-verification: prints a table + ChimeraX commands for N random
# mapped primary-set chains (seeded, so re-running with the same --seed
# picks the same chains)
.venv/bin/python scripts/spot_check_mapping.py --n 3 --seed 0
```

All download-based commands (Phases 1-3) cache every downloaded file
under `data/raw/`, so rerunning them doesn't re-fetch anything that's
already on disk. Phase 4 is a purely local computation (no network).

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
