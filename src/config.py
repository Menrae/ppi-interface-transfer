"""Central location for all filesystem paths and tunable thresholds.

No module outside this file should hardcode a path or a numeric threshold
used for filtering, binning, or cutoffs. Import from here instead.
"""

from pathlib import Path

# --- Paths -------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = REPO_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
INTERIM_DATA_DIR = DATA_DIR / "interim"
PROCESSED_DATA_DIR = DATA_DIR / "processed"

RESULTS_DIR = REPO_ROOT / "results"
NOTEBOOKS_DIR = REPO_ROOT / "notebooks"
LOGS_DIR = REPO_ROOT / "logs"
EXTERNAL_DIR = REPO_ROOT / "external"

# --- Structure quality / filtering thresholds ---------------------------

# Maximum crystallographic resolution (Angstrom) for an experimental
# structure to be included as a reference/training structure.
RESOLUTION_CUTOFF_ANGSTROM = 2.5

# AlphaFold2's training data cutoff, expressed as a PDB *release* date
# (rcsb_accession_info.initial_release_date), not deposition date — AF2 was
# trained on structures released on or before this date. Candidate
# complexes are selected with release_date > this cutoff, so AlphaFold has
# not seen them during training.
PDB_RELEASE_DATE_CUTOFF = "2018-04-30"

# AlphaFold per-residue confidence (pLDDT) band edges used to bucket
# predicted residues/interfaces into low/medium/high/very-high confidence.
PLDDT_BANDS = (50, 70, 90)

# Maximum heavy-atom distance (Angstrom) between residues on two chains
# for the pair to be considered part of the interface.
INTERFACE_DISTANCE_CUTOFF_ANGSTROM = 5.0

# Upper bound on distinct protein entities per PDB entry for candidate
# selection (Phase 1). Excludes very large assemblies (e.g. viral capsids,
# ribosomes) where "the interface" is no longer a single well-defined
# pairwise surface. Candidates have between 2 and this many protein
# entities, inclusive.
MAX_PROTEIN_ENTITIES = 10

# Minimum entity sequence length (entity_poly.rcsb_sample_sequence_length)
# for a chain to be considered a foldable domain rather than a peptide
# fragment. This is the polymer entity's full expressed/construct sequence
# length, NOT the number of residues actually resolved/modeled in the
# deposited coordinates -- confirmed on real data that these differ
# substantially and often (86% of Phase 1 candidates with a Phase 4 mapping
# have a different value; extreme case 7F90 chain B: 1,817-residue entity
# sequence vs. 45 residues actually observed). Comparable in spirit to
# PeSTo's own training filter (min_num_res=48 in their model/config.py);
# kept slightly more permissive here since this is candidate *selection*,
# not PeSTo's own training set.
MIN_CHAIN_LENGTH = 40

# Maximum pairwise sequence identity allowed between train and test
# entries to avoid redundancy/leakage.
SEQUENCE_IDENTITY_CUTOFF = 0.30

# Minimum alignment coverage (of the shorter sequence) required for an
# MMseqs2 hit to count, both for candidate redundancy clustering and for
# the PeSTo-overlap homology search. Guards against short local matches
# (e.g. a shared domain) being treated as full-chain redundancy/leakage.
CLUSTER_MIN_COVERAGE = 0.8

# How Phase 7 benchmarking may filter out chains that overlap PeSTo's own
# training/model-selection data. "homolog": exclude anything with >=
# SEQUENCE_IDENTITY_CUTOFF identity to PeSTo's train+test+validation splits
# (primary analysis). "exact_train": exclude only exact PDBID_CHAINID
# matches to PeSTo's published training-set file (sensitivity A).
# "none": no leakage filter at all (sensitivity B). See PLAN.md Phase 2/7.
LEAKAGE_FILTER_MODES = ("homolog", "exact_train", "none")

# Sample-size targets for the primary (leakage-filtered) benchmark set.
# See PLAN.md confound (e) / Phase 7: below MIN_TEST_CHAINS_FLOOR, don't
# loosen the leakage filter to hit the target -- stratify/report sensitivity
# modes instead.
MIN_TEST_CHAINS_TARGET = 100
MIN_TEST_CHAINS_FLOOR = 50

# Sequence-identity thresholds swept in Phase 2's leakage threshold-sweep
# report (data/interim/leakage_threshold_sweep.csv), to show how the
# survivor count against MIN_TEST_CHAINS_TARGET/FLOOR changes with the
# identity cutoff without pre-committing to a looser cutoff than
# SEQUENCE_IDENTITY_CUTOFF (0.30, the lowest value swept).
LEAKAGE_SWEEP_IDENTITY_THRESHOLDS = (0.30, 0.50, 0.70, 0.95)

# --- Phase 3: structure download -----------------------------------------

# Pre-flight estimate guardrails for src.data.fetch_structures: the module
# empirically calibrates a wall-clock/disk projection for the full download
# (see PLAN.md Phase 3) before committing to it. If the projection exceeds
# either budget, the run stops before the bulk download and asks for
# confirmation rather than proceeding unattended.
PHASE3_TIME_BUDGET_SECONDS = 2 * 60 * 60
PHASE3_DISK_BUDGET_BYTES = 20 * 1024**3

# Number of not-yet-cached unique PDB entries / UniProt accessions (each, in
# download-priority order) actually fetched to calibrate the pre-flight
# time/disk estimate before committing to the full run. These calibration
# fetches are real, cached downloads (not thrown away).
PHASE3_ESTIMATE_SAMPLE_SIZE = 15

# Polite fixed delay (seconds) between successive HTTP requests to PDBe,
# AlphaFold DB, and UniProt in src.data.fetch_structures.
PHASE3_REQUEST_DELAY_SECONDS = 0.1

# --- Phase 4: residue-numbering mapping ----------------------------------

# Minimum fraction of a chain's *observed* experimental residues that must
# end up "kept" (mapped to a UniProt position that is also within the
# AlphaFold model's range) for the chain to be retained. See PLAN.md Phase
# 4 for the justification and the coverage-threshold sweep this was chosen
# against -- 0.80 matches CLUSTER_MIN_COVERAGE's precedent elsewhere in
# this project for "a meaningful majority of the chain," while still
# tolerating some loss to expression tags/linkers and residues that fall
# outside the AlphaFold model's modeled range.
MIN_MAPPED_COVERAGE = 0.80

# Minimum fraction of a chain's UniProt-mapped-and-in-AlphaFold-range
# residues that must match the AlphaFold model's residue identity for a
# mapping *method* (SIFTS or the fallback alignment) to be trusted for
# that chain. Below this, isolated engineered point mutations (which are
# expected and should NOT fail the chain, see PLAN.md Phase 4) no longer
# explain the mismatch rate -- it's evidence of a numbering/frame error in
# that method, not biology, so the other method (or exclusion) is tried
# instead. 0.90 is deliberately generous to real single/few-residue
# mutations while still catching frame-shift-scale misalignment, which
# typically produces near-zero identity, not a borderline value.
RESIDUE_MAPPING_VALIDATION_IDENTITY = 0.90

# --- Phase 5: ground-truth interface labels ------------------------------

# Minimum ΔSASA (isolated chain minus the chain in its biological-assembly
# context), in square Angstrom, for a residue to count as buried by
# complex formation under the robustness (SASA-based) label. Given
# directly by the Phase 5 specification, not independently derived here.
SASA_BURIAL_CUTOFF_ANGSTROM2 = 1.0

# Minimum relative solvent accessibility (isolated chain) for a residue to
# count as "surface" when Phase 5/7 report interface fraction among
# surface residues specifically, as opposed to all residues. 0.25 is a
# widely used rule-of-thumb cutoff for "exposed" in the structural biology
# literature (residues below this are considered buried in the monomer
# even before any partner is added) -- see PLAN.md Phase 5 for the
# reasoning and citation.
SURFACE_RSA_THRESHOLD = 0.25

# --- Phase 6: PeSTo inference ---------------------------------------------

# Maximum heavy-atom count for a single PeSTo input (one chain, one
# variant). Above this, the job is skipped with reason
# "exceeds_memory_ceiling" rather than attempted -- profiled directly in
# this container (7.5 GB total RAM, no GPU): PeSTo's own
# src/data_encoding.extract_topology builds a full O(n^2) pairwise
# displacement/distance tensor, so peak memory grows worse than linearly
# with atom count. Measured peak RSS: 2338 atoms->1.0GB, 4855->2.1GB,
# 6041->3.0GB, 8833->6.2GB; 13822 atoms was OOM-killed outright. Solving
# for a peak-memory budget of ~4.8GB (this container's ~6.3GB "available"
# minus a ~1.5GB headroom, per session instruction) against that empirical
# curve gives ~7500 atoms -- used as the ceiling for every input variant
# (exp/af_full/af_trimmed) and also as the basis for the worker
# concurrency tiers in src/models/run_pesto.py (smaller inputs run with
# more concurrent workers, larger ones with fewer, so no combination of
# simultaneously-running jobs is projected to exceed this same budget).
PESTO_MAX_ATOMS = 7500

# Fixed, conservative concurrency schedule for src.models.run_pesto:
# (max_atoms_in_tier_exclusive, n_workers, n_torch_threads_per_worker).
# Tiers are processed in ascending order, one fully before the next, so
# memory from different tiers is never concurrent. Each tuple's
# n_workers * n_torch_threads_per_worker <= os.cpu_count() in this
# container (10), and n_workers * (worst-case peak memory for that tier's
# atom range, read off the PESTO_MAX_ATOMS profiling curve above) stays
# under the same ~4.8GB budget.
PESTO_CONCURRENCY_TIERS = (
    (2000, 4, 2),
    (5000, 2, 3),
    (PESTO_MAX_ATOMS, 1, 4),
)

# --- Phase 7: benchmarking -------------------------------------------------

# Number of chain-level bootstrap resamples (with replacement) for every
# 95% CI in src.analysis.benchmark. Chains, not residues, are resampled --
# residues within one chain aren't independent observations. 10,000 is
# comfortably enough for a stable percentile CI on ~500 chains and runs in
# well under a second.
BOOTSTRAP_N_RESAMPLES = 10_000

# Fixed seed for the chain-resampling bootstrap, so every CI is exactly
# reproducible.
PHASE7_BOOTSTRAP_SEED = 0

# Threshold on a chain's ground-truth interface fraction (Phase 5's
# interface_fraction_distance) used to split the pre-registered "small vs.
# large interface" descriptive stratum in Phase 7. 0.5 is the example
# given when this threshold was pre-registered (PLAN.md Phase 7 analysis
# plan) -- an even split of "less than half the chain is interface" vs.
# "at least half," not independently re-derived from the data.
INTERFACE_FRACTION_STRATUM_THRESHOLD = 0.5

# --- Phase 8: error analysis -----------------------------------------------

# Radius (Angstrom), in the experimental structure, defining a residue's
# local neighborhood for src.analysis.structural_metrics's local-RMSD
# feature: the set of mapped residues within this distance of a given
# residue's Calpha, each neighborhood superposed (Kabsch) on its own rather
# than reusing one whole-chain global fit. See PLAN.md Phase 8 analysis
# plan for why local (not global) superposition is used.
LOCAL_RMSD_RADIUS_ANGSTROM = 10.0

# Minimum number of Calpha-resolvable mapped residues a local-RMSD
# neighborhood must contain (including the residue itself) for the local
# Kabsch fit to be attempted -- below this a 3D rigid-body superposition is
# degenerate/underdetermined. Neighborhoods smaller than this get
# local_rmsd = NaN with a logged reason, not a silently-unstable fit.
LOCAL_RMSD_MIN_NEIGHBORS = 3

# Minimum number of residues a pLDDT band (or any other binned subset) must
# contain in src.analysis.error_analysis before its metrics are reported --
# below this, the bin is logged as too small rather than silently reported
# on a handful of residues. Matches the original Phase 8 sketch's success
# criterion.
PLDDT_BAND_MIN_RESIDUES = 30

# Variance-inflation-factor threshold above which a Phase 8 regression
# predictor is flagged as collinear in the output rather than silently
# reported as if independent. 5 is a standard rule-of-thumb VIF cutoff in
# the regression-diagnostics literature.
VIF_COLLINEARITY_THRESHOLD = 5.0

# Number of equal-sized chain-length quantile bins in Phase 8's
# chain-length-vs-AUPR-drop descriptive table. 4 (quartiles) gives a
# readable table without over-slicing the primary set's 554 scoreable
# chains.
LENGTH_ANALYSIS_N_QUANTILES = 4

# Chain-resampled bootstrap settings for Phase 8 (same method as Phase 7's
# BOOTSTRAP_N_RESAMPLES: resample chains, not residues, percentile CI). A
# separate seed from PHASE7_BOOTSTRAP_SEED so Phase 8's own reproducibility
# doesn't depend on Phase 7's call order.
PHASE8_BOOTSTRAP_SEED = 0

# Number of chain-resampled bootstrap draws for Phase 8's pLDDT-band pooled
# ROC-AUC/AUPR CIs specifically -- smaller than BOOTSTRAP_N_RESAMPLES
# (10,000, used for every other Phase 8 bootstrap, which resamples cheap
# per-chain scalars). Each band-metric draw here has to recompute a pooled
# O(n log n) ROC-AUC/AUPR over tens of thousands of pooled residues, not a
# median of ~500 scalars -- profiled empirically at ~0.02s/draw for a
# 50,000-residue band (both metrics), so 10,000 draws across every
# band x input combination would take on the order of 20-30 minutes; 2,000
# still gives a stable percentile CI at this n and keeps the full Phase 8
# run well under that.
PHASE8_BAND_BOOTSTRAP_N_RESAMPLES = 2_000

# --- Phase 8b: pLDDT-filtering post-hoc analysis (NOT pre-registered before
# Phase 7/8 were run -- see PLAN.md's post-hoc pre-registration note) ------

# Fixed probability threshold for turning PeSTo's continuous interface
# probability into a binary interface/non-interface call, used only for
# this precision/recall/F1 analysis (AUPR/ROC-AUC elsewhere are threshold-
# free). 0.5 is the natural midpoint of a probability output and is not
# tuned per input/cutoff -- deliberately a simple, fixed, un-optimized
# choice rather than each curve's own best operating point.
PLDDT_FILTER_PREDICTION_THRESHOLD = 0.5

# pLDDT cutoffs swept for the filtering analysis: keep only residues with
# pLDDT >= cutoff. 0 (keep everything, the no-filtering baseline) plus
# AlphaFold's own three confidence-band edges (PLDDT_BANDS).
PLDDT_FILTER_CUTOFFS = (0,) + PLDDT_BANDS
