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

# Minimum modeled residues (entity_poly.rcsb_sample_sequence_length) for a
# chain to be considered a foldable domain rather than a peptide fragment.
# Comparable in spirit to PeSTo's own training filter (min_num_res=48 in
# their model/config.py); kept slightly more permissive here since this is
# candidate *selection*, not PeSTo's own training set.
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
