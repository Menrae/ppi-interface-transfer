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

# How Phase 7 benchmarking may filter out chains that overlap PeSTo's own
# training/model-selection data. "homolog": exclude anything with >=
# SEQUENCE_IDENTITY_CUTOFF identity to PeSTo's train+test+validation splits
# (primary analysis). "exact_train": exclude only exact PDBID_CHAINID
# matches to PeSTo's published training-set file (sensitivity A).
# "none": no leakage filter at all (sensitivity B). See PLAN.md Phase 2/7.
LEAKAGE_FILTER_MODES = ("homolog", "exact_train", "none")
