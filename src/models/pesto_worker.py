"""Standalone PeSTo inference worker (Phase 6).

Deliberately **not** part of the `src` package's normal import graph: it
is always invoked as a plain script (`python pesto_worker.py ...`), never
via `-m` or `import`. Reason: PeSTo's own code does absolute imports of a
top-level package literally named `src` (`from src.model_operations import
...`), which collides with *this* project's own top-level `src` package if
both are ever resolvable in the same process. Running this file as a
bare script keeps `sys.path[0]` at this file's own directory (never the
project root), so the two `src` packages never collide -- confirmed
empirically during Phase 6 development (see PROGRESS.md). Do not import
anything from this project's `src` package in this file.

One invocation = one job (one chain, one input variant) = one OS process,
by design: PeSTo's forward pass on the largest chains can exceed this
container's available memory, and an OS-level OOM-kill can only be
contained (and reported by the caller, without losing other jobs) if each
job is its own process. See src/models/run_pesto.py for the orchestrator.

Usage:
    python pesto_worker.py <pesto_root> <checkpoint_dir> <input_pdb> \
        <keys_csv> <output_parquet> <n_threads>

<keys_csv> has a header row naming the final key column(s) (e.g.
"auth_seq_id,auth_ins_code" or "uniprot_resnum") and one data row per
residue, in the *same order* the residues appear in <input_pdb>. PeSTo's
own preprocessing renumbers residues sequentially in file-encounter order
(verified: this holds even across modified/HETATM-flagged residues, whose
atoms can be physically reordered by PeSTo's own subunit-tagging step, but
whose *residue index* still sorts back to original file order -- see
PROGRESS.md), so output row i corresponds to <keys_csv> row i.
"""

from __future__ import annotations

import csv
import os
import sys


def main() -> int:
    pesto_root, checkpoint_dir, input_pdb, keys_csv, output_parquet, n_threads = sys.argv[1:7]
    n_threads = int(n_threads)

    sys.path.insert(0, checkpoint_dir)
    sys.path.insert(0, pesto_root)

    import torch as pt

    pt.set_num_threads(n_threads)

    from src.dataset import StructuresDataset, collate_batch_features
    from src.data_encoding import encode_structure, encode_features, extract_topology
    from src.structure import concatenate_chains
    from config import config_model
    from model import Model

    device = pt.device("cpu")
    model = Model(config_model)
    model.load_state_dict(pt.load(os.path.join(checkpoint_dir, "model_ckpt.pt"), map_location=device))
    model = model.eval().to(device)

    with open(keys_csv, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        keys = [row for row in reader]

    dataset = StructuresDataset([input_pdb], with_preprocessing=True)
    subunits, _ = dataset[0]
    if subunits is None:
        print(f"ERROR: failed to parse {input_pdb}", file=sys.stderr)
        return 2
    structure = concatenate_chains(subunits)

    X, M = encode_structure(structure)
    q = encode_features(structure)[0]
    ids_topk, _, _, _, _ = extract_topology(X, 64)
    X2, ids_topk2, q2, M2 = collate_batch_features([[X, ids_topk, q, M]])

    with pt.no_grad():
        z = model(X2, ids_topk2, q2, M2.float())
    # Channel 0 = protein-protein interface (config.py's r_types[0] is the
    # protein category; confirmed against PeSTo's own precomputed example
    # output during Phase 6 development -- see PROGRESS.md).
    probs = pt.sigmoid(z[:, 0]).cpu().numpy()

    if len(probs) != len(keys):
        print(
            f"ERROR: residue count mismatch: model produced {len(probs)}, keys file has {len(keys)}",
            file=sys.stderr,
        )
        return 3

    # Write output without importing pandas' parquet stack unnecessarily
    # heavy for a worker process -- pyarrow directly, same on-disk result.
    import pyarrow as pa
    import pyarrow.parquet as pq

    columns = {name: [] for name in header}
    for i, name in enumerate(header):
        columns[name] = [row[i] for row in keys]
    columns["pesto_interface_prob"] = [float(p) for p in probs]

    # Numeric key columns come through as strings from the CSV; restore
    # int dtype where every value parses cleanly (auth_seq_id, uniprot_resnum).
    for name in header:
        values = columns[name]
        if all(v == "" for v in values):
            continue
        try:
            columns[name] = [int(v) for v in values]
        except ValueError:
            pass  # leave as string (e.g. insertion codes)

    os.makedirs(os.path.dirname(output_parquet), exist_ok=True)
    table = pa.table(columns)
    pq.write_table(table, output_parquet)
    return 0


if __name__ == "__main__":
    sys.exit(main())
