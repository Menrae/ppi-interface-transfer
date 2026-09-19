"""Renders Figure 6 (the example-chain structure figure) with a real
molecular renderer (PyMOL, cartoon representation) instead of
build_paper.py's matplotlib Calpha-trace fallback.

Three panels, same layout/colors/scale as the matplotlib fallback in
paper/build_paper.py's fig_example_structure: (a) true interface residues
(green) vs. non-interface (grey) on the experimental structure; (b) the
experimental structure colored by PeSTo's per-residue interface probability
from the experimental input (viridis, 0-1); (c) the AlphaFold model,
trimmed to the experimental chain's mapped UniProt range (same definition
src/models/run_pesto.py uses for the af_trimmed input) and superposed onto
the experimental structure's frame (same camera orientation as panel b, so
all three panels share one view), colored by PeSTo's per-residue interface
probability from the AlphaFold input.

Two-stage by necessity: PyMOL only exists in the isolated conda env at
external/pymol_render_env (see README.md), not in .venv, so this script
re-execs itself under that env's Python for the actual rendering once it
has prepared the small per-residue data file under .venv.

Usage: .venv/bin/python scripts/render_example_chain.py \
           --pdb-id 7DZ9 --chain-id C --output paper/figures/example_chain_render.png
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PYMOL_PYTHON = REPO_ROOT / "external" / "pymol_render_env" / "bin" / "python"

GREEN = "#009E73"
NON_INTERFACE_GREY = "#CCCCCC"
PAGE_WIDTH_IN = 6.6
RENDER_PX = (1500, 1300)  # per-panel PyMOL render size, downscaled into the composite


def _prepare_data(pdb_id: str, chain_id: str, work_dir: Path) -> Path:
    """Stage 1 (runs under .venv): reproduces build_paper.py's
    fig_example_structure data join, writes a per-residue CSV plus a small
    JSON of paths/metadata for stage 2 to consume."""
    sys.path.insert(0, str(REPO_ROOT))
    import numpy as np
    import pandas as pd
    from matplotlib import cm
    from src import config
    from src.analysis import structural_metrics as sme

    labels = pd.read_parquet(config.INTERIM_DATA_DIR / "interface_labels" / f"{pdb_id}_{chain_id}.parquet")
    exp_pred = pd.read_parquet(config.PROCESSED_DATA_DIR / "predictions" / "exp" / f"{pdb_id}_{chain_id}.parquet")
    af_pred = pd.read_parquet(config.PROCESSED_DATA_DIR / "predictions" / "af_trimmed" / f"{pdb_id}_{chain_id}.parquet")
    residue_map = pd.read_parquet(config.INTERIM_DATA_DIR / "residue_mappings" / f"{pdb_id}_{chain_id}.parquet")

    labels = labels.copy()
    labels["auth_ins_code"] = labels["auth_ins_code"].fillna("")
    exp_pred = exp_pred.copy()
    exp_pred["auth_ins_code"] = exp_pred["auth_ins_code"].fillna("")
    mapped = labels[labels["uniprot_resnum"].notna()].copy()
    mapped["uniprot_resnum"] = mapped["uniprot_resnum"].astype(int)

    merged = mapped.merge(exp_pred, on=["auth_seq_id", "auth_ins_code"], how="inner")
    merged = merged.merge(af_pred, on="uniprot_resnum", how="inner", suffixes=("_exp", "_af"))
    merged = merged.sort_values("auth_seq_id").reset_index(drop=True)

    if len(merged) == 0:
        raise ValueError(f"no residues survive the exp/af/labels join for {pdb_id}_{chain_id}")

    out = merged[["auth_seq_id", "auth_ins_code", "uniprot_resnum", "is_interface_contact"]].copy()
    out["exp_prob"] = merged["pesto_interface_prob_exp"]
    out["af_prob"] = merged["pesto_interface_prob_af"]
    out["is_interface_contact"] = out["is_interface_contact"].astype(int)
    # Viridis RGB computed here (matplotlib is not installed in the isolated
    # PyMOL env) so stage 2's coloring is pixel-identical to the colorbar
    # composited in stage 3 -- same colormap object, same 0-1 scale.
    exp_rgb = out["exp_prob"].clip(0, 1).map(lambda v: cm.viridis(v)[:3])
    af_rgb = out["af_prob"].clip(0, 1).map(lambda v: cm.viridis(v)[:3])
    out["exp_r"], out["exp_g"], out["exp_b"] = zip(*exp_rgb)
    out["af_r"], out["af_g"], out["af_b"] = zip(*af_rgb)
    residues_csv = work_dir / "residue_values.csv"
    out.to_csv(residues_csv, index=False)

    uniprot_acc = residue_map["uniprot_acc"].iloc[0]
    af_trim_lo = int(residue_map["uniprot_resnum"].min())
    af_trim_hi = int(residue_map["uniprot_resnum"].max())

    exp_cif = config.RAW_DATA_DIR / "pdb" / f"{pdb_id}_updated.cif.gz"
    af_cif = config.RAW_DATA_DIR / "alphafold" / f"{uniprot_acc}.cif.gz"
    if not exp_cif.exists():
        raise FileNotFoundError(exp_cif)
    if not af_cif.exists():
        raise FileNotFoundError(af_cif)

    # Kabsch rotation+translation (af onto exp) computed here, over the same
    # matched-residue CA set as the merged table, using the project's own
    # coordinate readers/superposition method (src/analysis/structural_metrics.py,
    # the same code build_paper.py's matplotlib fallback uses) -- PyMOL's own
    # `pair_fit` turns out to silently stop applying the transform past ~110
    # atom pairs in this build (an internal selection-string length limit),
    # so the fit is computed here instead and applied to all of afA's atoms
    # in stage 2 via a single get_coords/load_coords round-trip.
    auth_keys = list(zip(merged["auth_seq_id"].astype(int), merged["auth_ins_code"]))
    uniprot_keys = list(merged["uniprot_resnum"].astype(int))
    exp_coords_map = sme.read_ca_coords_by_auth(exp_cif, chain_id, auth_keys)
    af_coords_map = sme.read_ca_coords_by_uniprot(af_cif, uniprot_keys)
    exp_xyz = np.array([exp_coords_map[k] for k in auth_keys])
    af_xyz = np.array([af_coords_map[k] for k in uniprot_keys])

    af_mean = af_xyz.mean(axis=0)
    exp_mean = exp_xyz.mean(axis=0)
    mobile_c = af_xyz - af_mean
    target_c = exp_xyz - exp_mean
    h = mobile_c.T @ target_c
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    diag = np.diag([1.0, 1.0, d])
    r = vt.T @ diag @ u.T
    t = exp_mean - r @ af_mean
    rmsd = float(np.sqrt(np.mean(np.sum(((r @ af_xyz.T).T + t - exp_xyz) ** 2, axis=1))))
    print(f"Kabsch fit over {len(auth_keys)} matched CA pairs: RMSD {rmsd:.3f} A")

    meta = {
        "pdb_id": pdb_id, "chain_id": chain_id, "uniprot_acc": uniprot_acc,
        "exp_cif": str(exp_cif), "af_cif": str(af_cif),
        "af_trim_lo": af_trim_lo, "af_trim_hi": af_trim_hi,
        "residues_csv": str(residues_csv),
        "n_residues": len(out),
        "kabsch_r": r.tolist(), "kabsch_t": t.tolist(), "kabsch_rmsd": rmsd,
    }
    meta_path = work_dir / "meta.json"
    meta_path.write_text(json.dumps(meta, indent=2))
    print(f"Prepared {len(out)} residues (PDB {pdb_id} chain {chain_id}, UniProt {uniprot_acc}); "
          f"AlphaFold trim range [{af_trim_lo}, {af_trim_hi}]")
    return meta_path


def _kabsch_matrix_from_meta(meta: dict) -> tuple:
    import numpy as np
    r = np.array(meta["kabsch_r"])
    t = np.array(meta["kabsch_t"])
    return r, t


def _hex_to_rgb01(hex_color: str) -> tuple[float, float, float]:
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _render_panels(meta_path: Path, work_dir: Path) -> dict:
    """Stage 2 (runs under external/pymol_render_env's Python): loads the
    experimental and AlphaFold structures, superposes AlphaFold onto the
    experimental frame using the exact matched-residue CA pairs (same
    correspondence build_paper.py's Kabsch superposition uses), then
    ray-traces three cartoon panels sharing one camera orientation."""
    import csv
    import pymol2

    meta = json.loads(meta_path.read_text())
    rows = list(csv.DictReader(open(meta["residues_csv"])))
    for r in rows:
        r["auth_seq_id"] = int(r["auth_seq_id"])
        r["uniprot_resnum"] = int(r["uniprot_resnum"])
        r["is_interface_contact"] = int(r["is_interface_contact"])
        r["exp_prob"] = float(r["exp_prob"])
        r["af_prob"] = float(r["af_prob"])
        r["exp_rgb"] = [float(r["exp_r"]), float(r["exp_g"]), float(r["exp_b"])]
        r["af_rgb"] = [float(r["af_r"]), float(r["af_g"]), float(r["af_b"])]

    def resi_token(auth_seq_id: int, ins_code: str) -> str:
        return f"{auth_seq_id}{ins_code}" if ins_code else str(auth_seq_id)

    p = pymol2.PyMOL()
    p.start()
    cmd = p.cmd
    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.set("antialias", 2)
    cmd.set("cartoon_transparency", 0)
    cmd.set("ray_trace_mode", 0)
    # default specular puts a near-white hotspot on strand/helix tips that
    # face the light directly, which reads as a stray uncolored patch at
    # print size -- dial it down (verified this is a lighting effect, not a
    # coloring bug: every analyzed residue's cartoon color index was checked
    # and matched one of the 3 assigned colors, none left at PyMOL defaults)
    cmd.set("specular", 0.2)
    cmd.set("shininess", 20)

    all_resi = "+".join(resi_token(r["auth_seq_id"], r["auth_ins_code"]) for r in rows)

    cmd.load(meta["exp_cif"], "exp_full")
    # Restricted to exactly the analyzed (labels ∩ exp-pred ∩ af-pred)
    # residue set, not every observed residue in the crystal structure --
    # otherwise residues with no prediction/label (e.g. dropped upstream of
    # this join) show up as uncolored white cartoon with no explanation.
    cmd.create("expA_raw", f"exp_full and chain {meta['chain_id']} and polymer.protein")
    cmd.create("expA", f"expA_raw and resi {all_resi}")
    cmd.delete("exp_full")
    cmd.delete("expA_raw")

    cmd.load(meta["af_cif"], "af_full")
    cmd.create("afA_full", "af_full and polymer.protein")
    cmd.delete("af_full")
    cmd.create("afA", f"afA_full and resi {meta['af_trim_lo']}-{meta['af_trim_hi']}")
    cmd.delete("afA_full")

    if cmd.count_atoms("expA") == 0:
        raise RuntimeError("no atoms selected for the experimental chain")
    if cmd.count_atoms("afA") == 0:
        raise RuntimeError("no atoms selected for the trimmed AlphaFold model")

    # Apply the Kabsch fit computed in stage 1 to every atom of afA directly
    # (get_coords/load_coords preserve atom order within one round-trip) --
    # PyMOL's own `pair_fit` silently stops transforming past ~110 explicit
    # atom-pair arguments in this build, which a bisection confirmed hits an
    # internal selection-string length limit ("Selector-Error: Malformed
    # selection"), so it is not used for the full ~180-residue pair set.
    import numpy as np
    r_mat, t_vec = _kabsch_matrix_from_meta(meta)
    af_coords = np.array(cmd.get_coords("afA", state=1))
    af_coords_fit = af_coords @ r_mat.T + t_vec
    cmd.load_coords(af_coords_fit.tolist(), "afA", state=1)
    print(f"Applied stage-1 Kabsch fit to afA ({af_coords.shape[0]} atoms); "
          f"stage-1 CA RMSD was {meta['kabsch_rmsd']:.3f} A")

    interface_resi = [resi_token(r["auth_seq_id"], r["auth_ins_code"]) for r in rows if r["is_interface_contact"]]
    non_interface_resi = [resi_token(r["auth_seq_id"], r["auth_ins_code"]) for r in rows if not r["is_interface_contact"]]

    for obj in ("expA", "afA"):
        cmd.hide("everything", obj)
        cmd.show("cartoon", obj)
        cmd.color("grey70", obj)
    # afA is superposed onto expA (RMSD ~2-3 A over the whole chain, more at
    # flexible loops) but was staying visible+uncolored underneath expA for
    # panels (a)/(b), poking out as a stray grey/white sliver wherever the
    # two backbones diverge locally -- disable it now, only re-enable for
    # panel (c) where expA is disabled instead.
    cmd.disable("afA")
    # thicker loops and rounder helices read better at the printed panel
    # size (~1.7in wide) than PyMOL's thin-line cartoon defaults
    cmd.set("cartoon_loop_radius", 0.3)
    cmd.set("cartoon_tube_radius", 0.3)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_highlight_color", "grey50")

    cmd.set_color("true_interface_green", list(_hex_to_rgb01(GREEN)))
    cmd.set_color("non_interface_grey", list(_hex_to_rgb01(NON_INTERFACE_GREY)))
    if interface_resi:
        cmd.color("true_interface_green", "expA and resi " + "+".join(interface_resi))
    if non_interface_resi:
        cmd.color("non_interface_grey", "expA and resi " + "+".join(non_interface_resi))

    for r in rows:
        exp_name = f"expc_{r['auth_seq_id']}_{r['auth_ins_code'] or 'x'}"
        af_name = f"afc_{r['uniprot_resnum']}"
        cmd.set_color(exp_name, r["exp_rgb"])
        cmd.set_color(af_name, r["af_rgb"])

    # orient using the interface patch so it faces the camera, then zoom to
    # fit the whole chain -- this sets the one view shared by all 3 panels.
    if interface_resi:
        cmd.orient("expA and resi " + "+".join(interface_resi))
    else:
        cmd.orient("expA")
    cmd.zoom("expA", 4)
    cmd.turn("x", -15)

    # `orient` aligns the interface patch's principal axes with the screen
    # but doesn't guarantee which of the two resulting views has the patch
    # facing the camera vs. facing away (into the rest of the fold) -- check
    # explicitly and flip 180 degrees if the patch ended up on the far side.
    # Camera-space convention (larger Z = closer to viewer) confirmed
    # empirically: a probe point at model Z=+8 in front of an opaque
    # blocker at the origin rendered visible after cmd.reset(); after a
    # cmd.turn('y', 90), get_view()'s rotation R satisfies
    # camera_dir = R.reshape(3,3).T @ model_dir (note the transpose).
    if interface_resi:
        whole_com = np.array(cmd.centerofmass("expA"))
        interface_com = np.array(cmd.centerofmass("expA and resi " + "+".join(interface_resi)))
        outward = interface_com - whole_com
        rot = np.array(cmd.get_view()[:9]).reshape(3, 3)
        camera_dir_z = (rot.T @ outward)[2]
        if camera_dir_z < 0:
            cmd.turn("y", 180)
            print("Interface patch faced away from the camera after orient; flipped 180 deg.")

    view = cmd.get_view()
    cmd.set_view(view)  # afA shares this camera since it was superposed into expA's frame

    out_paths = {}

    # panel (a): true interface vs. non-interface
    png_a = work_dir / "panel_a.png"
    cmd.ray(RENDER_PX[0], RENDER_PX[1])
    cmd.png(str(png_a), dpi=300)
    out_paths["a"] = str(png_a)

    # panel (b): experimental structure colored by exp-input probability
    for r in rows:
        exp_name = f"expc_{r['auth_seq_id']}_{r['auth_ins_code'] or 'x'}"
        cmd.color(exp_name, f"expA and resi {resi_token(r['auth_seq_id'], r['auth_ins_code'])}")
    png_b = work_dir / "panel_b.png"
    cmd.ray(RENDER_PX[0], RENDER_PX[1])
    cmd.png(str(png_b), dpi=300)
    out_paths["b"] = str(png_b)

    # panel (c): AlphaFold (trimmed) colored by af-input probability, same view
    cmd.disable("expA")
    cmd.enable("afA")
    cmd.hide("everything", "expA")
    cmd.show("cartoon", "afA")
    for r in rows:
        af_name = f"afc_{r['uniprot_resnum']}"
        cmd.color(af_name, f"afA and resi {r['uniprot_resnum']}")
    cmd.set_view(view)
    png_c = work_dir / "panel_c.png"
    cmd.ray(RENDER_PX[0], RENDER_PX[1])
    cmd.png(str(png_c), dpi=300)
    out_paths["c"] = str(png_c)

    p.stop()
    return out_paths


def _compose_figure(panel_paths: dict, meta_path: Path, output_path: Path) -> None:
    """Runs under .venv: composites the three ray-traced PyMOL panels into
    one print-ready figure with a shared viridis colorbar and panel labels,
    matching build_paper.py's fig_example_structure layout exactly."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    from PIL import Image

    meta = json.loads(meta_path.read_text())

    fig, axes = plt.subplots(1, 3, figsize=(PAGE_WIDTH_IN, 2.6))
    titles = {"a": "(a) True interface", "b": "(b) Experimental", "c": "(c) AlphaFold (trimmed)"}
    for key, ax in zip("abc", axes):
        img = np.array(Image.open(panel_paths[key]).convert("RGB"))
        ax.imshow(img)
        ax.set_title(titles[key], fontsize=8.5 if key == "a" else 8)
        ax.axis("off")

    fig.subplots_adjust(wspace=0.05)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[axes[1], axes[2]], shrink=0.62, pad=0.03, label="Interface probability")
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label("Interface probability", fontsize=7.5)

    true_patch = mpatches.Patch(color=GREEN, label="True interface residue")
    non_patch = mpatches.Patch(color=NON_INTERFACE_GREY, label="Non-interface residue")
    axes[0].legend(handles=[true_patch, non_patch], loc="lower center", bbox_to_anchor=(0.5, -0.16),
                   fontsize=6.3, frameon=False, ncol=1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # dpi=450, not 300: bbox_inches="tight" crops this composite down to
    # ~5.3in wide, but the paper places it at 0.92*textwidth ~= 6.3in (see
    # template.tex), which would upscale a 300dpi save to an effective
    # ~250dpi at print size. 450 leaves that upscale still >=300 effective.
    fig.savefig(output_path, bbox_inches="tight", dpi=450, facecolor="white")
    plt.close(fig)
    print(f"Wrote {output_path} (PDB {meta['pdb_id']} chain {meta['chain_id']}, "
          f"{meta['n_residues']} residues, real PyMOL cartoon render)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--chain-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--work-dir", type=Path, default=None,
                         help="scratch dir for intermediate files (default: a temp dir)")
    parser.add_argument("--_stage2-meta", dest="stage2_meta", default=None,
                         help=argparse.SUPPRESS)  # internal: re-exec marker
    args = parser.parse_args()

    if args.stage2_meta is not None:
        # Running under external/pymol_render_env's Python (re-exec'd below).
        meta_path = Path(args.stage2_meta)
        work_dir = meta_path.parent
        panel_paths = _render_panels(meta_path, work_dir)
        (work_dir / "panels.json").write_text(json.dumps(panel_paths))
        return

    if not PYMOL_PYTHON.exists():
        raise SystemExit(
            f"PyMOL env not found at {PYMOL_PYTHON}. See README.md 'Substituting a rendered "
            f"structure figure' for how to create it (conda-forge pymol-open-source, isolated env)."
        )

    work_dir = args.work_dir or Path(tempfile.mkdtemp(prefix="render_example_chain_"))
    work_dir.mkdir(parents=True, exist_ok=True)

    meta_path = _prepare_data(args.pdb_id, args.chain_id, work_dir)

    subprocess.run(
        [str(PYMOL_PYTHON), str(Path(__file__).resolve()),
         "--pdb-id", args.pdb_id, "--chain-id", args.chain_id, "--output", str(args.output),
         "--_stage2-meta", str(meta_path)],
        check=True, cwd=str(REPO_ROOT),
    )

    panel_paths = json.loads((work_dir / "panels.json").read_text())
    _compose_figure(panel_paths, meta_path, args.output)


if __name__ == "__main__":
    main()
