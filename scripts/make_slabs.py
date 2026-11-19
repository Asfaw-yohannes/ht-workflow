#!/usr/bin/env python3
"""
Step 1: build slabs from relaxed bulk structures.

For every material in the materials CSV and every Miller index in config.yaml,
reads <bulk_dir>/<mp_id>.traj (or .cif), enumerates terminations and writes

    <output_dir>/<formula>_<mp_id>/<miller>/<termination>/
        <formula>_<mp_id>_<miller>_<termination>_initial.traj
        <formula>_<mp_id>_<miller>_<termination>.cif

plus a summary table (slabs_summary.csv) and, optionally, PNG previews.

Usage:
    python scripts/make_slabs.py --config config.yaml [--overwrite] [--png]
"""
import argparse
import csv
import sys
from pathlib import Path

from ase.io import read, write

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from htads.config import load_config, read_materials, material_name, miller_name  # noqa: E402
from htads.slabs import make_slabs  # noqa: E402


def find_bulk(bulk_dir, mp_id):
    for ext in (".traj", ".cif", ".vasp", ".json"):
        path = Path(bulk_dir) / f"{mp_id}{ext}"
        if path.is_file():
            return path
    raise FileNotFoundError(f"No bulk structure for {mp_id} in {bulk_dir} "
                            f"(expected {mp_id}.traj or .cif)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--overwrite", action="store_true",
                    help="rewrite slabs that already exist")
    ap.add_argument("--png", action="store_true",
                    help="also write PNG previews to <output_dir>/../slab_previews")
    args = ap.parse_args()

    cfg = load_config(args.config)
    scfg = cfg.get("slabs", {})
    out = Path(cfg["output_dir"])
    preview_dir = out.parent / "slab_previews"
    summary_rows = []

    for mp_id, formula in read_materials(cfg["materials_file"]):
        name = material_name(mp_id, formula)
        bulk = read(find_bulk(cfg["bulk_dir"], mp_id))
        for miller in cfg["miller_indices"]:
            mdir = out / name / miller_name(miller)
            try:
                slabs = make_slabs(bulk, miller,
                                   max_terminations=scfg.get("max_terminations", 4),
                                   min_slab_height=scfg.get("min_slab_height", 7.0),
                                   min_vacuum=scfg.get("min_vacuum", 17.0),
                                   min_atoms=scfg.get("min_atoms", 20))
            except Exception as err:          # keep going with other materials
                print(f"[FAIL] {name} {miller}: {err}")
                continue

            for s in slabs:
                tdir = mdir / s["label"]
                stem = f"{name}_{miller_name(miller)}_{s['label']}"
                traj = tdir / f"{stem}_initial.traj"
                if traj.exists() and not args.overwrite:
                    print(f"[skip] {traj} exists")
                else:
                    tdir.mkdir(parents=True, exist_ok=True)
                    s["atoms"].write(traj)
                    s["atoms"].write(tdir / f"{stem}.cif")
                    if args.png:
                        preview_dir.mkdir(parents=True, exist_ok=True)
                        write(preview_dir / f"{stem}.png", s["atoms"],
                              rotation="-80y,-80x")
                    print(f"[ok]   {traj}  ({s['info']['n_atoms']} atoms, "
                          f"{s['info']['n_fixed']} fixed)")
                summary_rows.append(dict(mp_id=mp_id, formula=formula,
                                         miller=miller_name(miller),
                                         termination=s["label"], **s["info"]))

    if summary_rows:
        summary = out.parent / "slabs_summary.csv"
        with open(summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary_rows[0]))
            w.writeheader()
            w.writerows(summary_rows)
        print(f"Summary written to {summary}")


if __name__ == "__main__":
    main()
