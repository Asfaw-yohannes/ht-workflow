#!/usr/bin/env python3
"""
Step 0: relax the bulk structures that make_slabs.py starts from.

Input structures (unrelaxed), per material in the materials CSV, in order:
  1. data/bulk_input/<mp_id>.cif / .traj / .vasp / .json, if present
  2. otherwise, with --fetch, the structure from the Materials Project
     (needs `pip install mp-api` and the MP_API_KEY environment variable)

Each bulk is relaxed in <bulk_runs_dir>/<formula>_<mp_id>/ with the
`bulk_vasp` settings from config.yaml (cell shape and volume relaxed, ISIF=3).

Usage:
    python scripts/relax_bulk.py --config config.yaml --dry-run [--fetch]
    python scripts/relax_bulk.py --config config.yaml --submit  [--fetch]
    python scripts/relax_bulk.py --config config.yaml --collect
        copies each finished final.traj to <bulk_dir>/<mp_id>.traj
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from htads.config import load_config, read_materials, material_name  # noqa: E402
from htads.kpoints import bulk_kpoints  # noqa: E402
from htads import jobs  # noqa: E402

INPUT_EXTS = (".cif", ".traj", ".vasp", ".json")


def local_input(input_dir, mp_id):
    for ext in INPUT_EXTS:
        path = Path(input_dir) / f"{mp_id}{ext}"
        if path.is_file():
            return path
    return None


def fetch_from_mp(mp_id, save_to):
    """Download a structure from the Materials Project and save it as CIF."""
    api_key = os.environ.get("MP_API_KEY")
    if not api_key:
        raise RuntimeError("Set the MP_API_KEY environment variable to use --fetch "
                           "(never write the key into a file in this repository).")
    try:
        from mp_api.client import MPRester
    except ImportError as err:
        raise RuntimeError("--fetch needs the mp-api package: pip install mp-api") from err
    with MPRester(api_key) as mpr:
        structure = mpr.get_structure_by_material_id(mp_id)
    save_to.parent.mkdir(parents=True, exist_ok=True)
    structure.to(filename=str(save_to))
    return save_to


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="write job folders only")
    mode.add_argument("--submit", action="store_true", help="write job folders and sbatch")
    mode.add_argument("--collect", action="store_true",
                      help="copy finished relaxed bulks into bulk_dir")
    ap.add_argument("--fetch", action="store_true",
                    help="download missing input structures from the Materials Project")
    ap.add_argument("--force", action="store_true",
                    help="rewrite/resubmit jobs that are done or already submitted")
    args = ap.parse_args()

    cfg = load_config(args.config)
    if "bulk_vasp" not in cfg:
        sys.exit("config.yaml has no `bulk_vasp` section (see config.example.yaml)")
    input_dir = Path(cfg.get("bulk_input_dir", "data/bulk_input"))
    runs_dir = Path(cfg.get("bulk_runs_dir", "runs_bulk"))
    bulk_dir = Path(cfg["bulk_dir"])
    density = cfg.get("bulk_kpoint_density", 30.0)

    if args.collect:
        bulk_dir.mkdir(parents=True, exist_ok=True)
        for mp_id, formula in read_materials(cfg["materials_file"]):
            job_dir = runs_dir / material_name(mp_id, formula)
            status = jobs.job_status(job_dir)
            if status != "done":
                print(f"[{status:9s}] {job_dir}")
                continue
            initial = read(next(job_dir.glob("*_initial.traj")))
            final = read(job_dir / "final.traj")
            vol_change = 100 * (final.get_volume() / initial.get_volume() - 1)
            dest = bulk_dir / f"{mp_id}.traj"
            if dest.exists() and not args.force:
                print(f"[exists   ] {dest} (use --force to overwrite)")
                continue
            shutil.copy(job_dir / "final.traj", dest)
            print(f"[collected] {dest}  volume change {vol_change:+.1f}%  "
                  f"a,b,c = {np.round(final.cell.lengths(), 3)}")
        return

    slurm_tmpl = jobs.load_template("job.slurm.tmpl")
    vasp_tmpl = jobs.load_template("vasp_relax.py.tmpl")
    counts = {"written": 0, "submitted": 0, "skipped": 0, "failed": 0}

    for mp_id, formula in read_materials(cfg["materials_file"]):
        name = material_name(mp_id, formula)
        job_dir = runs_dir / name
        if (bulk_dir / f"{mp_id}.traj").exists() and not args.force:
            print(f"[skip] {name}: relaxed bulk already in {bulk_dir}")
            counts["skipped"] += 1
            continue
        if jobs.job_status(job_dir) != "new" and not args.force:
            print(f"[skip] {name}: job {jobs.job_status(job_dir)}")
            counts["skipped"] += 1
            continue

        src = local_input(input_dir, mp_id)
        if src is None:
            if not args.fetch:
                print(f"[FAIL] {name}: no input structure in {input_dir} (or use --fetch)")
                counts["failed"] += 1
                continue
            try:
                src = fetch_from_mp(mp_id, input_dir / f"{mp_id}.cif")
            except Exception as err:
                print(f"[FAIL] {name}: {err}")
                counts["failed"] += 1
                continue

        atoms = read(src)
        atoms.pbc = True
        kpts = bulk_kpoints(atoms, density)
        jobs.write_job(job_dir, atoms, f"{name}-bulk", kpts, cfg, slurm_tmpl, vasp_tmpl,
                       vasp_settings=cfg["bulk_vasp"])
        counts["written"] += 1
        if args.submit:
            try:
                print(f"[sub]  {job_dir} -> job {jobs.submit(job_dir)}")
                counts["submitted"] += 1
            except (subprocess.CalledProcessError, FileNotFoundError) as err:
                print(f"[FAIL] sbatch in {job_dir}: {getattr(err, 'stderr', '') or err}")
                counts["failed"] += 1
        else:
            print(f"[dry]  {job_dir}  {len(atoms)} atoms  kpts={kpts}")

    print(", ".join(f"{k}: {v}" for k, v in counts.items()))


if __name__ == "__main__":
    main()
