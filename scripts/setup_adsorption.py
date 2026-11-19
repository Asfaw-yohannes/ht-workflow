#!/usr/bin/env python3
"""
Step 2: place adsorbates on every slab and (optionally) submit the VASP jobs.

For each slab <output_dir>/<material>/<miller>/<termination>/*_initial.traj
(written by make_slabs.py) this creates

    <termination>/<adsorbate>/<site index>/   relax.py, job.slurm, *_initial.traj
    <termination>/clean/1/                    the bare slab

Jobs that are finished (relaxed_pbe.e exists) or already submitted
(.submitted exists) are skipped unless --force is given.

Usage:
    python scripts/setup_adsorption.py --config config.yaml --dry-run   # write files only
    python scripts/setup_adsorption.py --config config.yaml --submit    # write + sbatch
"""
import argparse
import subprocess
import sys
from pathlib import Path

from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from htads.config import load_config, read_materials, material_name, miller_name  # noqa: E402
from htads.adsorbates import get_adsorbate  # noqa: E402
from htads.atoms_ops import find_adsorption_sites, add_adsorbate_onto_slab  # noqa: E402
from htads.kpoints import slab_kpoints  # noqa: E402
from htads import jobs  # noqa: E402


def find_slab_file(term_dir):
    files = sorted(term_dir.glob("*_initial.traj"))
    if len(files) != 1:
        raise FileNotFoundError(f"Expected one *_initial.traj in {term_dir}, found {len(files)}")
    return files[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="write job folders, do not submit")
    mode.add_argument("--submit", action="store_true", help="write job folders and run sbatch")
    ap.add_argument("--force", action="store_true",
                    help="rewrite/resubmit jobs that are done or already submitted")
    args = ap.parse_args()

    cfg = load_config(args.config)
    slurm_tmpl = jobs.load_template("job.slurm.tmpl")
    vasp_tmpl = jobs.load_template("vasp_relax.py.tmpl")
    out = Path(cfg["output_dir"])
    counts = {"written": 0, "submitted": 0, "skipped": 0, "failed": 0}

    def handle(job_dir, atoms, job_name):
        status = jobs.job_status(job_dir)
        if status != "new" and not args.force:
            counts["skipped"] += 1
            return
        kpts = slab_kpoints(atoms, cfg["kpoint_density"])
        jobs.write_job(job_dir, atoms, job_name, kpts, cfg, slurm_tmpl, vasp_tmpl)
        counts["written"] += 1
        if args.submit:
            try:
                job_id = jobs.submit(job_dir)
                counts["submitted"] += 1
                print(f"[sub]  {job_dir}  -> job {job_id}")
            except (subprocess.CalledProcessError, FileNotFoundError) as err:
                counts["failed"] += 1
                msg = getattr(err, "stderr", "") or err
                print(f"[FAIL] sbatch in {job_dir}: {msg}")
        else:
            print(f"[dry]  {job_dir}  kpts={kpts}")

    for mp_id, formula in read_materials(cfg["materials_file"]):
        name = material_name(mp_id, formula)
        for miller in cfg["miller_indices"]:
            mdir = out / name / miller_name(miller)
            if not mdir.is_dir():
                print(f"[warn] {mdir} not found - run make_slabs.py first")
                continue
            for term_dir in sorted(p for p in mdir.iterdir() if p.is_dir()):
                slab = read(find_slab_file(term_dir))
                tag = f"{name}-{miller_name(miller)}-{term_dir.name}"

                sites = find_adsorption_sites(slab)
                for ads_name in cfg["adsorbates"]:
                    adsorbate = get_adsorbate(ads_name)
                    for j, site in enumerate(sites):
                        adslab = add_adsorbate_onto_slab(adsorbate, slab, site)
                        handle(term_dir / ads_name / str(j), adslab, f"{tag}-{ads_name}-{j}")

                if cfg["include_clean_slab"]:
                    handle(term_dir / "clean" / "1", slab, f"{tag}-clean-1")

    print(", ".join(f"{k}: {v}" for k, v in counts.items()))


if __name__ == "__main__":
    main()
