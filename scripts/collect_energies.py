#!/usr/bin/env python3
"""
Step 3: collect energies and job status into one CSV.

Walks <output_dir>/<material>/<miller>/<termination>/<adsorbate>/<site>/ and
writes energies.csv with one row per job (status: done / submitted / new).

Usage:
    python scripts/collect_energies.py --config config.yaml [-o energies.csv]
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from htads.config import load_config  # noqa: E402
from htads import jobs  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("-o", "--output", default="energies.csv")
    args = ap.parse_args()

    out = Path(load_config(args.config)["output_dir"])
    rows = []
    for slurm in sorted(out.glob("*/*/*/*/*/job.slurm")):
        job_dir = slurm.parent
        material, miller, term, ads, site = job_dir.relative_to(out).parts
        status = jobs.job_status(job_dir)
        energy = ""
        if status == "done":
            energy = (job_dir / jobs.DONE_FILE).read_text().split()[0]
        rows.append(dict(material=material, miller=miller, termination=term,
                         adsorbate=ads, site=site, status=status, energy_eV=energy))

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["material", "miller", "termination",
                                           "adsorbate", "site", "status", "energy_eV"])
        w.writeheader()
        w.writerows(rows)
    n_done = sum(r["status"] == "done" for r in rows)
    print(f"{len(rows)} jobs, {n_done} done -> {args.output}")


if __name__ == "__main__":
    main()
