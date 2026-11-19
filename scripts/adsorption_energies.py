#!/usr/bin/env python3
"""
Step 4: adsorption energies.

For every finished adsorbate job
    <output_dir>/<material>/<miller>/<termination>/<adsorbate>/<site>/
and the matching clean slab  .../<termination>/clean/1/, computes

    E_ads = E(adslab) - E(slab) - sum_i n_i mu_i   (+ optional correction)

E_ads is the electronic adsorption energy (dE) when no corrections are set,
and the free energy (dG) when free_energy_corrections (dZPE - TdS per
adsorbate, see examples/free_energy_example.py) are given in config.yaml.
The applied correction is listed in the correction_eV column.

with mu_H = 1/2 E(H2), mu_O = E(H2O) - E(H2), mu_C = E(CO) - mu_O
(gas_references in config.yaml).

Each site is also checked for how far the adsorbate and the slab atoms moved
during relaxation (from final.traj, or the VASP output if that is missing);
sites beyond the limits are flagged (desorbed / dissociated / reconstructed).

Writes:
    adsorption_energies.csv   one row per site
    adsorption_best.csv       lowest unflagged energy per surface and adsorbate

Usage:
    python scripts/adsorption_energies.py --config config.yaml
"""
import argparse
import csv
import math
import sys
from pathlib import Path

from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from htads.config import load_config  # noqa: E402
from htads.adsorbates import get_adsorbate  # noqa: E402
from htads.energies import element_references, adsorption_energy, max_movements  # noqa: E402
from htads import jobs  # noqa: E402

FIELDS = ["material", "miller", "termination", "adsorbate", "site", "E_adslab_eV",
          "E_slab_eV", "correction_eV", "E_ads_eV", "max_ads_move_A", "max_slab_move_A", "flag"]


def read_energy(job_dir):
    return float((job_dir / jobs.DONE_FILE).read_text().split()[0])


def read_final(job_dir):
    """Relaxed structure in the original atom order, or None."""
    if (job_dir / "final.traj").exists():
        return read(job_dir / "final.traj")
    try:                                  # runs made with the original scripts
        from ase.calculators.vasp import Vasp
        return Vasp(restart=True, directory=str(job_dir)).get_atoms()
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("-o", "--output", default="adsorption_energies.csv")
    ap.add_argument("--best", default="adsorption_best.csv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    if "gas_references" not in cfg:
        sys.exit("config.yaml has no `gas_references` section (see config.example.yaml)")
    mu = element_references(cfg["gas_references"])
    corrections = cfg.get("free_energy_corrections") or {}
    limits = cfg.get("movement_limits", {})
    ads_max = limits.get("adsorbate", 1.5)
    slab_max = limits.get("slab", 1.5)
    out = Path(cfg["output_dir"])

    rows, n_skipped = [], 0
    for job_dir in sorted(p.parent for p in out.glob("*/*/*/*/*/" + jobs.DONE_FILE)):
        material, miller, term, ads_name, site = job_dir.relative_to(out).parts
        if ads_name == "clean":
            continue
        slab_dir = out / material / miller / term / "clean" / "1"
        if jobs.job_status(slab_dir) != "done":
            n_skipped += 1
            continue

        try:
            symbols = get_adsorbate(ads_name).get_chemical_symbols()
            e_adslab, e_slab = read_energy(job_dir), read_energy(slab_dir)
            corr = float(corrections.get(ads_name, 0.0))
            e_ads = adsorption_energy(e_adslab, e_slab, symbols, mu, corr)
        except (KeyError, ValueError) as err:
            print(f"[FAIL] {job_dir}: {err}")
            continue

        ads_move = slab_move = math.nan
        flag = ""
        initial = next(iter(job_dir.glob("*_initial.traj")), None)
        final = read_final(job_dir)
        if initial is not None and final is not None:
            ads_move, slab_move = max_movements(read(initial), final, len(symbols))
            problems = []
            if ads_move > ads_max:
                problems.append("adsorbate_moved")
            if slab_move > slab_max:
                problems.append("slab_moved")
            flag = ";".join(problems)
        else:
            flag = "no_final_structure"

        rows.append(dict(material=material, miller=miller, termination=term,
                         adsorbate=ads_name, site=site,
                         E_adslab_eV=f"{e_adslab:.6f}", E_slab_eV=f"{e_slab:.6f}",
                         correction_eV=f"{corr:.4f}",
                         E_ads_eV=f"{e_ads:.4f}",
                         max_ads_move_A=f"{ads_move:.3f}", max_slab_move_A=f"{slab_move:.3f}",
                         flag=flag))

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    # Lowest energy per surface and adsorbate, among sites without movement flags
    best = {}
    for r in rows:
        if r["flag"] in ("adsorbate_moved", "slab_moved", "adsorbate_moved;slab_moved"):
            continue
        key = (r["material"], r["miller"], r["termination"], r["adsorbate"])
        if key not in best or float(r["E_ads_eV"]) < float(best[key]["E_ads_eV"]):
            best[key] = r
    with open(args.best, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(best[k] for k in sorted(best))

    n_flag = sum(bool(r["flag"]) for r in rows)
    print(f"{len(rows)} sites ({n_flag} flagged) -> {args.output}; "
          f"{len(best)} surface/adsorbate minima -> {args.best}")
    if n_skipped:
        print(f"{n_skipped} sites skipped because their clean slab has not finished")


if __name__ == "__main__":
    main()
