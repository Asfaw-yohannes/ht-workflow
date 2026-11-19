# ht-workflow

High-throughput workflow for DFT adsorption calculations on surfaces of bulk
materials, using ASE, pymatgen, VASP and SLURM.

```
bulk structure ─► relax_bulk.py ─► make_slabs.py ─► setup_adsorption.py ─► adsorption_energies.py
(CIF or MP id)    ISIF=3 relax     slabs, fixed     adsorbates on every    E_ads per site,
                  (VASP/SLURM)     bottom layers    site (VASP/SLURM)      lowest per surface
```

## Requirements

- Python ≥ 3.10, ASE ≥ 3.23, pymatgen ≥ 2024.1.27, NumPy ≥ 1.26, SciPy ≥ 1.11, PyYAML
  (tested with those minimum versions and with Python 3.13, ASE 3.29,
  pymatgen 2026.10.31, NumPy 2.5, SciPy 1.18)
- Optional: matplotlib (slab previews), mp-api (downloading bulks)
- VASP and its PBE pseudopotentials (not included; licensed separately)
- A SLURM cluster

```bash
conda env create -f environment.yml
conda activate htads
python -m pytest tests          # quick checks, no VASP needed
```

## Setup

1. `cp config.example.yaml config.yaml` and fill in your partition, account,
   modules, VASP command and POTCAR path. `config.yaml` is git-ignored.
2. List the materials in `data/mp_ids.csv` (`mp_id,formula`).
3. Put the unrelaxed bulk structures in `data/bulk_input/<mp_id>.cif`, or let
   `relax_bulk.py --fetch` download them from the Materials Project
   (`pip install mp-api`, and `export MP_API_KEY=...`; never put the key in a file
   in this repository). If you already have relaxed bulks, put them in
   `data/bulk/<mp_id>.traj` and skip step 0.
4. Check `gas_references` in `config.yaml`: gas-phase energies of CO, H2O and H2
   computed with the same VASP settings as the slabs.

## Running

```bash
# 0. Bulk relaxation, then copy the relaxed cells to data/bulk/
python scripts/relax_bulk.py --config config.yaml --submit [--fetch]
python scripts/relax_bulk.py --config config.yaml --collect

# 1. Slabs: runs/<formula>_<mp_id>/<miller>/<termination>/*_initial.traj
python scripts/make_slabs.py --config config.yaml --png

# 2. Adsorbate placement + job folders. Check first, then submit.
python scripts/setup_adsorption.py --config config.yaml --dry-run
python scripts/setup_adsorption.py --config config.yaml --submit

# 3. Job status and total energies
python scripts/collect_energies.py --config config.yaml

# 4. Adsorption energies
python scripts/adsorption_energies.py --config config.yaml
```

### Adsorption energies

`E_ads = E(slab+adsorbate) − E(slab) − Σ nᵢ μᵢ`, with μ_H = ½E(H₂),
μ_O = E(H₂O) − E(H₂) and μ_C = E(CO) − μ_O. So CO* is referenced to CO(g),
OH* to H₂O(g) − ½H₂(g), H* to ½H₂(g) and CHO* to CO(g) + ½H₂(g). These are
electronic energies (ΔE). For free energies (ΔG = ΔE + ΔZPE − TΔS), add the
corrections per adsorbate under `free_energy_corrections`;
`examples/free_energy_example.py` shows how to compute them from vibrational
frequencies and gas-phase ZPE/TS values. The applied correction is listed in
the `correction_eV` column.

Sites where the adsorbate or the slab atoms moved more than 1.5 Å
(`movement_limits`) are flagged as desorbed, dissociated or reconstructed and
left out of `adsorption_best.csv`, which lists the lowest energy per surface
and adsorbate.

Re-running step 2 is safe: jobs that are finished (`relaxed_pbe.e`) or
already submitted (`.submitted`) are skipped unless you pass `--force`.

Check `slabs_summary.csv` after step 1: it lists the atom count, number of
fixed atoms and in-plane size of every slab.

## Folder layout produced

```
runs/Cu3N_mp-1933/100/A/
├── Cu3N_mp-1933_100_A_initial.traj    # slab (bottom layers fixed)
├── CO/0/  relax.py  job.slurm  *_initial.traj  → relaxed_pbe.e, final.traj
├── CO/1/ ...
└── clean/1/                           # bare slab
```

## Repository contents

| Path | Purpose |
|---|---|
| `src/htads/slabs.py` | slab thickness/vacuum heuristics, terminations, tiling, fixed layers |
| `src/htads/atoms_ops.py` | orientation, adsorption sites, adsorbate placement, fingerprints |
| `src/htads/adsorbates.py` | starting geometries of adsorbates |
| `src/htads/kpoints.py` | k-point meshes for slabs and bulks |
| `src/htads/energies.py` | reference energies, adsorption energies, movement checks |
| `src/htads/jobs.py` | writes job folders, submits with `sbatch` |
| `templates/` | VASP relaxation script and SLURM job templates |
| `scripts/` | the workflow steps |
| `examples/free_energy_example.py` | ΔE → ΔG corrections (placeholder numbers) |

## Credits and license

Released under the GNU LGPL v3 (`COPYING.LESSER`, with `COPYING`).
Parts of the code are modified from
[GASpy](https://github.com/ulissigroup/GASpy) (Ulissi group, CMU), which is
also LGPL v3; see `NOTICE.md`.

If you use this workflow, please cite:

A. G. Yohannes, C. Lee, P. Talebi, D. H. Mok, M. Karamad, S. Back and
S. Siahrostami, "Combined High-Throughput DFT and ML Screening of Transition
Metal Nitrides for Electrochemical CO2 Reduction", *ACS Catalysis* **13**,
9007–9017 (2023). https://doi.org/10.1021/acscatal.3c01249

and GASpy, on which parts of the code are based:
K. Tran and Z. W. Ulissi, *Nature Catalysis* 1, 696–703 (2018).

**Data:** the adsorption energies from this paper are available on
[Catalysis-Hub](https://www.catalysis-hub.org/publications) (dataset
`YohannesCombined2023`).
