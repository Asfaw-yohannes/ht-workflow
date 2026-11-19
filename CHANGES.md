# Changes from the original scripts

Private note for you; delete it before publishing if you prefer.

| Original file | New location |
|---|---|
| `make_slabs.py` | `scripts/make_slabs.py` + `src/htads/slabs.py` |
| `relax.py` | `scripts/setup_adsorption.py` + `src/htads/jobs.py`, `kpoints.py` |
| `adsorbate.template` | `templates/vasp_relax.py.tmpl` |
| `slurm.template` | `templates/job.slurm.tmpl` |
| `atoms_operators.py` | `src/htads/atoms_ops.py` |
| `defaults.py` | `src/htads/adsorbates.py`, `settings.py` (GASpy database settings dropped, unused) |
| `utils.py` | dropped (only `unfreeze_dict` was used; replaced by a 2-line helper) |
| `mp_ids.csv` / `mp_ids_main.csv` | `data/mp_ids.csv` (one file for both steps) |
| hard-coded partition, paths, VASP settings | `config.yaml` |

## New steps (not in the original scripts)

- `scripts/relax_bulk.py`: bulk relaxation (ISIF=3) from local CIFs or the
  Materials Project, and `--collect` to copy relaxed cells to `data/bulk/`.
  The `bulk_vasp` settings are suggested defaults, not taken from your paper.
- `scripts/adsorption_energies.py`: adsorption energies with your CO/H2O/H2
  references, movement flags, and the lowest-energy site per surface.
  It finds the adsorbate by atom order (last atoms), not tags, so it also
  works on runs made with the original scripts. For those runs it reads the
  relaxed structure from the VASP output, since they have no final.traj.

## Bugs fixed

- **Adsorbate tags reversed**: slab atoms were tagged 1 and adsorbate atoms 0.
  Energies were not affected, but anything that uses tags was.
- **SLURM exit code**: `exit $exit_code-` made every job exit with status 2.
- **SLURM settings ignored**: the template had no placeholders, so the
  partition/CPU settings in relax.py never reached the job.
- **NCORE=8 with 20 MPI tasks**: now 4 and configurable.
- **make_slabs crashed with more than one Miller index** (loop variable `k`
  was reused for the termination counter).
- **k-point rounding**: `round(x + 0.5)` replaced with `ceil(x)`. Results are
  identical except when 20/|a| is exactly an odd integer (e.g. 3.0), where the
  old code gave one extra k-point (4 instead of 3).
- Imports that fail on Python ≥ 3.10 / current ASE / SciPy 2.0.
- Hard-coded Materials Project API key removed. The formula is now taken from
  the CSV instead of a Materials Project query (the script used the legacy
  `rest/v1` API, which Materials Project has deprecated). **Revoke the old key.**

## Behaviour changes to be aware of

- Directory names now use the formula from `data/mp_ids.csv`. Previously
  make_slabs used the Materials Project formula; they should match, but check.
- Jobs that are finished or already submitted are skipped on re-runs.
- `sbatch` failures are reported and counted.
- The generated VASP script is always named `relax.py`; the job runs it with
  the configured Python instead of `./*.py`.
- Hookean springs on adsorbates were removed: they never reached VASP.
- Each job also writes `final.traj` (relaxed structure).

## Kept on purpose (so new slabs match your existing ones)

- The slab-thickness rules, the one-extra-step thickness and vacuum loops,
  the CIF round-trips, `max_normal_search=1`, the first 4 terminations, and
  tiling to at least 20 atoms. Verified on Cu3N (100), (110), (111): output
  identical to the original script.
- The height below which atoms are fixed is taken from the oriented unit cell,
  not the slab itself. This was correct for the Cu3N test slabs, but check
  `n_fixed` in `slabs_summary.csv` for terminations B, C, D.
- Your VASP settings, including `ediffg = -0.055555`, no `sigma`, and no
  dipole correction (commented suggestions are in config.example.yaml).
