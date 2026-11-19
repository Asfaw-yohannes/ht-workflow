# License and third-party code

This repository is distributed under the **GNU Lesser General Public License
v3** (`COPYING.LESSER`), which builds on the GNU General Public License v3
(`COPYING`).

## GASpy

`src/htads/atoms_ops.py`, `src/htads/adsorbates.py` and
`src/htads/settings.py` are modified versions of code from **GASpy**
(Generalized Adsorption Simulator for Python), Ulissi group, Carnegie Mellon
University: https://github.com/ulissigroup/GASpy

GASpy is licensed under the GNU LGPL v3, so these modified files are
distributed under the same license. Each file states at the top that it was
modified, and how.

GASpy reference (as given in the GASpy README):
K. Tran and Z. W. Ulissi, "Active learning across intermetallics to guide
discovery of electrocatalysts for CO2 reduction and H2 evolution",
*Nature Catalysis* 1, 696–703 (2018).

Main modifications: removed the database/Luigi/FireWorks dependencies, fixed
the atom-tag order in `add_adsorbate_onto_slab`, updated imports for current
Python, SciPy and ASE, and moved the adsorbate geometries and slab settings
into separate modules.

## Dependencies (not included in this repository)

ASE (LGPL), pymatgen (MIT), NumPy and SciPy (BSD). VASP is commercial
software and is not included; neither are its pseudopotentials (POTCARs).
