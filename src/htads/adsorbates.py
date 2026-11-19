"""
Adsorbate geometries.

Each adsorbate is an `ase.Atoms` object whose FIRST atom is the binding atom,
placed at (or near) the origin and "pointing" upwards in +z. The adsorbate is
translated onto the adsorption site found by pymatgen, which already sits
roughly 2 A above the surface.

Note: earlier versions attached `ase.constraints.Hookean` springs to some
adsorbates (OOH, CHO, NO2, ...). Those were dropped here because (1) they were
lost when the adsorbate was added to the slab, and (2) the relaxation is done
inside VASP (IBRION=1), which only receives FixAtoms constraints from ASE.
If you need to keep adsorbates intact, relax with an ASE optimizer instead.

Copyright (C) the GASpy contributors (Ulissi group, Carnegie Mellon University)
Modifications Copyright (C) 2026 Asfaw Yohannes

This file is a modified version of code from GASpy
(https://github.com/ulissigroup/GASpy). Modified November 2026: moved adsorbate geometries out of defaults.py, removed Hookean constraints, added NO/NOH/HNO/NO2 variants.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation, version 3. It is distributed WITHOUT ANY WARRANTY; see
COPYING.LESSER and COPYING for details.
"""

from math import sqrt

from ase import Atoms

# Bond lengths (Angstrom) used to build some starting geometries
S_NO = 1.3     # N-O single bond (used for NOH)
S_OH = 0.94    # O-H bond
D_NO2R = 1.0   # height parameter used for the "r" (bridging) NO2 variants


def get_adsorbates():
    """Return a dict {name: ase.Atoms}. A fresh copy is built on every call."""
    ads = {}

    # Uranium is a place-holder used when fingerprinting sites
    ads["U"] = Atoms("U")

    # Atomic adsorbates (some placed slightly lower to help them bind)
    ads["H"] = Atoms("H", positions=[[0.0, 0.0, -0.5]])
    ads["N"] = Atoms("N", positions=[[0.0, 0.0, -1.0]])
    ads["O"] = Atoms("O")
    ads["C"] = Atoms("C")

    # Di- and polyatomics: binding atom at the origin
    ads["CO"] = Atoms("CO", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.2]])
    ads["OH"] = Atoms("OH", positions=[[0.0, 0.0, 0.0], [0.92, 0.0, 0.32]])
    ads["OOH"] = Atoms("OOH", positions=[[0.0, 0.0, 0.0],
                                         [1.28, 0.0, 0.67],
                                         [1.44, -0.96, 0.81]])
    # CHO: C binds (index 0), then H (index 1), then O (index 2).
    # Check: the starting C-O distance here is ~1.40 A (formyl C=O is ~1.2 A).
    ads["CHO"] = Atoms("CHO", positions=[[0.0, 0.0, 1.0],
                                         [-0.94, 0.2, 1.7],
                                         [0.986, 0.6, 1.8]])

    # Nitrogen-containing species
    ads["NO"] = Atoms("NO", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.23]])
    ads["NOH"] = Atoms("NOH", positions=[[0.0, 0.0, 0.0],
                                         [0.0, 0.0, S_NO],
                                         [0.0, S_OH / sqrt(2.0), S_NO + S_OH / sqrt(2.0)]])
    ads["HNO"] = Atoms("HNO", positions=[[0.0, 0.0, 0.0],
                                         [0.0, -0.55, 0.9526279441628825],
                                         [0.0, 1.0539529164056618, 0.6085]])
    ads["NO2"] = Atoms("NOO", positions=[[0.0, 0.0, 0.0],
                                         [0.8, 0.0, 0.8],
                                         [-0.8, 0.0, 0.8]])
    ads["NO2r"] = Atoms("NOO", positions=[[1.5 / 2, 0.0, D_NO2R * sqrt(2) / 2],
                                          [0.0, 0.0, 0.0],
                                          [1.5, 0.0, 0.0]])
    ads["NO2H"] = Atoms("NOOH", positions=[[0.0, 0.0, 0.0],
                                           [0.8, 0.0, 0.8],
                                           [-0.8, 0.0, 0.8],
                                           [0.8, 0.0, 1.6]])
    ads["HNO2"] = Atoms("NOOH", positions=[[1.5 / 2, 0.0, D_NO2R * sqrt(2) / 2],
                                           [0.0, 0.0, 0.0],
                                           [1.5, 0.0, 0.0],
                                           [1.5 / 2, 0.0, D_NO2R * sqrt(2) / 2 + 1]])
    ads["NO2Hr"] = Atoms("NOOH", positions=[[0.0, 0.0, 0.0],
                                            [0.8, 0.0, 0.8],
                                            [-0.8, 0.0, 0.8],
                                            [1.2, 0.0, 0.4]])
    return ads


def get_adsorbate(name):
    """Return a copy of one adsorbate, with a helpful error for unknown names."""
    ads = get_adsorbates()
    if name not in ads:
        raise KeyError(f"Unknown adsorbate '{name}'. Available: {sorted(ads)}")
    return ads[name].copy()
