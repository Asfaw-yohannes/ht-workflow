"""Default settings for slab enumeration (passed to pymatgen's SlabGenerator).

Copyright (C) the GASpy contributors (Ulissi group, Carnegie Mellon University)
Modifications Copyright (C) 2026 Asfaw Yohannes

This file is a modified version of code from GASpy
(https://github.com/ulissigroup/GASpy). Modified November 2026: kept only the slab-generation settings from defaults.py.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation, version 3. It is distributed WITHOUT ANY WARRANTY; see
COPYING.LESSER and COPYING for details.
"""

from collections import OrderedDict


def slab_generator_settings():
    """Keyword arguments for `pymatgen.core.surface.SlabGenerator`."""
    return OrderedDict(min_slab_size=4,
                       min_vacuum_size=4,
                       lll_reduce=False,
                       center_slab=False,
                       primitive=True,
                       max_normal_search=10)


def get_slab_settings():
    """Keyword arguments for `SlabGenerator.get_slabs`."""
    return OrderedDict(tol=0.3,
                       bonds=None,
                       max_broken_bonds=0,
                       symmetrize=False)
