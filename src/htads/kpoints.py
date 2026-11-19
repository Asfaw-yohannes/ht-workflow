"""k-point mesh for slabs."""

import math

import numpy as np


def slab_kpoints(atoms, density=20.0):
    """
    Gamma-centred k-point mesh for a slab: ceil(density / |a_i|) along the two
    in-plane lattice vectors and 1 along the surface normal.

    Args:
        atoms    ase.Atoms slab (surface normal along the third cell vector)
        density  k-point density in A (k_i * |a_i| >= density)
    Returns:
        [k1, k2, 1]
    """
    cell = atoms.get_cell()
    k = [max(1, math.ceil(density / np.linalg.norm(cell[i]))) for i in range(2)]
    return [k[0], k[1], 1]


def bulk_kpoints(atoms, density=30.0):
    """Gamma-centred mesh for a bulk: ceil(density / |a_i|) along all three vectors."""
    cell = atoms.get_cell()
    return [max(1, math.ceil(density / np.linalg.norm(cell[i]))) for i in range(3)]
