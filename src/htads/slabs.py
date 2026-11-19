"""
Slab generation: choose slab thickness and vacuum, enumerate terminations,
orient, tile in-plane and fix the bottom layers.

This reproduces the logic of the original make_slabs script, so new slabs
match the ones already in your runs/ folder. Quirks that were kept on
purpose are marked "KEPT".
"""
import math
import os
import tempfile
import warnings

from ase.constraints import FixAtoms
from ase.io import read
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

from .atoms_ops import orient_atoms_upwards

TERMINATION_LABELS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _cif_roundtrip(obj):
    """
    Write a pymatgen Structure or ase.Atoms to CIF and read it back with ASE.
    KEPT: the original script did this, and it sets the cell orientation
    that the fixed-layer height is measured in.
    """
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "tmp.cif")
        if hasattr(obj, "to"):          # pymatgen Structure / Slab
            obj.to(filename=path)
        else:                           # ase.Atoms
            obj.write(path)
        return read(path)


def _slab_generator(structure, miller, min_slab, min_vac):
    # KEPT: max_normal_search=True in the original, which pymatgen treats as 1.
    return SlabGenerator(structure, miller,
                         min_slab_size=min_slab, min_vacuum_size=min_vac,
                         in_unit_planes=True, primitive=True,
                         max_normal_search=1)


def choose_min_slab_planes(p, n_atoms_ouc, h_ouc, n_desired_layers=2, default=3):
    """
    Heuristic for `min_slab_size` (in unit planes) from the oriented unit cell.

    Args:
        p            projected height / d_hkl (number of planes in the OUC)
        n_atoms_ouc  atoms in the oriented unit cell
        h_ouc        length of the OUC c vector (A)
    Returns:
        N (float), or `default` when no rule applies (a warning is issued).
    """
    def is_(x, y):
        return math.isclose(x, y, abs_tol=1e-6)

    if p > 2 and not is_(p, 2):
        base = round(n_desired_layers / (p * 2), 8)
        if 10 < n_atoms_ouc < 20:
            return base + (1 if h_ouc > 8 else 2)
        if n_atoms_ouc <= 10:
            return base + (2 if h_ouc > 8 else 3)
        return base + (1 if h_ouc > 8 else 2)          # n_atoms_ouc >= 20

    if is_(p, 2):
        base = round(n_desired_layers / p, 8)
        if 10 < n_atoms_ouc < 20:
            return base + (1 if h_ouc > 8 else 2)
        if n_atoms_ouc <= 10:
            return base + (2 if h_ouc > 8 else 3)
        return base + (1 if h_ouc > 8 else 2)

    if is_(p, 1):
        base = round(n_desired_layers / p, 8)
        if 10 < n_atoms_ouc < 20 or 5 < n_atoms_ouc < 10:
            return base
        if n_atoms_ouc <= 5:
            return base + (1 if h_ouc > 8 else 2)
        if n_atoms_ouc >= 20:
            n = base - 2
            return 1 if n == 0 else n
        # n_atoms_ouc == 10 falls through (also in the original script)

    warnings.warn(f"No slab-thickness rule for p={p}, n_atoms_ouc={n_atoms_ouc}; "
                  f"using min_slab_size={default}.")
    return default


def make_slabs(bulk_atoms, miller, max_terminations=4, min_slab_height=7.0,
               min_vacuum=17.0, min_atoms=20, max_atoms_untiled=25):
    """
    Build fixed, tiled, upward-oriented slabs for one bulk and Miller index.

    Returns:
        list of dicts with keys: label, atoms (ase.Atoms with FixAtoms), info
    """
    structure = AseAtomsAdaptor.get_structure(bulk_atoms)
    d = structure.lattice.d_hkl(miller)

    # 1) Oriented unit cell (OUC) and number of planes in it
    gen0 = _slab_generator(structure, miller, 3, 3)
    slabs0 = gen0.get_slabs()
    if not slabs0:
        raise RuntimeError(f"pymatgen found no slabs for {miller}")
    ouc = slabs0[0].oriented_unit_cell
    n_atoms_ouc = len(ouc)
    p = round(gen0._proj_height / d, 8)
    ouc_atoms = _cif_roundtrip(ouc)
    h_ouc = ouc_atoms.cell.lengths()[2]

    n_planes = choose_min_slab_planes(p, n_atoms_ouc, h_ouc)
    n_planes_initial = n_planes

    # 2) Slab thickness: grow until the OUC stack is taller than min_slab_height
    #    KEPT: the loop takes one extra step after the condition is met.
    n_ouc_slab = int(math.ceil(n_planes / p))
    h_stack = ouc_atoms.repeat((1, 1, n_ouc_slab)).cell.lengths()[2]
    n_slab = n_planes
    if n_atoms_ouc < 20:
        while h_stack <= min_slab_height:
            h_stack = ouc_atoms.repeat((1, 1, n_ouc_slab)).cell.lengths()[2]
            n_slab += 1
            n_ouc_slab = int(math.ceil(n_slab / p))

    # 3) Vacuum: grow until the OUC stack height (z) exceeds min_vacuum
    #    KEPT: same extra-step behaviour.
    n_vac = n_planes
    n_ouc_vac = int(math.ceil(n_vac / p))
    vac_size = ouc_atoms.repeat((1, 1, n_ouc_vac)).cell[2][2]
    while vac_size <= min_vacuum:
        n_ouc_vac = int(math.ceil(n_vac / p))
        vac_size = ouc_atoms.repeat((1, 1, n_ouc_vac)).cell[2][2]
        n_vac += 1

    # 4) Height below which atoms are fixed: top of (n-1) or (n-2) OUC repeats.
    #    Note: measured in the OUC frame, not on the slab itself. Check the
    #    fixed fraction printed in slabs_summary.csv, especially for
    #    terminations other than A.
    if n_ouc_slab == 1:
        z_fix = ouc_atoms.positions[:, 2].max() / 2
    elif n_ouc_slab <= 3:
        z_fix = ouc_atoms.repeat((1, 1, n_ouc_slab - 1)).positions[:, 2].max()
    else:
        z_fix = ouc_atoms.repeat((1, 1, n_ouc_slab - 2)).positions[:, 2].max()

    # 5) Enumerate terminations
    gen = _slab_generator(structure, miller, n_slab, n_vac)
    all_slabs = gen.get_slabs()
    if len(all_slabs) > max_terminations:
        warnings.warn(f"{len(all_slabs)} terminations found for {miller}; "
                      f"keeping the first {max_terminations}.")

    results = []
    repeats = [(1, 1), (1, 2), (2, 2), (2, 3), (3, 3)]
    for m, slab in enumerate(all_slabs[:max_terminations]):
        atoms = _cif_roundtrip(slab)
        atoms = _cif_roundtrip(orient_atoms_upwards(atoms))

        # In-plane tiling: smallest repeat giving >= min_atoms atoms
        n1, n2 = 1, 1
        if len(atoms) < max_atoms_untiled:
            for n1, n2 in repeats:
                if len(atoms) * n1 * n2 >= min_atoms:
                    break
        atoms = atoms.repeat((n1, n2, 1))

        mask = [z <= z_fix for z in atoms.positions[:, 2]]
        atoms.set_constraint(FixAtoms(mask=mask))

        n_fixed = sum(mask)
        inplane = atoms.cell.lengths()[:2]
        if n_fixed == 0 or n_fixed == len(atoms):
            warnings.warn(f"{miller} {TERMINATION_LABELS[m]}: {n_fixed}/{len(atoms)} "
                          f"atoms fixed - check this slab.")
        results.append(dict(
            label=TERMINATION_LABELS[m],
            atoms=atoms,
            info=dict(n_atoms=len(atoms), n_atoms_ouc=n_atoms_ouc, n_fixed=n_fixed,
                      p=p, n_planes_initial=n_planes_initial, n_planes_slab=n_slab,
                      n_planes_vacuum=n_vac, repeat=f"{n1}x{n2}",
                      a=round(float(inplane[0]), 3), b=round(float(inplane[1]), 3))))
    return results
