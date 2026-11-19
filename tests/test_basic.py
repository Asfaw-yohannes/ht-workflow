"""Quick checks that run without VASP or SLURM:  python -m pytest tests"""
import sys
from pathlib import Path

import numpy as np
from ase import Atoms

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from htads.atoms_ops import add_adsorbate_onto_slab, remove_adsorbate  # noqa: E402
from htads.adsorbates import get_adsorbate, get_adsorbates  # noqa: E402
from htads.kpoints import slab_kpoints  # noqa: E402
from htads.slabs import make_slabs  # noqa: E402


def cu3n_bulk(a=3.815):
    return Atoms("NCu3", scaled_positions=[[0, 0, 0], [.5, 0, 0], [0, .5, 0], [0, 0, .5]],
                 cell=[a, a, a], pbc=True)


def test_tags_slab_first():
    slab = Atoms("Cu4", positions=[[0, 0, i] for i in range(4)], cell=[5, 5, 20], pbc=True)
    adslab = add_adsorbate_onto_slab(get_adsorbate("CO"), slab, [0, 0, 6])
    assert list(adslab.get_tags()) == [0, 0, 0, 0, 1, 1]
    assert adslab.get_chemical_symbols()[-2:] == ["C", "O"]
    bare, binding = remove_adsorbate(adslab)
    assert len(bare) == 4 and np.allclose(binding[1], [0, 0, 6])


def test_kpoints_ceil():
    for length, k in [(10.0, 2), (5.0, 4), (6.0, 4), (20.0, 1)]:
        atoms = Atoms(cell=[length, length, 30])
        assert slab_kpoints(atoms, 20) == [k, k, 1]


def test_adsorbates_have_atoms():
    for name, ads in get_adsorbates().items():
        assert len(ads) >= 1, name


def test_make_slabs_cu3n_100():
    slabs = make_slabs(cu3n_bulk(), (1, 0, 0))
    assert slabs
    for s in slabs:
        n, fixed = s["info"]["n_atoms"], s["info"]["n_fixed"]
        assert n >= 20 and 0 < fixed < n
        assert s["atoms"].cell[2, 2] > 0


def test_references_and_adsorption_energy():
    from htads.energies import element_references, reference_energy, adsorption_energy
    gas = {"CO": -14.7850961, "H2O": -14.22925213, "H2": -6.77325473}
    mu = element_references(gas)
    # CO* is referenced to CO(g); CHO* to CO(g) + 1/2 H2(g); OH* to H2O - 1/2 H2
    assert np.isclose(reference_energy(["C", "O"], mu), gas["CO"])
    assert np.isclose(reference_energy(["C", "H", "O"], mu), gas["CO"] + 0.5 * gas["H2"])
    assert np.isclose(reference_energy(["O", "H"], mu), gas["H2O"] - 0.5 * gas["H2"])
    assert np.isclose(adsorption_energy(-100.0, -85.0, ["C", "O"], mu), -100 + 85 + 14.7850961)


def test_movements_use_atom_order_not_tags():
    from htads.energies import max_movements
    init = Atoms("Cu2CO", positions=[[0, 0, 0], [2, 0, 0], [0, 0, 2], [0, 0, 3.2]],
                 cell=[10, 10, 20], pbc=True)
    init.set_tags([1, 1, 0, 0])          # reversed tags, as in the original runs
    fin = init.copy()
    fin.positions[2:] += [0, 0, 3.0]     # adsorbate moved 3 A
    ads, slab = max_movements(init, fin, n_adsorbate=2)
    assert np.isclose(ads, 3.0) and np.isclose(slab, 0.0)


def test_bulk_kpoints():
    from htads.kpoints import bulk_kpoints
    assert bulk_kpoints(cu3n_bulk(), 30) == [8, 8, 8]
