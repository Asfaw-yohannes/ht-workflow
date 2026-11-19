"""
Adsorption energies from total DFT energies.

    E_ads = E(slab + adsorbate) - E(slab) - sum_i n_i * mu_i

with elemental references built from gas-phase molecules:

    mu_H = 1/2 E(H2)
    mu_O = E(H2O) - E(H2)
    mu_C = E(CO) - mu_O          = E(CO) - E(H2O) + E(H2)
    mu_N = 1/2 E(N2)             (only if N2 is given)

So, for example, E_ads(CO*) is relative to CO(g), E_ads(OH*) to
H2O(g) - 1/2 H2(g), and E_ads(CHO*) to CO(g) + 1/2 H2(g).
These are electronic energies; zero-point and entropy corrections can be
added per adsorbate (free_energy_corrections in config.yaml).
"""
from collections import Counter

import numpy as np
from ase.geometry import find_mic


def element_references(gas):
    """Elemental chemical potentials (eV) from gas-phase total energies (eV)."""
    for key in ("H2", "H2O", "CO"):
        if key not in gas:
            raise KeyError(f"Missing gas reference '{key}' in config.yaml")
    mu = {"H": 0.5 * gas["H2"],
          "O": gas["H2O"] - gas["H2"]}
    mu["C"] = gas["CO"] - mu["O"]
    if "N2" in gas:
        mu["N"] = 0.5 * gas["N2"]
    return mu


def reference_energy(symbols, mu):
    """sum_i n_i mu_i for a list of chemical symbols."""
    counts = Counter(symbols)
    missing = [el for el in counts if el not in mu]
    if missing:
        raise KeyError(f"No reference energy for element(s) {missing}")
    return sum(n * mu[el] for el, n in counts.items())


def adsorption_energy(e_adslab, e_slab, adsorbate_symbols, mu, correction=0.0):
    return e_adslab - e_slab - reference_energy(adsorbate_symbols, mu) + correction


def max_movements(initial, final, n_adsorbate):
    """
    Largest displacement (A, minimum-image) of the adsorbate atoms and of the
    slab atoms between the initial and relaxed structures.

    The adsorbate is taken to be the LAST `n_adsorbate` atoms (adslab =
    slab + adsorbate). Tags are not used, because runs made with the original
    scripts have the slab/adsorbate tags reversed.
    """
    disp = final.positions - initial.positions
    _, dist = find_mic(disp, final.cell, final.pbc)
    n_slab = len(initial) - n_adsorbate
    ads, slab = dist[n_slab:], dist[:n_slab]
    return (float(ads.max()) if ads.size else 0.0,
            float(slab.max()) if slab.size else 0.0)
