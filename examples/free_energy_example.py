#!/usr/bin/env python3
"""
Example: from adsorption energies (dE) to free energies (dG).

    dG = dE + dZPE - T dS

The correction for each adsorbate uses the same reference scheme as
adsorption_energies.py, applied to (ZPE - TS) instead of the DFT energy:

    corr(X*) = [ZPE - TS](X*) - sum_i n_i g_i
    g_H = 1/2 g(H2),  g_O = g(H2O) - g(H2),  g_C = g(CO) - g_O
    with g(gas) = [ZPE - TS](gas)

Adsorbates: ZPE and TS from their vibrational frequencies (harmonic
approximation, ase.thermochemistry.HarmonicThermo), e.g. from a VASP
IBRION=5 run with only the adsorbate atoms free.

*** ALL NUMBERS BELOW ARE ILLUSTRATIVE PLACEHOLDERS. ***
Replace them with your own frequencies and gas-phase ZPE/TS values
(the ones used in your paper) before using the output.

Run:
    python examples/free_energy_example.py
It prints dZPE - TdS per adsorbate and a `free_energy_corrections:` block
to paste into config.yaml. adsorption_energies.py then reports dG.
"""
import sys
from collections import Counter
from pathlib import Path

from ase.thermochemistry import HarmonicThermo

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from htads.energies import element_references  # noqa: E402

T = 298.15                       # K
CM_TO_EV = 1.239841984e-4        # 1 cm^-1 in eV

# Gas-phase ZPE and T*S at T and the chosen partial pressures (eV).
# PLACEHOLDERS: use the values from your own calculations / paper.
GAS = {
    "H2":  {"zpe": 0.27, "ts": 0.40},
    "H2O": {"zpe": 0.56, "ts": 0.67},
    "CO":  {"zpe": 0.13, "ts": 0.61},
}

# Adsorbate vibrational frequencies in cm^-1 (real modes only).
# PLACEHOLDERS: replace with frequencies from your vibrational calculations.
ADSORBATE_FREQS = {
    "H":   [1100, 700, 600],
    "OH":  [3600, 700, 650, 400, 350, 300],
    "CO":  [2000, 350, 300, 250, 150, 150],
    "CHO": [2800, 1700, 1300, 900, 600, 450, 350, 250, 150],
}
ADSORBATE_ATOMS = {"H": "H", "OH": "OH", "CO": "CO", "CHO": "CHO"}


def harmonic_zpe_ts(freqs_cm, temperature):
    """ZPE and T*S (eV) of an adsorbate in the harmonic approximation."""
    vib = [f * CM_TO_EV for f in freqs_cm]
    thermo = HarmonicThermo(vib_energies=vib)
    zpe = 0.5 * sum(vib)
    ts = temperature * thermo.get_entropy(temperature, verbose=False)
    return zpe, ts


def main():
    # Per-element (ZPE - TS), built exactly like the energy references
    g_gas = {name: v["zpe"] - v["ts"] for name, v in GAS.items()}
    g_el = element_references(g_gas)

    print(f"T = {T} K\n")
    print(f"{'adsorbate':10s}{'ZPE':>8s}{'TS':>8s}{'ref':>9s}{'dZPE-TdS':>11s}")
    corrections = {}
    for name, freqs in ADSORBATE_FREQS.items():
        zpe, ts = harmonic_zpe_ts(freqs, T)
        ref = sum(n * g_el[el] for el, n in Counter(ADSORBATE_ATOMS[name]).items())
        corrections[name] = (zpe - ts) - ref
        print(f"{name:10s}{zpe:8.3f}{ts:8.3f}{ref:9.3f}{corrections[name]:11.3f}")

    # Worked example: one site
    dE = -0.50                                    # e.g. CO* from adsorption_energies.csv
    print(f"\nExample: dE(CO*) = {dE:.2f} eV  ->  dG(CO*) = {dE + corrections['CO']:.2f} eV")

    print("\n# Paste into config.yaml (placeholder values!):")
    print("free_energy_corrections:")
    for name, c in corrections.items():
        print(f"  {name}: {c:.3f}")


if __name__ == "__main__":
    main()
