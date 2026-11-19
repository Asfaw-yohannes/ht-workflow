"""Load the YAML configuration and the list of materials."""

import csv
from pathlib import Path

import yaml

REQUIRED_KEYS = ("materials_file", "output_dir", "miller_indices", "adsorbates",
                 "slurm", "vasp")


def load_config(path):
    """Read config.yaml and check that the required sections are present."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"Config file '{path}' not found. Copy config.example.yaml to "
            f"config.yaml and edit it for your cluster.")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    missing = [k for k in REQUIRED_KEYS if k not in cfg]
    if missing:
        raise KeyError(f"Missing keys in {path}: {missing}")
    cfg["miller_indices"] = [tuple(m) for m in cfg["miller_indices"]]
    cfg.setdefault("kpoint_density", 20.0)
    cfg.setdefault("include_clean_slab", True)
    return cfg


def read_materials(path):
    """
    Read the materials list: CSV with columns mp_id,formula.
    Blank lines, '#' comments and a header row are skipped.

    Returns a list of (mp_id, formula) tuples.
    """
    materials = []
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if not row or not row[0].strip() or row[0].strip().startswith("#"):
                continue
            if row[0].strip().lower() in ("mp_id", "mpid", "id"):
                continue
            if len(row) < 2:
                raise ValueError(f"Bad line in {path}: {row} (expected mp_id,formula)")
            materials.append((row[0].strip(), row[1].strip()))
    return materials


def material_name(mp_id, formula):
    """Directory name for a material, e.g. 'Cu3N_mp-1933'."""
    return f"{formula}_{mp_id}"


def miller_name(miller):
    """Directory name for a Miller index, e.g. (1, 1, 0) -> '110'."""
    return "".join(str(i) for i in miller)
