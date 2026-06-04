"""VASP I/O utilities for the Pd DMC project.

Handles POTCAR generation, KPOINTS, and ASE<->VASP interop.
Always use VASP_PP_PATH=/home/hyunjin/POTENTIAL (potpaw_PBE).
"""

import os
import shutil
from pathlib import Path

VASP_PP_PATH = Path("/home/hyunjin/POTENTIAL/potpaw_PBE")

POTCAR_MAP = {
    "Pd": "Pd_pv",
    "O": "O",
}


def write_potcar(elements, outdir):
    potcar_path = Path(outdir) / "POTCAR"
    with open(potcar_path, "w") as out:
        for el in elements:
            pp_name = POTCAR_MAP.get(el, el)
            src = VASP_PP_PATH / pp_name / "POTCAR"
            if not src.exists():
                raise FileNotFoundError(f"POTCAR not found: {src}")
            with open(src) as f:
                out.write(f.read())
    return potcar_path


def write_kpoints_gamma(kpts, outdir, label="Automatic mesh"):
    lines = [
        label,
        "0",
        "Gamma",
        f"  {kpts[0]}  {kpts[1]}  {kpts[2]}",
        "  0  0  0",
    ]
    path = Path(outdir) / "KPOINTS"
    path.write_text("\n".join(lines) + "\n")
    return path


def write_kpoints_mp(kpts, outdir, label="Automatic mesh"):
    lines = [
        label,
        "0",
        "Monkhorst-Pack",
        f"  {kpts[0]}  {kpts[1]}  {kpts[2]}",
        "  0  0  0",
    ]
    path = Path(outdir) / "KPOINTS"
    path.write_text("\n".join(lines) + "\n")
    return path


def get_unique_elements(atoms):
    """Get unique elements in order from ASE Atoms, preserving POSCAR order."""
    seen = []
    for s in atoms.get_chemical_symbols():
        if s not in seen:
            seen.append(s)
    return seen
