"""Render G2 slabs with ASE POV-Ray textures (ase2, ase3).

Usage:
    conda run -n pddmc python scripts/render_povray.py
"""

import numpy as np
import subprocess
from pathlib import Path
from ase.io import read
from ase.io.pov import write_pov

PROJECT = Path(__file__).resolve().parent.parent
G2_DIR = PROJECT / "calculations" / "G2_slab"

SLABS = [
    ("S1_Pd100",        "S1: Pd(100)"),
    ("S2_PdO101_Pd100", "S2: PdO(101)/Pd(100)"),
    ("S3_PdO100",       "S3: PdO(100)"),
    ("S4_PdO2_110",     "S4: PdO₂(110)"),
]

TEXTURES = ["ase2", "ase3"]


def render_one(atoms, name, tex, rotation="-75x,15y,3z"):
    pov_path = G2_DIR / f"{name}_{tex}.pov"
    ini_path = G2_DIR / f"{name}_{tex}.ini"
    png_path = G2_DIR / f"{name}_{tex}.png"

    write_pov(
        str(pov_path),
        atoms,
        rotation=rotation,
        radii=0.85,
        show_unit_cell=0,
        povray_settings={
            "textures": [tex] * len(atoms),
            "canvas_width": 800,
            "transparent": False,
        },
    )

    # Fix POV-Ray 3.7 orthographic camera bug:
    # orthographic + right/up triggers "angle >= 180" parse error in 3.7
    # Replace with perspective camera
    txt = pov_path.read_text()
    import re
    txt = re.sub(
        r'camera \{orthographic\n\s*right ([^\n]+)\n\s*(direction [^\n]*\n)?\s*location ([^\}]+)\}',
        lambda m: f'camera {{perspective\n  angle 20\n  right {m.group(1)}\n  location {m.group(3)}}}',
        txt,
    )
    pov_path.write_text(txt)

    ret = subprocess.run(
        ["povray", ini_path.name],
        cwd=str(G2_DIR),
        capture_output=True, text=True,
    )
    if ret.returncode != 0:
        err_lines = ret.stderr.strip().split("\n")
        for l in err_lines[-5:]:
            print(f"  ERR: {l}")
        return None

    if png_path.exists():
        return png_path
    return None


def make_comparison():
    """Render S1 with both ase2 and ase3 for comparison."""
    poscar = G2_DIR / "S1_Pd100" / "POSCAR"
    atoms = read(str(poscar))

    results = {}
    for tex in TEXTURES:
        print(f"  Rendering S1 with {tex}...")
        png = render_one(atoms, "S1_Pd100", tex)
        if png:
            print(f"    OK: {png} ({png.stat().st_size/1024:.0f} KB)")
            results[tex] = png
        else:
            print(f"    FAILED")

    return results


def make_all_slabs(tex="ase3"):
    """Render all 4 slabs with chosen texture."""
    results = {}
    for dirname, label in SLABS:
        poscar = G2_DIR / dirname / "POSCAR"
        if not poscar.exists():
            print(f"  {dirname}: POSCAR not found, skipping")
            continue
        atoms = read(str(poscar))
        print(f"  Rendering {label} with {tex}...")
        png = render_one(atoms, dirname, tex)
        if png:
            print(f"    OK: {png.stat().st_size/1024:.0f} KB")
            results[dirname] = png
        else:
            print(f"    FAILED")
    return results


if __name__ == "__main__":
    print("=== Style comparison (S1 only) ===")
    make_comparison()
    print("\nDone.")
