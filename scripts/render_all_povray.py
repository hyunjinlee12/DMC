"""Render all G2 slabs with ase2 and ase3 POV-Ray textures, make comparison.

Usage:
    conda run -n pddmc python scripts/render_all_povray.py
"""

import re
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path
from ase.io import read
from ase.io.pov import write_pov

PROJECT = Path(__file__).resolve().parent.parent
G2 = PROJECT / "calculations" / "G2_slab"

SLABS = [
    ("S1_Pd100",        "S1: Pd(100)"),
    ("S2_PdO101_Pd100", "S2: PdO(101)/Pd(100)"),
    ("S3_PdO100",       "S3: PdO(100)"),
    ("S4_PdO2_110",     "S4: PdO₂(110)"),
]
TEXTURES = ["ase2", "ase3"]


def render_one(atoms, name, tex, rotation="-75x,15y,3z"):
    pov_path = G2 / f"{name}_{tex}.pov"
    ini_path = G2 / f"{name}_{tex}.ini"
    png_path = G2 / f"{name}_{tex}.png"

    write_pov(
        str(pov_path), atoms,
        rotation=rotation, radii=0.85, show_unit_cell=0,
        povray_settings={
            "textures": [tex] * len(atoms),
            "canvas_width": 600,
            "transparent": False,
        },
    )

    txt = pov_path.read_text()
    txt = re.sub(
        r'camera \{orthographic\n\s*right ([^\n]+)\n\s*(direction [^\n]*\n)?\s*location ([^\}]+)\}',
        lambda m: f'camera {{perspective\n  angle 20\n  right {m.group(1)}\n  location {m.group(3)}}}',
        txt,
    )
    pov_path.write_text(txt)

    ret = subprocess.run(
        ["povray", ini_path.name],
        cwd=str(G2),
        capture_output=True, text=True,
    )
    if ret.returncode != 0:
        err = ret.stderr.strip().split("\n")
        for l in err[-5:]:
            print(f"    ERR: {l}")
        return None
    return png_path if png_path.exists() else None


def main():
    results = {}
    for dirname, label in SLABS:
        poscar = G2 / dirname / "POSCAR"
        if not poscar.exists():
            print(f"  {dirname}: POSCAR not found")
            continue
        atoms = read(str(poscar))
        for tex in TEXTURES:
            print(f"  {dirname} [{tex}]...", end=" ", flush=True)
            png = render_one(atoms, dirname, tex)
            if png:
                results[(dirname, tex)] = png
                print(f"OK ({png.stat().st_size/1024:.0f} KB)")
            else:
                print("FAILED")

    # Comparison figure
    fig, axes = plt.subplots(4, 2, figsize=(10, 20))
    fig.suptitle("G2 Slab POV-Ray: ase2 vs ase3",
                 fontsize=16, fontweight="bold", y=0.98)

    for col, tex in enumerate(TEXTURES):
        axes[0, col].set_title(f"Texture: {tex}", fontsize=14,
                               fontweight="bold", pad=10)

    for row, (dirname, label) in enumerate(SLABS):
        for col, tex in enumerate(TEXTURES):
            ax = axes[row, col]
            key = (dirname, tex)
            if key in results and results[key].exists():
                img = mpimg.imread(str(results[key]))
                ax.imshow(img)
            else:
                ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                        ha="center", va="center", fontsize=14)
            ax.axis("off")
            if col == 0:
                ax.text(-0.05, 0.5, label, transform=ax.transAxes,
                        fontsize=12, fontweight="bold", rotation=90,
                        va="center", ha="right")

    plt.tight_layout(rect=[0.06, 0.01, 1, 0.96])
    out = G2 / "G2_ase2_vs_ase3.png"
    plt.savefig(str(out), dpi=150, bbox_inches="tight", facecolor="white")
    print(f"\nComparison: {out}")
    plt.close()


if __name__ == "__main__":
    print("=== Rendering all G2 slabs (ase2 + ase3) ===")
    main()
    print("Done.")
