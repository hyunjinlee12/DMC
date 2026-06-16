"""Render all 4 slabs with custom colors: Pd=#2a9d8f, O=#f07167."""

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

G2 = Path(__file__).resolve().parent.parent / "calculations" / "G2_slab"

COLORS = {
    "Pd": [0.165, 0.616, 0.561],  # #2a9d8f
    "O":  [0.941, 0.443, 0.404],  # #f07167
}

SLABS = [
    ("S1_Pd100",        "S1: Pd(100)"),
    ("S2_PdO101_Pd100", "S2: PdO(101)/Pd(100)"),
    ("S3_PdO100",       "S3: PdO(100)"),
    ("S4_PdO2_110",     "S4: PdO₂(110)"),
]


def fix_pov(pov_path):
    txt = pov_path.read_text()
    txt = re.sub(
        r'camera \{orthographic\n\s*right ([^\n]+)\n\s*(direction [^\n]*\n)?\s*location ([^\}]+)\}',
        lambda m: f'camera {{perspective\n  angle 20\n  right {m.group(1)}\n  location {m.group(3)}}}',
        txt,
    )
    pov_path.write_text(txt)


def render(atoms, name, rotation="-75x,15y,3z", width=600):
    syms = atoms.get_chemical_symbols()
    colors_arr = [COLORS.get(s, [0.5, 0.5, 0.5]) for s in syms]

    pov_path = G2 / f"_custom_{name}.pov"
    ini_path = G2 / f"_custom_{name}.ini"
    png_path = G2 / f"_custom_{name}.png"

    write_pov(
        str(pov_path), atoms,
        rotation=rotation, radii=0.85, show_unit_cell=0,
        colors=colors_arr,
        povray_settings={
            "textures": ["ase3"] * len(atoms),
            "canvas_width": width,
            "transparent": False,
        },
    )
    fix_pov(pov_path)

    ret = subprocess.run(["povray", ini_path.name], cwd=str(G2),
                         capture_output=True, text=True)
    for ext in [".pov", ".ini"]:
        (G2 / f"_custom_{name}{ext}").unlink(missing_ok=True)

    if ret.returncode == 0 and png_path.exists():
        return png_path
    return None


def main():
    results = {}
    for dirname, label in SLABS:
        poscar = G2 / dirname / "POSCAR"
        if not poscar.exists():
            continue
        atoms = read(str(poscar))
        print(f"  {label}...", end=" ", flush=True)
        png = render(atoms, dirname)
        if png:
            results[dirname] = (label, png)
            print("OK")
        else:
            print("FAILED")

    # 2x2 grid
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    fig.suptitle("G2 Slab Structures — Custom Colors (Pd: #2a9d8f, O: #f07167)",
                 fontsize=14, fontweight="bold", y=0.98)

    for idx, (dirname, label) in enumerate(SLABS):
        ax = axes[idx // 2, idx % 2]
        if dirname in results:
            lbl, png_path = results[dirname]
            img = mpimg.imread(str(png_path))
            ax.imshow(img)
            ax.set_title(lbl, fontsize=12, fontweight="bold")
            png_path.unlink()
        ax.axis("off")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = G2 / "G2_custom_colors.png"
    plt.savefig(str(out), dpi=200, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {out}")
    plt.close()


if __name__ == "__main__":
    print("=== Custom color rendering ===")
    main()
