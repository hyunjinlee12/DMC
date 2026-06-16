"""Test 6 color schemes on S3 PdO(100) with ase3 texture."""

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
from matplotlib.patches import Circle

G2 = Path(__file__).resolve().parent.parent / "calculations" / "G2_slab"

COLOR_SCHEMES = {
    "A: VESTA classic": {
        "Pd": [0.43, 0.48, 0.54],
        "O":  [1.00, 0.05, 0.05],
    },
    "B: Silver + Red": {
        "Pd": [0.75, 0.75, 0.78],
        "O":  [0.85, 0.15, 0.15],
    },
    "C: Blue-gray + Orange": {
        "Pd": [0.40, 0.50, 0.65],
        "O":  [0.90, 0.35, 0.10],
    },
    "D: Dark teal + Coral": {
        "Pd": [0.20, 0.45, 0.55],
        "O":  [0.95, 0.30, 0.30],
    },
    "E: Gunmetal + Crimson": {
        "Pd": [0.33, 0.33, 0.38],
        "O":  [0.80, 0.10, 0.20],
    },
    "F: Warm silver + Scarlet": {
        "Pd": [0.65, 0.63, 0.60],
        "O":  [0.90, 0.10, 0.10],
    },
}


def fix_pov(pov_path):
    txt = pov_path.read_text()
    txt = re.sub(
        r'camera \{orthographic\n\s*right ([^\n]+)\n\s*(direction [^\n]*\n)?\s*location ([^\}]+)\}',
        lambda m: f'camera {{perspective\n  angle 20\n  right {m.group(1)}\n  location {m.group(3)}}}',
        txt,
    )
    pov_path.write_text(txt)


def main():
    atoms = read(str(G2 / "S3_PdO100" / "POSCAR"))
    syms = atoms.get_chemical_symbols()

    results = {}
    for scheme_name, cmap in COLOR_SCHEMES.items():
        tag = scheme_name[0]
        colors_arr = [cmap.get(s, [0.5, 0.5, 0.5]) for s in syms]

        pov_path = G2 / f"_color_{tag}.pov"
        ini_path = G2 / f"_color_{tag}.ini"
        png_path = G2 / f"_color_{tag}.png"

        write_pov(
            str(pov_path), atoms,
            rotation="-75x,15y,3z", radii=0.85, show_unit_cell=0,
            colors=colors_arr,
            povray_settings={
                "textures": ["ase3"] * len(atoms),
                "canvas_width": 500,
                "transparent": False,
            },
        )
        fix_pov(pov_path)

        ret = subprocess.run(
            ["povray", ini_path.name], cwd=str(G2),
            capture_output=True, text=True,
        )
        if ret.returncode == 0 and png_path.exists():
            results[scheme_name] = png_path
            print(f"  {scheme_name}: OK")
        else:
            print(f"  {scheme_name}: FAILED")
            err = [l for l in ret.stderr.split("\n") if "Error" in l]
            if err:
                print(f"    {err[-1]}")

    # 2x3 comparison
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("S3 PdO(100) — Color Scheme Options (ase3)",
                 fontsize=14, fontweight="bold")

    for idx, (name, png_path) in enumerate(results.items()):
        ax = axes[idx // 3, idx % 3]
        img = mpimg.imread(str(png_path))
        ax.imshow(img)
        cmap = COLOR_SCHEMES[name]
        ax.set_title(name, fontsize=11, fontweight="bold")
        ax.add_patch(Circle((50, img.shape[0] - 40), 15,
                            fc=cmap["Pd"], ec="black", lw=1))
        ax.text(75, img.shape[0] - 40, "Pd", fontsize=10, va="center",
                fontweight="bold")
        ax.add_patch(Circle((140, img.shape[0] - 40), 12,
                            fc=cmap["O"], ec="black", lw=1))
        ax.text(162, img.shape[0] - 40, "O", fontsize=10, va="center",
                fontweight="bold")
        ax.axis("off")

    for idx in range(len(results), 6):
        axes[idx // 3, idx % 3].axis("off")

    plt.tight_layout()
    out = G2 / "color_comparison.png"
    plt.savefig(str(out), dpi=150, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {out}")
    plt.close()

    # cleanup
    for tag in "ABCDEF":
        for ext in [".pov", ".ini"]:
            f = G2 / f"_color_{tag}{ext}"
            if f.exists():
                f.unlink()


if __name__ == "__main__":
    print("=== Color scheme comparison ===")
    main()
