"""Convert k-mesh convergence data to KSPACING equivalents and plot.

VASP formula: N_i = max(1, ceil(2*pi / (a_i * KSPACING)))
Inverse:      KSPACING = 2*pi / (a_i * N_i)

Usage:
    conda run -n pddmc python scripts/kspacing_analysis.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent.parent / "calculations" / "G1_bulk" / "convergence" / "results"

# G1 optimized lattice parameters
a_pd = 3.8907
a_pdo, c_pdo = 3.0536, 5.4058
a_pdo2, c_pdo2 = 4.5424, 3.1772

# k-mesh convergence data from G1 (dE vs densest ref)
pd_data = [
    ("6x6x6",     6, 6,  -9.81),
    ("8x8x8",     8, 8,  -1.68),
    ("10x10x10", 10, 10,  +1.96),
    ("12x12x12", 12, 12,  -1.00),
    ("14x14x14", 14, 14,  +0.18),
    ("16x16x16", 16, 16,   0.00),
]

pdo_data = [
    ("4x4x3",    4,  3, -0.84),
    ("6x6x4",    6,  4, -2.26),
    ("8x8x6",    8,  6, -1.21),
    ("10x10x8", 10,  8, -0.31),
    ("12x12x10",12, 10,  0.00),
]

pdo2_data = [
    ("4x4x4",    4,  4, -7.73),
    ("4x4x6",    4,  6, +4.04),
    ("6x6x8",    6,  8, -0.04),
    ("8x8x10",   8, 10, -0.66),
    ("10x10x12",10, 12,  0.00),
]


def calc_kspacing(a_vals, n_vals):
    """KSPACING = min over directions of 2*pi/(a_i * N_i)."""
    return min(2 * np.pi / (a * n) for a, n in zip(a_vals, n_vals))


def main():
    # Compute KSPACING for each k-mesh
    print("=" * 78)
    print(f"{'Material':8s} {'k-mesh':12s} {'KS_a':>8s} {'KS_c':>8s} "
          f"{'KS_min':>8s} {'dE(meV)':>8s}")
    print("-" * 78)

    pd_ks, pd_de = [], []
    for label, na, nc, de in pd_data:
        ks_a = 2 * np.pi / (a_pd * na)
        ks_min = ks_a  # cubic
        pd_ks.append(ks_min)
        pd_de.append(de)
        print(f"{'Pd':8s} {label:12s} {ks_a:8.4f} {'—':>8s} {ks_min:8.4f} {de:+8.2f}")
    print()

    pdo_ks, pdo_de = [], []
    for label, na, nc, de in pdo_data:
        ks_a = 2 * np.pi / (a_pdo * na)
        ks_c = 2 * np.pi / (c_pdo * nc)
        ks_min = min(ks_a, ks_c)
        pdo_ks.append(ks_min)
        pdo_de.append(de)
        print(f"{'PdO':8s} {label:12s} {ks_a:8.4f} {ks_c:8.4f} {ks_min:8.4f} {de:+8.2f}")
    print()

    pdo2_ks, pdo2_de = [], []
    for label, na, nc, de in pdo2_data:
        ks_a = 2 * np.pi / (a_pdo2 * na)
        ks_c = 2 * np.pi / (c_pdo2 * nc)
        ks_min = min(ks_a, ks_c)
        pdo2_ks.append(ks_min)
        pdo2_de.append(de)
        print(f"{'PdO2':8s} {label:12s} {ks_a:8.4f} {ks_c:8.4f} {ks_min:8.4f} {de:+8.2f}")

    print("=" * 78)
    print()
    print("Adopted meshes (red circles in plot):")
    print(f"  Pd   12x12x12  → KSPACING = {pd_ks[3]:.4f}")
    print(f"  PdO  8x8x6     → KSPACING = {pdo_ks[2]:.4f}")
    print(f"  PdO2 6x6x8     → KSPACING = {pdo2_ks[2]:.4f}")

    # --- Plot ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    fig.suptitle("k-mesh Convergence vs KSPACING (VASP equivalent)", fontsize=14)

    datasets = [
        ("Pd (fcc)", pd_ks, pd_de,
         [d[0] for d in pd_data], 3, "C0", "o"),
        ("PdO (P4$_2$/mmc)", pdo_ks, pdo_de,
         [d[0] for d in pdo_data], 2, "C1", "s"),
        ("PdO$_2$ (P4$_2$/mnm)", pdo2_ks, pdo2_de,
         [d[0] for d in pdo2_data], 2, "C2", "D"),
    ]

    for ax, (title, ks, de, labels, aidx, color, mkr) in zip(axes, datasets):
        ax.plot(ks, de, mkr + "-", color=color, ms=8, lw=1.5, zorder=3)

        # ±1 meV band
        ax.axhline(1.0, color="red", ls="--", alpha=0.5, lw=1)
        ax.axhline(-1.0, color="red", ls="--", alpha=0.5, lw=1)
        ax.axhspan(-1, 1, alpha=0.08, color="green", label="±1 meV/atom")

        # Adopted point
        ax.plot(ks[aidx], de[aidx], mkr, color="red", ms=14,
                mfc="none", mew=2.5, zorder=5, label=f"adopted ({labels[aidx]})")

        # Annotations
        for i in range(len(ks)):
            yoff = 12 if de[i] >= 0 else -16
            ax.annotate(labels[i], (ks[i], de[i]), fontsize=7,
                        textcoords="offset points", xytext=(0, yoff),
                        ha="center", rotation=25)

        ax.set_xlabel("KSPACING (Å⁻¹)", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.invert_xaxis()
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("ΔE (meV/atom)", fontsize=11)
    plt.tight_layout()

    out_path = OUT_DIR / "kspacing_convergence.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=150, bbox_inches="tight")
    print(f"\nPlot saved: {out_path}")


if __name__ == "__main__":
    main()
