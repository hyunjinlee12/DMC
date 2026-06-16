"""PdO2 lattice parameter comparison: Exp / MP(GGA-PBE) / this work PBE / PBE+D3.

Usage:
    conda run -n pddmc python scripts/plot_pdo2_pbe_comparison.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "calculations" / "G1_bulk"

# Data
# Exp:      Shaplygin et al. 1978 (ICSD 647283, high-pressure XRD)
# MP PBE:   Materials Project mp-1018886 (GGA-PBE relaxed)
# PBE:      this work, VASP PBE, ENCUT=520, k=6x6x8
# PBE+D3:   this work, VASP PBE+D3-BJ, ENCUT=520, k=6x6x8
labels  = ["Exp.\n(ICSD)", "MP PBE\n(mp-1018886)", "PBE\n(this work)", "PBE+D3\n(this work)"]
a_vals  = [4.4862,  4.5222,  4.5740,  4.5424]
c_vals  = [3.1032,  3.1463,  3.1914,  3.1772]
colors  = ["#555555", "#66C2A5", "#4C72B0", "#DD8452"]

fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
fig.suptitle("PdO$_2$ (rutile, P4$_2$/mnm): Lattice Parameter Comparison",
             fontsize=14, y=1.02)

# ── Panel 1: a, c absolute values ────────────────────────────────────
ax = axes[0]
x = np.arange(len(labels))
w = 0.3
bars_a = ax.bar(x - w/2, a_vals, w, color=colors, edgecolor="black", lw=0.8)
bars_c = ax.bar(x + w/2, c_vals, w, color=colors, edgecolor="black", lw=0.8,
                alpha=0.65, hatch="//")

for bar, val in zip(bars_a, a_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f"{val:.4f}", ha="center", va="bottom", fontsize=8, rotation=45)
for bar, val in zip(bars_c, c_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f"{val:.4f}", ha="center", va="bottom", fontsize=8, rotation=45)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("Lattice parameter (Å)", fontsize=11)
ax.set_title("a, c (Å)", fontsize=12)
ax.legend(["a-axis", "c-axis"], fontsize=9)
ax.set_ylim(2.5, 5.2)
ax.grid(axis="y", alpha=0.3)

# ── Panel 2: % error vs Exp ─────────────────────────────────────────
ax = axes[1]
comp_labels = ["MP PBE", "PBE\n(this work)", "PBE+D3\n(this work)"]
err_a = [(v - a_vals[0]) / a_vals[0] * 100 for v in a_vals[1:]]
err_c = [(v - c_vals[0]) / c_vals[0] * 100 for v in c_vals[1:]]
x2 = np.arange(len(comp_labels))
c_bar = colors[1:]

bars1 = ax.bar(x2 - w/2, err_a, w, color=c_bar, edgecolor="black", lw=0.8)
bars2 = ax.bar(x2 + w/2, err_c, w, color=c_bar, edgecolor="black", lw=0.8,
               alpha=0.65, hatch="//")

for bar, val in zip(bars1, err_a):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
            f"+{val:.2f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
for bar, val in zip(bars2, err_c):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
            f"+{val:.2f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold")

ax.axhline(2.0, color="red", ls="--", alpha=0.6, lw=1.2)
ax.text(2.55, 2.08, "2%", color="red", fontsize=9, alpha=0.8)
ax.set_xticks(x2)
ax.set_xticklabels(comp_labels, fontsize=8)
ax.set_ylabel("Error vs Exp. (%)", fontsize=11)
ax.set_title("Lattice Error vs Experiment", fontsize=12)
ax.legend(["2% threshold", "a-axis", "c-axis"], fontsize=8, loc="upper right")
ax.set_ylim(0, 3.5)
ax.grid(axis="y", alpha=0.3)

# ── Panel 3: Volume ──────────────────────────────────────────────────
ax = axes[2]
v_vals = [a**2 * c for a, c in zip(a_vals, c_vals)]
bars_v = ax.bar(x, v_vals, 0.45, color=colors, edgecolor="black", lw=0.8)

for bar, val in zip(bars_v, v_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.15,
            f"{val:.2f}", ha="center", va="bottom", fontsize=9)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("Volume (Å³)", fontsize=11)
ax.set_title("Cell Volume (Å³)", fontsize=12)
ax.set_ylim(58, 70)
ax.grid(axis="y", alpha=0.3)

# Volume error annotation
v_exp = v_vals[0]
box_text = "\n".join([f"{labels[i].replace(chr(10),' ')}: {(v_vals[i]-v_exp)/v_exp*100:+.1f}%"
                      for i in range(1, len(v_vals))])
ax.text(0.98, 0.95, box_text, transform=ax.transAxes, fontsize=8,
        va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", ec="gray", alpha=0.9))

# Footer
fig.text(0.5, -0.02,
         "Exp: Shaplygin et al., Zh. Neorg. Khim. 1978 (ICSD 647283)  |  "
         "MP: mp-1018886 (GGA-PBE)  |  "
         "This work: VASP ENCUT=520 eV, k=6×6×8",
         ha="center", fontsize=8, color="gray")

plt.tight_layout()
out_path = OUT / "PdO2_PBE_vs_D3.png"
plt.savefig(str(out_path), dpi=150, bbox_inches="tight")
print(f"Saved: {out_path}")


if __name__ == "__main__":
    pass
