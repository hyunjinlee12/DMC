"""VESTA-style slab structure visualization from POSCAR files.

Renders all G2 slab structures with 3D perspective, shaded spheres,
and VESTA-like coloring. Produces both individual and combined figures.

Usage:
    conda run -n pddmc python scripts/visualize_slabs.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import to_rgba
from pathlib import Path
from ase.io import read
from ase.data import covalent_radii, atomic_numbers
from ase.data.colors import jmol_colors
from ase.constraints import FixAtoms

PROJECT = Path(__file__).resolve().parent.parent
G2_DIR = PROJECT / "calculations" / "G2_slab"
OUT_DIR = G2_DIR

ATOM_COLORS = {
    "Pd": np.array([0.0, 0.70, 0.72]),   # teal/cyan (VESTA-like)
    "O":  np.array([0.85, 0.15, 0.15]),   # red
}
ATOM_RADII = {"Pd": 0.85, "O": 0.55}

SLABS = [
    ("S1_Pd100",        "S1: Pd(100)\np(4×4) 5L",           "Pd⁰"),
    ("S2_PdO101_Pd100", "S2: PdO(101)/Pd(100)\n(√5×√5)R27° 2×2", "Pd⁰/Pd²⁺"),
    ("S3_PdO100",       "S3: PdO(100)\nO-rich 2×2",         "Pd²⁺"),
    ("S4_PdO2_110",     "S4: PdO₂(110)\nStoich 3×2",        "Pd⁴⁺"),
]


def render_atoms_3d(ax, atoms, rotation="-75x,15y,3z", scale=0.45,
                    show_cell=True, show_bonds=True, bond_cutoff=2.5):
    """Render atoms as depth-sorted shaded spheres with 3D perspective."""

    # Parse rotation
    angles = {}
    for part in rotation.split(","):
        part = part.strip()
        for axis in "xyz":
            if part.endswith(axis):
                angles[axis] = float(part[:-1])
                break
    rx = np.radians(angles.get("x", 0))
    ry = np.radians(angles.get("y", 0))
    rz = np.radians(angles.get("z", 0))

    Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
    Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
    Rz = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
    R = Rz @ Ry @ Rx

    pos = atoms.positions.copy()
    syms = atoms.get_chemical_symbols()

    # Get fixed atoms
    fix_idx = set()
    for con in atoms.constraints:
        if isinstance(con, FixAtoms):
            fix_idx.update(con.index)

    # Center positions
    center = pos.mean(axis=0)
    pos -= center

    # Rotate
    rpos = pos @ R.T

    # Sort by depth (back to front)
    depth = rpos[:, 2]
    order = np.argsort(depth)

    # Depth range for shading
    z_min, z_max = depth.min(), depth.max()
    z_range = max(z_max - z_min, 1.0)

    # Draw bonds first (behind atoms)
    if show_bonds:
        from ase.geometry import get_distances
        D = atoms.get_all_distances(mic=True)
        drawn_bonds = set()
        for i in order:
            for j in order:
                if i >= j:
                    continue
                bond_key = (min(i, j), max(i, j))
                if bond_key in drawn_bonds:
                    continue
                if D[i, j] < bond_cutoff and syms[i] != syms[j]:
                    drawn_bonds.add(bond_key)
                    x1, y1 = rpos[i, 0] * scale, rpos[i, 1] * scale
                    x2, y2 = rpos[j, 0] * scale, rpos[j, 1] * scale
                    mid_depth = (depth[i] + depth[j]) / 2
                    alpha = 0.15 + 0.25 * (mid_depth - z_min) / z_range
                    ax.plot([x1, x2], [y1, y2], color="#555555",
                            lw=0.4, alpha=alpha, zorder=0)

    # Draw atoms
    for idx in order:
        x, y = rpos[idx, 0] * scale, rpos[idx, 1] * scale
        sym = syms[idx]

        base_color = ATOM_COLORS.get(sym, np.array([0.5, 0.5, 0.5]))
        radius = ATOM_RADII.get(sym, 0.7) * scale * 1.8

        # Depth-based shading
        depth_frac = (depth[idx] - z_min) / z_range
        brightness = 0.6 + 0.4 * depth_frac
        alpha = 0.75 + 0.25 * depth_frac

        # Dimmer for fixed atoms
        if idx in fix_idx:
            brightness *= 0.75
            alpha *= 0.9

        color = np.clip(base_color * brightness, 0, 1)

        # Outer sphere (main body)
        circle = plt.Circle((x, y), radius, fc=color, ec="none",
                            alpha=alpha, zorder=int(depth_frac * 1000) + 1)
        ax.add_patch(circle)

        # Edge
        edge = plt.Circle((x, y), radius, fc="none",
                          ec=color * 0.4, lw=0.3, alpha=alpha * 0.8,
                          zorder=int(depth_frac * 1000) + 2)
        ax.add_patch(edge)

        # Specular highlight (small bright spot — VESTA style)
        hx = x - radius * 0.22
        hy = y + radius * 0.22
        hr = radius * 0.30
        highlight_color = np.clip(color + 0.5, 0, 1)
        hl = plt.Circle((hx, hy), hr, fc=highlight_color, ec="none",
                        alpha=alpha * 0.45, zorder=int(depth_frac * 1000) + 3)
        ax.add_patch(hl)

        # Secondary soft glow
        hx2 = x - radius * 0.10
        hy2 = y + radius * 0.10
        hr2 = radius * 0.55
        glow_color = np.clip(color + 0.2, 0, 1)
        gl = plt.Circle((hx2, hy2), hr2, fc=glow_color, ec="none",
                        alpha=alpha * 0.15, zorder=int(depth_frac * 1000) + 2)
        ax.add_patch(gl)

    # Cell outline
    if show_cell:
        cell = atoms.cell.copy()
        cell_r = (cell - center) @ R.T * scale  # not correct but approximate
        # Just draw a simple box outline at the slab boundary
        x_all = rpos[:, 0] * scale
        y_all = rpos[:, 1] * scale
        margin = 0.5
        xmin, xmax = x_all.min() - margin, x_all.max() + margin
        ymin, ymax = y_all.min() - margin, y_all.max() + margin

    ax.set_xlim(rpos[:, 0].min() * scale - 2, rpos[:, 0].max() * scale + 2)
    ax.set_ylim(rpos[:, 1].min() * scale - 2, rpos[:, 1].max() * scale + 2)
    ax.set_aspect("equal")
    ax.axis("off")


def make_combined_figure():
    """4-panel combined figure of all slabs."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    fig.suptitle("G2 Clean Slab Structures — 3D Perspective View",
                 fontsize=16, fontweight="bold", y=0.98)

    rotations = [
        "-70x,20y,5z",    # S1: slightly tilted
        "-70x,15y,3z",    # S2: similar angle
        "-75x,10y,0z",    # S3: more top-down
        "-70x,20y,5z",    # S4: tilted
    ]

    for idx, ((dirname, label, ox_state), rot) in enumerate(zip(SLABS, rotations)):
        ax = axes[idx // 2][idx % 2]
        poscar = G2_DIR / dirname / "POSCAR"
        if not poscar.exists():
            ax.text(0.5, 0.5, f"{dirname}\nNot found", transform=ax.transAxes,
                    ha="center", va="center")
            continue

        atoms = read(str(poscar))
        syms = atoms.get_chemical_symbols()
        n_pd = syms.count("Pd")
        n_o = syms.count("O")
        a = np.linalg.norm(atoms.cell[0])
        b = np.linalg.norm(atoms.cell[1])
        z = atoms.positions[:, 2]
        thick = z.max() - z.min()

        show_bonds = n_o > 0
        render_atoms_3d(ax, atoms, rotation=rot, scale=0.42,
                       show_cell=False, show_bonds=show_bonds,
                       bond_cutoff=2.3)

        ax.set_title(f"{label}\n{ox_state}", fontsize=11, fontweight="bold", pad=10)

        info = (f"{len(atoms)} atoms (Pd={n_pd}"
                + (f", O={n_o}" if n_o > 0 else "")
                + f")\n{a:.1f}×{b:.1f} Å, thick={thick:.1f} Å")
        ax.text(0.02, 0.02, info, transform=ax.transAxes, fontsize=8,
                va="bottom", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.85))

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=ATOM_COLORS["Pd"], markersize=12,
               markeredgecolor="gray", markeredgewidth=0.5, label="Pd"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=ATOM_COLORS["O"], markersize=9,
               markeredgecolor="gray", markeredgewidth=0.5, label="O"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=2,
               fontsize=11, frameon=True, bbox_to_anchor=(0.5, 0.01))

    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    out = OUT_DIR / "G2_slab_3D_view.png"
    plt.savefig(str(out), dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")
    plt.close()


def make_individual_figures():
    """Individual high-res figures for each slab with multiple views."""
    for dirname, label, ox_state in SLABS:
        poscar = G2_DIR / dirname / "POSCAR"
        if not poscar.exists():
            continue

        atoms = read(str(poscar))
        syms = atoms.get_chemical_symbols()
        show_bonds = syms.count("O") > 0

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        clean_label = label.split("\n")[0]
        fig.suptitle(f"{clean_label} — {ox_state}", fontsize=14, fontweight="bold")

        views = [
            ("Perspective", "-70x,20y,5z"),
            ("Top-down", "0x,0y,0z"),
            ("Side view", "-90x,0y,0z"),
        ]

        for ax, (view_name, rot) in zip(axes, views):
            render_atoms_3d(ax, atoms, rotation=rot, scale=0.40,
                           show_cell=False, show_bonds=show_bonds,
                           bond_cutoff=2.3)
            ax.set_title(view_name, fontsize=11)

        plt.tight_layout()
        out = OUT_DIR / f"{dirname}_views.png"
        plt.savefig(str(out), dpi=200, bbox_inches="tight", facecolor="white")
        print(f"Saved: {out}")
        plt.close()


if __name__ == "__main__":
    print("Generating 3D slab visualizations...")
    make_combined_figure()
    make_individual_figures()
    print("Done.")
