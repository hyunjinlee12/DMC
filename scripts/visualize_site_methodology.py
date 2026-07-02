"""Visualize the site classification METHODOLOGY (not the picks).

For 4 representative cases (atop / bridge / hollow / physisorbed), show:
  - top-down structure
  - anchor atom (C of CO or O of methoxy) highlighted
  - 2.60 Å cutoff circle around anchor
  - substrate atoms colored:
       within cutoff  → red ring
       outside cutoff → gray
  - text annotation with neighbor count → label

Helps explain HOW site_type is determined.
"""
import glob, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from ase.io import read

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.size'] = 11

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
V2 = ROOT / 'calculations/G3_adsorption/DFT_shortlist_v2'
OUT = ROOT / 'reports/predft_advisor_figures/v2_site_renders'
OUT.mkdir(parents=True, exist_ok=True)

PD_BOND_CUTOFF = 2.60

# Pick 4 representative examples (different site types from v2 picks)
CASES = [
    ('S1/single_CO/00_single_CO_atop_Pd_idx00064.vasp',     'atop_Pd',     'CO'),
    ('S1/single_CO/01_single_CO_bridge_Pd-Pd_idx00008.vasp','bridge_Pd-Pd','CO'),
    ('S1/single_CO/02_single_CO_hollow_3Pd_idx00043.vasp', 'hollow_3Pd', 'CO'),
    ('S2/single_CO/01_single_CO_physisorbed_idx00130.vasp', 'physisorbed', 'CO'),
]


def find_anchor(atoms, kind):
    syms = atoms.get_chemical_symbols()
    if kind == 'CO':
        c_idx = [i for i, s in enumerate(syms) if s == 'C']
        return c_idx[-1]   # the C atom (anchor for CO)
    elif kind == 'CH3O':
        c_idx = [i for i, s in enumerate(syms) if s == 'C']
        c = c_idx[-1]
        o_to_c = [i for i, s in enumerate(syms) if s == 'O' and
                  atoms.get_distance(i, c, mic=True) < 1.7]
        return o_to_c[0]
    return None


def render_top_view(ax, atoms, anchor_idx, kind, label):
    """Top view with 2.6 Å cutoff circle. anchor atom red/larger."""
    syms = atoms.get_chemical_symbols()
    z = atoms.positions[:, 2]
    # Visualization box
    cell = atoms.cell
    a_vec, b_vec = cell[0][:2], cell[1][:2]
    # MIC distance from anchor for all atoms
    d_anchor = atoms.get_distances(anchor_idx, list(range(len(atoms))), mic=True)
    # Substrate = Pd + lattice O (not adsorbate H, not adsorbate C/O of CO or methoxy O)
    # For CO: substrate excludes C and the O bonded to C
    # For methoxy: substrate excludes C, O bonded to C, and H atoms
    n = len(atoms)
    if kind == 'CO':
        c_idx = [i for i in range(n) if syms[i] == 'C']
        c = c_idx[-1]
        o_to_c = [i for i in range(n) if syms[i] == 'O' and atoms.get_distance(i, c, mic=True) < 1.4]
        ads_set = {c} | set(o_to_c)
    else:
        c_idx = [i for i in range(n) if syms[i] == 'C']
        c = c_idx[-1]
        h_idx = [i for i in range(n) if syms[i] == 'H']
        o_to_c = [i for i in range(n) if syms[i] == 'O' and atoms.get_distance(i, c, mic=True) < 1.7]
        ads_set = {c, o_to_c[0]} | set(h_idx)

    sub_idx = [i for i in range(n) if i not in ads_set]

    # Draw all substrate atoms (top view = xy)
    # Project around anchor's position
    anchor_xy = atoms.positions[anchor_idx, :2]

    # Use MIC vector to anchor for visualization
    for i in sub_idx:
        # MIC vector from anchor to i
        d_vec = atoms.positions[i] - atoms.positions[anchor_idx]
        for k in [0, 1]:
            while d_vec[k] > cell[k][k]/2: d_vec[k] -= cell[k][k]
            while d_vec[k] < -cell[k][k]/2: d_vec[k] += cell[k][k]
        plot_x = anchor_xy[0] + d_vec[0]
        plot_y = anchor_xy[1] + d_vec[1]
        d_total = np.sqrt(d_vec[0]**2 + d_vec[1]**2 + d_vec[2]**2)

        if syms[i] == 'Pd':
            color = '#3070b0' if d_total < PD_BOND_CUTOFF else '#bbbbbb'
            size = 380 if d_total < PD_BOND_CUTOFF else 230
            edge = 'red' if d_total < PD_BOND_CUTOFF else 'black'
            ew = 2.5 if d_total < PD_BOND_CUTOFF else 0.5
        else:  # O
            color = '#e74c3c' if d_total < PD_BOND_CUTOFF else '#fadbd8'
            size = 220 if d_total < PD_BOND_CUTOFF else 150
            edge = 'red' if d_total < PD_BOND_CUTOFF else 'black'
            ew = 2.5 if d_total < PD_BOND_CUTOFF else 0.5
        ax.scatter(plot_x, plot_y, s=size, c=color, edgecolor=edge,
                   linewidth=ew, zorder=3 if d_total < PD_BOND_CUTOFF else 2)
        # Distance annotation for nearby atoms
        if d_total < PD_BOND_CUTOFF:
            ax.annotate(f'{d_total:.2f}', (plot_x, plot_y), xytext=(3, 3),
                       textcoords='offset points', fontsize=8, color='red',
                       fontweight='bold', zorder=5)

    # Draw adsorbate atoms (CO or methoxy)
    for i in ads_set:
        d_vec = atoms.positions[i] - atoms.positions[anchor_idx]
        for k in [0, 1]:
            while d_vec[k] > cell[k][k]/2: d_vec[k] -= cell[k][k]
            while d_vec[k] < -cell[k][k]/2: d_vec[k] += cell[k][k]
        plot_x = anchor_xy[0] + d_vec[0]
        plot_y = anchor_xy[1] + d_vec[1]
        if i == anchor_idx:
            color = 'gold'; size = 500; edge = 'black'; ew = 2.5
        elif syms[i] == 'C':
            color = '#1a1a1a'; size = 200; edge = 'black'; ew = 1
        elif syms[i] == 'O':
            color = '#c0392b'; size = 200; edge = 'black'; ew = 1
        else:  # H
            color = 'white'; size = 100; edge = 'black'; ew = 0.5
        ax.scatter(plot_x, plot_y, s=size, c=color, edgecolor=edge,
                   linewidth=ew, zorder=6)

    # Draw 2.60 Å cutoff circle around anchor
    circle = Circle(anchor_xy, PD_BOND_CUTOFF, fill=False, edgecolor='red',
                    linewidth=2, linestyle='--', zorder=4)
    ax.add_patch(circle)
    # Indicate cutoff radius
    ax.annotate('', xy=(anchor_xy[0]+PD_BOND_CUTOFF, anchor_xy[1]),
                xytext=anchor_xy,
                arrowprops=dict(arrowstyle='->', color='red', lw=1.2))
    ax.text(anchor_xy[0]+PD_BOND_CUTOFF/2, anchor_xy[1]-0.25,
            '2.60 Å', fontsize=9, color='red', ha='center', fontweight='bold')

    # Set view window: ±5 Å around anchor
    ax.set_xlim(anchor_xy[0]-5, anchor_xy[0]+5)
    ax.set_ylim(anchor_xy[1]-5, anchor_xy[1]+5)
    ax.set_aspect('equal')
    ax.set_title(label, fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.2)
    ax.set_xticks([]); ax.set_yticks([])


# ===== Build 4-panel methodology figure =====
fig, axes = plt.subplots(1, 4, figsize=(18, 5.2))

for i, (rel, lbl, kind) in enumerate(CASES):
    f = V2 / rel
    atoms = read(f)
    anchor = find_anchor(atoms, kind)
    if anchor is None: continue

    # Count neighbors in cutoff
    syms = atoms.get_chemical_symbols()
    sub_idx = [j for j in range(len(atoms)) if syms[j] in ('Pd', 'O') and j != anchor]
    if kind == 'CO':
        # exclude the O bonded to anchor C
        bonded_o = [j for j in sub_idx if syms[j] == 'O' and
                    atoms.get_distance(j, anchor, mic=True) < 1.4]
        sub_idx = [j for j in sub_idx if j not in bonded_o]
    d_sub = atoms.get_distances(anchor, sub_idx, mic=True)
    n_in = sum(1 for d in d_sub if d < PD_BOND_CUTOFF)
    pd_in = sum(1 for j, d in zip(sub_idx, d_sub) if d < PD_BOND_CUTOFF and syms[j] == 'Pd')
    o_in = sum(1 for j, d in zip(sub_idx, d_sub) if d < PD_BOND_CUTOFF and syms[j] == 'O')

    panel_title = f'({chr(ord("a")+i)}) {lbl}\n→ {n_in} neighbors  (n_Pd={pd_in}, n_O={o_in})'
    render_top_view(axes[i], atoms, anchor, kind, panel_title)

# Legend (custom)
legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='gold',
           markeredgecolor='black', markersize=15, label='anchor (C of CO)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#3070b0',
           markeredgecolor='red', markersize=14, label='Pd within 2.6 Å (counted)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#e74c3c',
           markeredgecolor='red', markersize=11, label='O within 2.6 Å (counted)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#bbbbbb',
           markeredgecolor='black', markersize=12, label='Pd outside (ignored)'),
    Line2D([0],[0], color='red', linestyle='--', linewidth=2, label='2.60 Å cutoff'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=5,
           bbox_to_anchor=(0.5, -0.02), frameon=True, fontsize=10)

fig.suptitle('Site Classification Methodology — 2.60 Å neighbor counting',
             fontsize=14, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig(OUT / 'site_methodology.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {OUT/"site_methodology.png"}')


# ===== Second figure: classification rule table =====
fig, ax = plt.subplots(figsize=(11, 7))
ax.axis('off')

rules = [
    ['n_neighbors', 'Composition', 'Label', 'Geometry'],
    ['0',           '—',           'physisorbed',   'isolated, no contact'],
    ['1',           'n_Pd=1',      'atop_Pd',       'on top of single Pd'],
    ['1',           'n_O=1',       'atop_O',        'on top of single O (lattice)'],
    ['2',           'n_Pd=2',      'bridge_Pd-Pd',  'between two Pd atoms'],
    ['2',           'n_Pd=1, n_O=1','bridge_Pd-O',  'between Pd and lattice O'],
    ['2',           'n_O=2',       'bridge_O-O',    'between two lattice O'],
    ['3',           'n_Pd=3',      'hollow_3Pd',    '3-fold Pd pocket (fcc/hcp)'],
    ['3',           'n_O=3',       'hollow_3O',     '3-fold O pocket'],
    ['3',           'mixed',       'hollow_3(xPd yO)', 'mixed 3-fold'],
    ['4+',          'any',         '4f / 5f (xPd yO)','high-coord (oxide cus)'],
]
table = ax.table(cellText=rules, colWidths=[0.12, 0.22, 0.25, 0.41],
                  cellLoc='left', loc='center')
table.auto_set_font_size(False); table.set_fontsize(11); table.scale(1, 2)
# Style header row
for j in range(4):
    table[(0, j)].set_facecolor('#1f4e79')
    table[(0, j)].set_text_props(color='white', fontweight='bold')

# Title
ax.set_title('Site Classification Rule Table  (anchor = C of CO or O of methoxy)',
             fontsize=13, fontweight='bold', pad=20)

plt.savefig(OUT / 'site_rule_table.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {OUT/"site_rule_table.png"}')
