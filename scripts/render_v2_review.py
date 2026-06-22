"""Render all 48 v2 DFT candidates for visual review.

Outputs:
  reports/G3/v2_review/
    A_v2_site_distribution.png       — site distribution
    {surface}_single_CO.png          — per-surface CO grid (5 panels)
    {surface}_single_CH3O.png        — CH3O grid
    {surface}_coads.png              — coads grid
    summary_v2.md                    — text summary
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ase.io import read
from ase.visualize.plot import plot_atoms
from collections import Counter

plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10})

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
SHORTLIST_V2 = ROOT / 'calculations/G3_adsorption/DFT_shortlist_v2'
OUT = ROOT / 'reports/G3/v2_review'
OUT.mkdir(parents=True, exist_ok=True)

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']


def render_grid(vasps, title, out_path, ncols=3):
    """Render a grid of side+top views for a list of POSCAR paths."""
    n = len(vasps)
    if n == 0:
        return
    nrows = (n + ncols - 1) // ncols
    fig = plt.figure(figsize=(5 * ncols, 6 * nrows))
    gs = GridSpec(nrows * 2, ncols, figure=fig, hspace=0.3, wspace=0.1)

    for i, vasp in enumerate(vasps):
        atoms = read(vasp)
        row = i // ncols
        col = i % ncols
        # Extract info from filename
        name = vasp.stem
        # parse: e.g., "00_single_CO_atop_Pd_idx00064"
        parts = name.split('_')
        rank = parts[0]
        # E from info?
        # Top view
        ax_top = fig.add_subplot(gs[row * 2, col])
        plot_atoms(atoms, ax_top, rotation='0x,0y,0z', radii=0.85, show_unit_cell=2)
        ax_top.set_xticks([]); ax_top.set_yticks([])
        ax_top.set_title(f'rank {rank}: {name[:55]}', fontsize=9)
        # Side view
        ax_side = fig.add_subplot(gs[row * 2 + 1, col])
        plot_atoms(atoms, ax_side, rotation='-80x,5y,0z', radii=0.85, show_unit_cell=2)
        ax_side.set_xticks([]); ax_side.set_yticks([])

    fig.suptitle(title, fontsize=13, weight='bold', y=0.995)
    plt.savefig(out_path, dpi=110, bbox_inches='tight')
    plt.close()


for sid in SURFACES:
    sid_dir = SHORTLIST_V2 / sid
    if not sid_dir.exists(): continue
    for kind in ['single_CO', 'single_CH3O', 'coads_SetA']:
        kdir = sid_dir / kind
        if not kdir.exists(): continue
        vasps = sorted(kdir.glob('*.vasp'))
        if not vasps: continue
        title = f'{sid} — {kind}  (n={len(vasps)})'
        out = OUT / f'{sid}_{kind}.png'
        render_grid(vasps, title, out, ncols=min(len(vasps), 3))
        print(f'rendered {sid}/{kind}: {len(vasps)} → {out.name}')

# Site distribution figure (re-use stats already there)
print('\nBuilding site distribution figure...')
records = []
for vasp in SHORTLIST_V2.rglob('*.vasp'):
    sid = vasp.parts[-3]
    kind = vasp.parts[-2]
    # parse site from filename
    name = vasp.stem
    # Find site labels (atop_Pd, atop_O, bridge_Pd-Pd, hollow_3Pd, etc.)
    if 'physisorbed' in name: site = 'physisorbed'
    elif 'atop_Pd' in name and 'OMe' not in name: site = 'atop_Pd'
    elif 'atop_O' in name: site = 'atop_O'
    elif 'bridge_Pd-Pd' in name: site = 'bridge_Pd-Pd'
    elif 'bridge_Pd-O' in name: site = 'bridge_Pd-O'
    elif 'bridge_O-O' in name: site = 'bridge_O-O'
    elif 'hollow_3Pd' in name: site = 'hollow_3Pd'
    elif '4f' in name: site = '4-fold'
    else: site = 'other'
    records.append({'surface': sid, 'kind': kind, 'site': site, 'name': name})

# Plot
fig, axes = plt.subplots(1, 3, figsize=(17, 6))
titles = ['(a) CO* anchor', '(b) CH₃O* anchor', '(c) co-ads (both anchors)']
kinds = ['single_CO', 'single_CH3O', 'coads_SetA']
all_sites = sorted(set(r['site'] for r in records))
cmap = plt.cm.tab20
site_color = {s: cmap(i / max(len(all_sites), 1)) for i, s in enumerate(all_sites)}

for ax, kind, title in zip(axes, kinds, titles):
    surf_data = {sid: Counter() for sid in SURFACES}
    for r in records:
        if r['kind'] != kind: continue
        surf_data[r['surface']][r['site']] += 1
    x = np.arange(len(SURFACES))
    bottoms = np.zeros(len(SURFACES))
    for site in all_sites:
        vals = [surf_data[sid].get(site, 0) for sid in SURFACES]
        if sum(vals) == 0: continue
        ax.bar(x, vals, bottom=bottoms, color=site_color[site], edgecolor='black', label=site, linewidth=0.5)
        bottoms += vals
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('Count'); ax.set_title(title)
    ax.legend(fontsize=8, loc='upper right', ncol=2)
    ax.grid(True, alpha=0.3, axis='y')
plt.suptitle('v2 DFT Shortlist Site Distribution (48 candidates, guide-strict)',
             fontsize=13, weight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUT / 'A_v2_site_distribution.png', dpi=140, bbox_inches='tight')
plt.close()
print(f'\nAll figures in {OUT}')

# Summary text
summary = OUT / 'v2_summary.md'
summary.write_text(
    "# v2 DFT Shortlist 시각 검토 자료\n\n"
    "48 candidates (guide-strict re-pick).\n\n"
    "## 파일\n"
    "- A_v2_site_distribution.png : site type 분포 (3 panels)\n"
    "- {S1,S2,S3,S3b,S4}_single_CO.png : 표면별 CO* 후보 grid (top+side views)\n"
    "- {S1,S2,S3,S3b,S4}_single_CH3O.png : CH3O* 후보\n"
    "- {S1,S2,S3,S3b,S4}_coads_SetA.png : co-ads 후보\n"
)
print(f'Summary: {summary}')
