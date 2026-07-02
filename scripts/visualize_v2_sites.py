"""Render v2 DFT shortlist picks — one row per surface, columns = candidates.
Shows top-down view with adsorbate highlighted. ASE plot, no labels on figure
(filename label shown in title only)."""
import glob, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read
from ase.visualize.plot import plot_atoms

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
V2 = ROOT / 'calculations/G3_adsorption/DFT_shortlist_v2'
OUT = ROOT / 'reports/predft_advisor_figures/v2_site_renders'
OUT.mkdir(parents=True, exist_ok=True)

SURF = ['S1', 'S2', 'S3', 'S3b', 'S4']


def parse_label(fname):
    name = Path(fname).stem
    m = re.match(r'\d+_single_(CO|CH3O)_(.+)_idx\d+', name)
    if m: return m.group(1), m.group(2)
    m = re.match(r'\d+_coads_CO-(.+?)_OMe-(.+?)_d', name)
    if m: return 'coads', f'CO[{m.group(1)}]/OMe[{m.group(2)}]'
    return None, None


def render_one(ax, atoms, label, view='top'):
    if view == 'top':
        rot = '0x,0y,0z'
    else:
        rot = '-90x'
    plot_atoms(atoms, ax, rotation=rot, radii=0.55, show_unit_cell=2)
    ax.set_axis_off()
    ax.set_title(label, fontsize=9, pad=3)


# ===== single CO panels (one figure per surface) =====
for sid in SURF:
    co_files = sorted(glob.glob(str(V2 / sid / 'single_CO' / '*.vasp')))
    ch_files = sorted(glob.glob(str(V2 / sid / 'single_CH3O' / '*.vasp')))
    coads_files = sorted(glob.glob(str(V2 / sid / 'coads_SetA' / '*.vasp')))

    n_max = max(len(co_files), len(ch_files), len(coads_files))
    if n_max == 0: continue
    fig, axes = plt.subplots(3, n_max, figsize=(2.4*n_max, 7.5))
    if n_max == 1: axes = axes.reshape(-1, 1)

    for col, (files, row_label) in enumerate([(co_files, 'CO'),
                                              (ch_files, 'CH3O'),
                                              (coads_files, 'co-ads')]):
        pass  # placeholder, we iterate per row below

    for row, files in enumerate([co_files, ch_files, coads_files]):
        for col in range(n_max):
            ax = axes[row, col]
            if col < len(files):
                atoms = read(files[col])
                _, lbl = parse_label(files[col])
                render_one(ax, atoms, f'{lbl}' if lbl else '?', view='top')
            else:
                ax.set_axis_off()

    # Add row labels on first column
    row_titles = ['CO*', 'CH3O*', 'co-ads (SetA)']
    for row, title in enumerate(row_titles):
        axes[row, 0].text(-0.15, 0.5, title, transform=axes[row, 0].transAxes,
                          fontsize=12, fontweight='bold', va='center', ha='right',
                          rotation=90)

    fig.suptitle(f'{sid} DFT shortlist v2 — top view', fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0.02, 0.02, 1, 0.98])
    plt.savefig(OUT / f'{sid}_v2_sites.png', dpi=120, bbox_inches='tight')
    plt.close()
    print(f'Saved: {OUT/(sid+"_v2_sites.png")}')

# ===== Combined summary: one site type per surface (representative) =====
fig, axes = plt.subplots(3, 5, figsize=(15, 9))
for col, sid in enumerate(SURF):
    co_files = sorted(glob.glob(str(V2 / sid / 'single_CO' / '*.vasp')))
    ch_files = sorted(glob.glob(str(V2 / sid / 'single_CH3O' / '*.vasp')))
    coads_files = sorted(glob.glob(str(V2 / sid / 'coads_SetA' / '*.vasp')))

    for row, (files, kind) in enumerate([(co_files, 'CO'),
                                          (ch_files, 'CH3O'),
                                          (coads_files, 'coads')]):
        ax = axes[row, col]
        if files:
            atoms = read(files[0])   # rank 00 = lowest E
            _, lbl = parse_label(files[0])
            render_one(ax, atoms, lbl or '?', view='top')
        else:
            ax.set_axis_off()
            ax.text(0.5, 0.5, 'n/a', ha='center', va='center', fontsize=10,
                    color='gray', transform=ax.transAxes)

for col, sid in enumerate(SURF):
    axes[0, col].set_title(f'{sid}\n{axes[0, col].get_title()}', fontsize=11, fontweight='bold')

row_titles = ['CO*', 'CH3O*', 'co-ads (SetA)']
for row, title in enumerate(row_titles):
    axes[row, 0].text(-0.18, 0.5, title, transform=axes[row, 0].transAxes,
                      fontsize=13, fontweight='bold', va='center', ha='right',
                      rotation=90)

plt.tight_layout()
plt.savefig(OUT / 'ALL_top1_sites.png', dpi=120, bbox_inches='tight')
plt.close()
print(f'Saved: {OUT/"ALL_top1_sites.png"}')
