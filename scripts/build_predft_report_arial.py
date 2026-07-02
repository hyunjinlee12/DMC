"""Pre-DFT advisor report — Arial-style figures (no titles, axes/ticks/legend only).

Figures (Liberation Sans = Arial metric-compatible):
  F01_oxidation_trend.png       — Pd⁰→Pd⁴⁺ chemisorption strength
  F02_descriptor_preview.png    — E_CO vs E_CH3O scatter
  F03_d_min_distributions.png   — Pd-C, Pd-O distributions
  F04_d_reactive_coads.png      — co-ads distance per surface
  F05_site_distribution.png     — v2 site type bar
  F06_representative_S1.png     — S1 representative structures
  F07_representative_S2.png
  F08_representative_S3.png
  F09_representative_S3b.png
  F10_representative_S4.png
  F11_committee_timeline.png    — 11 committee cycles

Style: no titles, Arial-equiv (Liberation Sans), axes/labels/ticks/legend only.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from ase.io import read
from ase.visualize.plot import plot_atoms
from collections import Counter

# === ARIAL-EQUIV STYLE (Liberation Sans, Arial metric-compatible) ===
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.titlesize'] = 0   # never use titles
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab'
G3 = ROOT / 'calculations/G3_adsorption'
SHORTLIST = G3 / 'DFT_shortlist_v2'
OUT_FIG = ROOT / 'reports/predft_advisor_figures'
OUT_FIG.mkdir(parents=True, exist_ok=True)

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
COLORS = {'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}
E_SLAB = {'S1':-434.380,'S2':-618.565,'S3':-724.103,'S3b':-570.772,'S4':-788.493}


# Load Phase 1 unique data
def load_phase1():
    data = {}
    for sid in SURFACES:
        d = G3 / SDIRS[sid] / 'MLIP_phase1'
        data[sid] = {
            'CO': json.load(open(d/'unique_CO.json')),
            'CH3O': json.load(open(d/'unique_CH3O.json')),
        }
    return data


def load_phase2():
    data = {}
    for sid in SURFACES:
        f = G3 / SDIRS[sid] / 'MLIP_phase2_filtered/unique_SetA.json'
        data[sid] = json.load(open(f)) if f.exists() else []
    return data


p1 = load_phase1()
p2 = load_phase2()

# ==============================
# F01: Oxidation trend (top-1 d_min for CO and CH3O)
# ==============================
fig, ax = plt.subplots(figsize=(8, 5.5))
order = ['S1', 'S2', 'S3b', 'S3', 'S4']
xpos = list(range(len(order)))
co_d = [p1[s]['CO'][0]['d_min'] for s in order]
ch_d = [p1[s]['CH3O'][0]['d_min'] for s in order]
ax.plot(xpos, co_d, 'o-', markersize=14, lw=2.0, color='#1f4e79',
        markeredgecolor='black', label=r'CO$^*$ (Pd–C)')
ax.plot(xpos, ch_d, 's-', markersize=14, lw=2.0, color='#e76f51',
        markeredgecolor='black', label=r'CH$_3$O$^*$ (Pd–O)')
for x, y in zip(xpos, co_d):
    ax.annotate(f'{y:.2f}', (x, y), xytext=(6, 6), textcoords='offset points', fontsize=11)
for x, y in zip(xpos, ch_d):
    ax.annotate(f'{y:.2f}', (x, y), xytext=(6, -14), textcoords='offset points', fontsize=11)
ax.axhspan(1.85, 2.15, alpha=0.12, color='green', label='chemisorbed band')
ax.axhline(3.0, ls='--', color='red', alpha=0.4, lw=1.2, label='physisorption (~3 Å)')
ax.set_xticks(xpos)
ax.set_xticklabels(['S1\n'+r'Pd$^0$', 'S2\n'+r'Pd$^0$+Pd$^{2+}$',
                    'S3b\n'+r'Pd$^{2+}$ Pd-top', 'S3\n'+r'Pd$^{2+}$ O-top',
                    'S4\n'+r'Pd$^{4+}$'])
ax.set_ylabel(r'top-1 $d_{\mathrm{min}}$ / Å')
ax.legend(loc='upper left', frameon=True)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_FIG/'F01_oxidation_trend.png', dpi=160, bbox_inches='tight')
plt.close()
print('F01 ✓')

# ==============================
# F02: Descriptor map preview
# ==============================
fig, ax = plt.subplots(figsize=(8.5, 7))
for sid in SURFACES:
    co_top = p1[sid]['CO'][0]['E']
    ch_top = p1[sid]['CH3O'][0]['E']
    dE_CO = co_top - E_SLAB[sid]
    dE_CH = ch_top - E_SLAB[sid]
    ax.scatter(dE_CO, dE_CH, s=350, color=COLORS[sid], edgecolor='black',
               linewidth=2, zorder=3, label=sid)
    ax.annotate(f'  {sid}', (dE_CO, dE_CH), xytext=(10, 5),
                textcoords='offset points', fontsize=13, fontweight='bold')
ax.set_xlabel(r'$E_{\mathrm{MACE+D3}}^{\mathrm{slab+CO}} - E_{\mathrm{slab}}$  /  eV')
ax.set_ylabel(r'$E_{\mathrm{MACE+D3}}^{\mathrm{slab+CH_3O}} - E_{\mathrm{slab}}$  /  eV')
ax.invert_xaxis(); ax.invert_yaxis()
ax.grid(True, alpha=0.3)
ax.legend(loc='upper left', ncol=1, frameon=True)
plt.tight_layout()
plt.savefig(OUT_FIG/'F02_descriptor_preview.png', dpi=160, bbox_inches='tight')
plt.close()
print('F02 ✓')

# ==============================
# F03: d_min distributions (per surface)
# ==============================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
for ax, ads in [(axes[0], 'CO'), (axes[1], 'CH3O')]:
    for sid in SURFACES:
        recs = p1[sid][ads]
        d_mins = [r['d_min'] for r in recs]
        ax.hist(d_mins, bins=25, alpha=0.55, color=COLORS[sid],
                edgecolor='black', label=f'{sid} (n={len(d_mins)})', linewidth=0.6)
    band = (1.85, 2.10) if ads == 'CO' else (2.00, 2.15)
    ax.axvspan(*band, alpha=0.13, color='green', label='chemisorbed')
    ax.axvline(3.0, ls='--', color='red', alpha=0.4, lw=1.2, label='physisorption')
    ax.set_xlabel(r'$d_{\mathrm{min}}$ / Å')
    ax.set_ylabel('count' if ads == 'CO' else '')
    ax.legend(loc='upper right', fontsize=10, frameon=True)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1.5, 5.0)
plt.tight_layout()
plt.savefig(OUT_FIG/'F03_d_min_distributions.png', dpi=160, bbox_inches='tight')
plt.close()
print('F03 ✓')

# ==============================
# F04: d_reactive (co-ads) per surface
# ==============================
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()
for i, sid in enumerate(SURFACES):
    ax = axes[i]
    recs = p2.get(sid, [])
    if not recs:
        ax.text(0.5, 0.5, f'{sid}\n(no data)', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_xticks([]); ax.set_yticks([])
        continue
    E0 = recs[0]['E']
    d = [r['d_reactive'] for r in recs]
    dE = [(r['E'] - E0)*1000 for r in recs]
    ax.scatter(d, dE, s=14, alpha=0.4, color=COLORS[sid])
    ax.axvspan(2.1, 4.0, alpha=0.12, color='green', label='SetA band')
    ax.set_xlabel(r'$d(\mathrm{C}_{\mathrm{CO}} - \mathrm{O}_{\mathrm{CH_3O}})$ / Å')
    ax.set_ylabel(r'$\Delta E_{\mathrm{MACE}}$ / meV')
    ax.legend(loc='upper right', fontsize=10, frameon=True)
    ax.grid(True, alpha=0.3)
    # surface label as text inside panel (no title)
    ax.text(0.04, 0.94, f'{sid} (n={len(recs)})', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='gray'))
axes[-1].axis('off')
plt.tight_layout()
plt.savefig(OUT_FIG/'F04_d_reactive_coads.png', dpi=160, bbox_inches='tight')
plt.close()
print('F04 ✓')

# ==============================
# F05: v2 site distribution per surface (stacked bar)
# ==============================
# Collect from v2 POSCARs (filename has site info)
records = []
for vasp in SHORTLIST.rglob('*.vasp'):
    sid = vasp.parts[-3]
    kind = vasp.parts[-2]
    name = vasp.stem
    # parse site
    if 'physisorbed' in name: site = 'physisorbed'
    elif 'hollow_3Pd' in name: site = 'hollow_3Pd'
    elif 'hollow_3' in name: site = 'hollow_3(Pd-O)'
    elif '4f' in name: site = '4-fold'
    elif 'bridge_Pd-Pd' in name: site = 'bridge_Pd-Pd'
    elif 'bridge_Pd-O' in name: site = 'bridge_Pd-O'
    elif 'bridge_O-O' in name: site = 'bridge_O-O'
    elif 'atop_Pd' in name and 'OMe' not in name: site = 'atop_Pd'
    elif 'atop_O' in name: site = 'atop_O'
    else: site = 'other'
    records.append({'surface': sid, 'kind': kind, 'site': site})

all_sites = sorted(set(r['site'] for r in records))
cmap = plt.cm.tab20
site_color = {s: cmap(i / max(len(all_sites), 1)) for i, s in enumerate(all_sites)}

fig, axes = plt.subplots(1, 3, figsize=(17, 6))
kinds = ['single_CO', 'single_CH3O', 'coads_SetA']
labels_kind = [r'(a) CO$^*$', r'(b) CH$_3$O$^*$', r'(c) co-ads (both anchors)']
for ax, kind, lbl in zip(axes, kinds, labels_kind):
    surf_data = {sid: Counter() for sid in SURFACES}
    for r in records:
        if r['kind'] != kind: continue
        surf_data[r['surface']][r['site']] += 1
    x = np.arange(len(SURFACES))
    bottoms = np.zeros(len(SURFACES))
    for site in all_sites:
        vals = [surf_data[sid].get(site, 0) for sid in SURFACES]
        if sum(vals) == 0: continue
        ax.bar(x, vals, bottom=bottoms, color=site_color[site], edgecolor='black',
               label=site, linewidth=0.5)
        bottoms += vals
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('count' if kind == 'single_CO' else '')
    ax.text(0.04, 0.94, lbl, transform=ax.transAxes, fontsize=12, fontweight='bold',
            va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
    ax.legend(loc='upper right', fontsize=8, ncol=2, frameon=True)
    ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(OUT_FIG/'F05_site_distribution.png', dpi=160, bbox_inches='tight')
plt.close()
print('F05 ✓')

# ==============================
# F06-F10: Representative structures per surface (top + side)
# ==============================
def render_surface(sid, out_path):
    sid_dir = SHORTLIST / sid
    if not sid_dir.exists(): return
    # Collect all candidates per kind
    candidates = {}
    for kind in ['single_CO', 'single_CH3O', 'coads_SetA']:
        kd = sid_dir / kind
        if kd.exists():
            candidates[kind] = sorted(kd.glob('*.vasp'))

    # Layout: 3 columns (kinds), N rows (max candidates per kind × 2 for top+side)
    max_per_kind = max(len(v) for v in candidates.values()) if candidates else 0
    if max_per_kind == 0: return
    n_cols = sum(1 for v in candidates.values() if v)
    fig = plt.figure(figsize=(5.5 * n_cols, 3.0 * max_per_kind))
    gs = GridSpec(max_per_kind, n_cols * 2, figure=fig, hspace=0.30, wspace=0.05)

    col = 0
    kind_labels = {'single_CO': r'CO$^*$', 'single_CH3O': r'CH$_3$O$^*$', 'coads_SetA': r'CO$^*$ + CH$_3$O$^*$'}
    for kind in ['single_CO', 'single_CH3O', 'coads_SetA']:
        vasps = candidates.get(kind, [])
        if not vasps: continue
        for row, vasp in enumerate(vasps):
            atoms = read(vasp)
            # top view
            ax_top = fig.add_subplot(gs[row, col*2])
            plot_atoms(atoms, ax_top, rotation='0x,0y,0z', radii=0.85, show_unit_cell=2)
            ax_top.set_xticks([]); ax_top.set_yticks([])
            # side view
            ax_side = fig.add_subplot(gs[row, col*2 + 1])
            plot_atoms(atoms, ax_side, rotation='-80x,5y,0z', radii=0.85, show_unit_cell=2)
            ax_side.set_xticks([]); ax_side.set_yticks([])

            # parse short info
            name = vasp.stem
            # extract site label
            for marker in ['atop_Pd', 'atop_O', 'bridge_Pd-Pd', 'bridge_Pd-O',
                           'bridge_O-O', 'hollow_3Pd', 'hollow_3', '4f', 'physisorbed']:
                if marker in name:
                    site = marker; break
            else:
                site = '?'
            # extract rank
            rank = name.split('_')[0]
            ax_top.text(0.04, 0.96, f'rank {rank}: {site}', transform=ax_top.transAxes,
                        fontsize=10, va='top',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
            if row == 0:
                ax_top.text(0.5, 1.08, kind_labels[kind], transform=ax_top.transAxes,
                            fontsize=12, fontweight='bold', ha='center')
        col += 1

    plt.savefig(out_path, dpi=130, bbox_inches='tight')
    plt.close()
    print(f'  rendered {sid} → {out_path.name}')


for sid in SURFACES:
    render_surface(sid, OUT_FIG / f'F0{6 + SURFACES.index(sid)}_representative_{sid}.png')

# ==============================
# F11: Committee timeline (already have, just regen with Arial style)
# ==============================
fig, ax = plt.subplots(figsize=(13, 6))
cycles = ['Phase1\n(no-D3)', 'Phase1\n(D3)', 'Phase2\n(raw)', 'Phase2\n(filt)',
          'Phase3\n(raw)', 'Phase3\n(filt)', 'G1\nre-audit', 'G2\nre-audit',
          'G2\n(lit)', 'T1.10-15\naudit', 'T1.15 v2\nsite-strict']
judges = ['methods', 'physics', 'statistics', 'silent-error', 'malicious']
# 4=Pass, 3=Pass-w-c, 2=Concern, 1=Reject, None=N/A
matrix = [
    [3, 4, 2, 1, None],
    [4, 2, 2, 4, None],
    [4, 1, 2, 1, 4],
    [4, 2, 2, 4, 4],
    [4, 1, 2, 1, 4],
    [4, 2, 2, 4, 4],
    [4, 2, 4, 2, 4],
    [None, 4, 4, 2, 4],
    [4, 4, None, None, None],   # lit check
    [4, 1, 3, 2, 4],   # 10
    [4, 4, 4, 4, 4],   # 11 (v2 strict)
]
cmap_v = {1:'#c0392b', 2:'#f39c12', 3:'#27ae60', 4:'#16a085', None:'#cccccc'}
labels_short = {1:'R', 2:'C', 3:'P-c', 4:'P', None:'—'}
labels_full = {4:'Pass', 3:'Pass-w-c', 2:'Concern', 1:'Reject', None:'N/A'}
for i, cycle in enumerate(cycles):
    for j, jud in enumerate(judges):
        v = matrix[i][j] if j < len(matrix[i]) else None
        ax.add_patch(plt.Rectangle((j, len(cycles)-1-i), 1, 1, facecolor=cmap_v[v], edgecolor='black'))
        ax.text(j+0.5, len(cycles)-1-i+0.5, labels_short[v], ha='center', va='center',
                fontsize=11, fontweight='bold', color='white')
ax.set_xlim(0, len(judges)); ax.set_ylim(0, len(cycles))
ax.set_xticks([i+0.5 for i in range(len(judges))]); ax.set_xticklabels(judges, rotation=15)
ax.set_yticks([i+0.5 for i in range(len(cycles))]); ax.set_yticklabels(reversed(cycles))
ax.tick_params(axis='both', length=0)   # remove tick marks since this is a heatmap
handles = [mpatches.Patch(color=cmap_v[v], label=labels_full[v]) for v in [4,3,2,1,None]]
ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=True)
plt.tight_layout()
plt.savefig(OUT_FIG/'F11_committee_timeline.png', dpi=160, bbox_inches='tight')
plt.close()
print('F11 ✓')

print(f'\nAll figures saved to {OUT_FIG}')
