"""Build F12_binding_distribution.png — per-surface binding energy distribution.

Replaces / complements F01 (top-1 d_min only) with:
  - box plot + violin + scatter overlay (all unique candidates)
  - per surface (S1, S2, S3, S3b, S4) × per adsorbate (CO, CH3O)
  - binding proxy = E(slab+ads) - E(slab)  (relative, MACE+D3)

Note: MLIP absolute E 는 DFT 와 다름. distribution shape + relative ordering 이 의미.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['legend.fontsize'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
OUT = ROOT / 'reports/predft_advisor_figures'

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
COLORS = {'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}
E_SLAB = {'S1':-434.380,'S2':-618.565,'S3':-724.103,'S3b':-570.772,'S4':-788.493}


def load_unique(sid, ads):
    f = G3 / SDIRS[sid] / 'MLIP_phase1' / f'unique_{ads}.json'
    return json.load(open(f)) if f.exists() else []


def load_unique_coads(sid):
    f = G3 / SDIRS[sid] / 'MLIP_phase2_filtered' / 'unique_SetA.json'
    return json.load(open(f)) if f.exists() else []


# Collect E_binding (relative) for each surface × adsorbate
fig, axes = plt.subplots(1, 3, figsize=(18, 6.5))

# === (a) CO* ===
ax = axes[0]
data_CO = []
for sid in SURFACES:
    recs = load_unique(sid, 'CO')
    if recs:
        E_top1 = min(r['E'] for r in recs)
        E_relmev = [(r['E'] - E_top1)*1000 for r in recs]
    else:
        E_relmev = []
    data_CO.append(E_relmev)

# Box + violin + scatter
positions = np.arange(len(SURFACES))
parts = ax.violinplot(data_CO, positions=positions, widths=0.7, showmeans=False,
                     showmedians=False, showextrema=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(COLORS[SURFACES[i]]); pc.set_alpha(0.4); pc.set_edgecolor('black')

bp = ax.boxplot(data_CO, positions=positions, widths=0.3, patch_artist=True,
                showfliers=False, medianprops={'color':'red','lw':2})
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(COLORS[SURFACES[i]]); patch.set_alpha(0.7)

# Individual points overlay (jittered)
np.random.seed(42)
for i, d in enumerate(data_CO):
    if not d: continue
    x_jit = np.random.normal(i, 0.04, size=len(d))
    ax.scatter(x_jit, d, s=4, color='black', alpha=0.3, zorder=2)
    # mark top-1 (most stable = most negative)
    # top-1 = 0 by construction (we subtracted E_top1). highlight as gold star at y=0
    ax.scatter([i], [0], s=140, color='gold', edgecolor='black', linewidth=2, zorder=5,
               marker='*', label='top-1 (lowest E)' if i == 0 else None)

ax.set_xticks(positions); ax.set_xticklabels(SURFACES)
ax.set_ylabel(r'$\Delta E$ vs top-1 / meV')
ax.text(0.04, 0.96, r'(a) CO$^*$', transform=ax.transAxes, fontsize=13, fontweight='bold',
        va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
ax.legend(loc='lower right', frameon=True)
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=-50)  # top-1 at 0, rest above

# === (b) CH3O* ===
ax = axes[1]
data_CH = []
for sid in SURFACES:
    recs = load_unique(sid, 'CH3O')
    if recs:
        E_top1 = min(r['E'] for r in recs)
        E_relmev = [(r['E'] - E_top1)*1000 for r in recs]
    else:
        E_relmev = []
    data_CH.append(E_relmev)

parts = ax.violinplot(data_CH, positions=positions, widths=0.7, showmeans=False,
                     showmedians=False, showextrema=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(COLORS[SURFACES[i]]); pc.set_alpha(0.4); pc.set_edgecolor('black')

bp = ax.boxplot(data_CH, positions=positions, widths=0.3, patch_artist=True,
                showfliers=False, medianprops={'color':'red','lw':2})
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(COLORS[SURFACES[i]]); patch.set_alpha(0.7)

for i, d in enumerate(data_CH):
    if not d: continue
    x_jit = np.random.normal(i, 0.04, size=len(d))
    ax.scatter(x_jit, d, s=4, color='black', alpha=0.3, zorder=2)
    # top-1 = 0 by construction (we subtracted E_top1). highlight as gold star at y=0
    ax.scatter([i], [0], s=140, color='gold', edgecolor='black', linewidth=2, zorder=5,
               marker='*', label='top-1 (lowest E)' if i == 0 else None)

ax.set_xticks(positions); ax.set_xticklabels(SURFACES)
ax.set_ylabel(r'$E_{\mathrm{MACE+D3}}^{\mathrm{slab+CH_3O^*}} - E_{\mathrm{slab}}$ / eV')
ax.text(0.04, 0.96, r'(b) CH$_3$O$^*$', transform=ax.transAxes, fontsize=13, fontweight='bold',
        va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
ax.legend(loc='lower right', frameon=True)
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=-50)

# === (c) co-ads SetA ===
ax = axes[2]
data_coads = []
for sid in SURFACES:
    recs = load_unique_coads(sid)
    if recs:
        E_top1 = min(r['E'] for r in recs)
        E_relmev = [(r['E'] - E_top1)*1000 for r in recs]
    else:
        E_relmev = []
    data_coads.append(E_relmev)

parts = ax.violinplot(data_coads, positions=positions, widths=0.7, showmeans=False,
                     showmedians=False, showextrema=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(COLORS[SURFACES[i]]); pc.set_alpha(0.4); pc.set_edgecolor('black')

bp = ax.boxplot(data_coads, positions=positions, widths=0.3, patch_artist=True,
                showfliers=False, medianprops={'color':'red','lw':2})
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(COLORS[SURFACES[i]]); patch.set_alpha(0.7)

for i, d in enumerate(data_coads):
    if not d: continue
    x_jit = np.random.normal(i, 0.04, size=len(d))
    ax.scatter(x_jit, d, s=4, color='black', alpha=0.2, zorder=2)
    # top-1 = 0 by construction (we subtracted E_top1). highlight as gold star at y=0
    ax.scatter([i], [0], s=140, color='gold', edgecolor='black', linewidth=2, zorder=5,
               marker='*', label='top-1 (lowest E)' if i == 0 else None)

ax.set_xticks(positions); ax.set_xticklabels(SURFACES)
ax.set_ylabel(r'$E_{\mathrm{MACE+D3}}^{\mathrm{slab+CO^*+CH_3O^*}} - E_{\mathrm{slab}}$ / eV')
ax.text(0.04, 0.96, r'(c) co-ads (SetA filtered)', transform=ax.transAxes, fontsize=13,
        fontweight='bold', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
ax.legend(loc='lower right', frameon=True)
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=-50)

plt.tight_layout()
plt.savefig(OUT/'F12_binding_distribution.png', dpi=160, bbox_inches='tight')
plt.close()
print(f'F12 saved: {OUT/"F12_binding_distribution.png"}')

# Print summary stats
print('\n=== Summary: median / IQR / range per surface × ads ===')
print(f"{'Sur':<5} {'ads':<8} {'n':<5} {'median':<10} {'IQR':<15} {'top-1 (min)':<12} {'max':<10}")
for sid in SURFACES:
    for ads, dd in [('CO', data_CO), ('CH3O', data_CH), ('coads', data_coads)]:
        idx = SURFACES.index(sid)
        d = dd[idx]
        if not d: continue
        d = np.array(d)
        med = np.median(d)
        q1, q3 = np.percentile(d, [25, 75])
        print(f"{sid:<5} {ads:<8} {len(d):<5} {med:<10.3f} {q1:.3f}–{q3:.3f}    {min(d):<12.3f} {max(d):<10.3f}")
