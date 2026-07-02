"""F01_binding: oxidation trend in BINDING ENERGY (not d_min).

Replaces F01 d_min with proper MACE+D3 binding energy.
- Top-1 (most stable) per surface × adsorbate
- Median also shown (lighter marker) for fluke robustness
- NO LEGEND — labels placed on lines directly
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
plt.rcParams['axes.linewidth'] = 1.2

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
OUT = ROOT / 'reports/predft_advisor_figures'

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
order = ['S1', 'S2', 'S3b', 'S3', 'S4']

refs = json.load(open(G3 / 'mace_d3_references.json'))
E_SLAB = refs['slab']
E_CO = refs['gas']['CO']
E_CH3O = refs['gas']['CH3O_ref']

def load_unique(sid, ads):
    f = G3 / SDIRS[sid] / 'MLIP_phase1' / f'unique_{ads}.json'
    recs = json.load(open(f)) if f.exists() else []
    return sorted(recs, key=lambda r: r['E'])   # ascending E (top-1 first)

# Compute top-1 and median binding energies
co_top1, co_med = [], []
ch_top1, ch_med = [], []
for sid in order:
    co = load_unique(sid, 'CO')
    ch = load_unique(sid, 'CH3O')
    co_E = [r['E'] - E_SLAB[sid] - E_CO for r in co]
    ch_E = [r['E'] - E_SLAB[sid] - E_CH3O for r in ch]
    co_top1.append(min(co_E))
    co_med.append(float(np.median(co_E)))
    ch_top1.append(min(ch_E))
    ch_med.append(float(np.median(ch_E)))

fig, ax = plt.subplots(figsize=(8.5, 5.8))
xpos = list(range(len(order)))

# Top-1 (primary): solid lines, big markers
l_co, = ax.plot(xpos, co_top1, 'o-', markersize=13, lw=2.2, color='#1f4e79',
                markeredgecolor='black', zorder=4)
l_ch, = ax.plot(xpos, ch_top1, 's-', markersize=13, lw=2.2, color='#e76f51',
                markeredgecolor='black', zorder=4)

# Median (secondary): dashed, smaller, semi-transparent
ax.plot(xpos, co_med, 'o--', markersize=9, lw=1.2, color='#1f4e79',
        alpha=0.5, markerfacecolor='white', markeredgecolor='#1f4e79', zorder=3)
ax.plot(xpos, ch_med, 's--', markersize=9, lw=1.2, color='#e76f51',
        alpha=0.5, markerfacecolor='white', markeredgecolor='#e76f51', zorder=3)

# Annotations: top-1 values
for x, y in zip(xpos, co_top1):
    ax.annotate(f'{y:+.2f}', (x, y), xytext=(7, 8), textcoords='offset points',
                fontsize=11, color='#1f4e79', fontweight='bold')
for x, y in zip(xpos, ch_top1):
    ax.annotate(f'{y:+.2f}', (x, y), xytext=(7, -16), textcoords='offset points',
                fontsize=11, color='#e76f51', fontweight='bold')

# Reference shading
ax.axhline(0, ls='--', color='gray', alpha=0.6, lw=1.2)
ax.axhspan(-100, 0, alpha=0.07, color='green')   # bound region
ax.axhspan(0, 100, alpha=0.07, color='red')      # unbound region

# Direct labels on lines (no legend)
ax.text(4.15, co_top1[-1], r'CO$^*$', color='#1f4e79', fontsize=13,
        fontweight='bold', va='center')
ax.text(4.15, ch_top1[-1], r'CH$_3$O$^*$', color='#e76f51', fontsize=13,
        fontweight='bold', va='center')

# Text annotations for bound/unbound regions
ax.text(-0.42, 0.4, 'unbound', color='red', alpha=0.6, fontsize=11,
        rotation=90, va='center', fontweight='bold')
ax.text(-0.42, -0.4, 'bound', color='green', alpha=0.6, fontsize=11,
        rotation=90, va='center', fontweight='bold')

# Note about fluke
ax.text(4.0, -3.4, 'S4: MLIP fluke\n(decomposition)', color='gray', fontsize=9,
        ha='center', style='italic', alpha=0.7)

ax.set_xticks(xpos)
ax.set_xticklabels(['S1\n'+r'Pd$^0$', 'S2\n'+r'Pd$^0$+Pd$^{2+}$',
                    'S3b\n'+r'Pd$^{2+}$ Pd-top', 'S3\n'+r'Pd$^{2+}$ O-top',
                    'S4\n'+r'Pd$^{4+}$'])
ax.set_ylabel(r'top-1 $E_{\mathrm{bind}}$ / eV  (MACE+D3)')
ax.set_xlim(-0.5, 4.6)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT/'F01_oxidation_binding.png', dpi=160, bbox_inches='tight')
plt.close()

print(f'Saved: {OUT/"F01_oxidation_binding.png"}')
print('\nTop-1 binding energy / eV:')
print(f"{'Surface':<8} {'CO':<10} {'CH3O':<10}")
for sid, c, h in zip(order, co_top1, ch_top1):
    print(f"  {sid:<8} {c:+.3f}     {h:+.3f}")
print('\nMedian binding energy / eV:')
print(f"{'Surface':<8} {'CO':<10} {'CH3O':<10}")
for sid, c, h in zip(order, co_med, ch_med):
    print(f"  {sid:<8} {c:+.3f}     {h:+.3f}")
