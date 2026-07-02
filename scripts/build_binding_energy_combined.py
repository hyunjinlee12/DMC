"""F14: 2 ref schemes stacked (top = thermo CH3OH-1/2H2, bottom = radical CH3O)."""
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

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
COLORS = {'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs = json.load(open(G3 / 'mace_d3_references.json'))
E_SLAB = refs['slab']
E_CO = refs['gas']['CO']
E_CH3O_thermo = refs['gas']['CH3O_ref']        # CH3OH - 1/2 H2
E_CH3O_radical = refs['gas']['CH3O_radical']   # direct CH3O

def load_uniq(sid, ads):
    f = G3 / SDIRS[sid] / 'MLIP_phase1' / f'unique_{ads}.json'
    return json.load(open(f)) if f.exists() else []
def load_coads(sid):
    f = G3 / SDIRS[sid] / 'MLIP_phase2_filtered' / 'unique_SetA.json'
    return json.load(open(f)) if f.exists() else []

def plot_panel(ax, data_list, letter, ads_label):
    positions = np.arange(len(SURFACES))
    parts = ax.violinplot(data_list, positions=positions, widths=0.7, showmeans=False,
                         showmedians=False, showextrema=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(COLORS[SURFACES[i]]); pc.set_alpha(0.4); pc.set_edgecolor('black')
    bp = ax.boxplot(data_list, positions=positions, widths=0.3, patch_artist=True,
                    showfliers=False, medianprops={'color':'red','lw':2})
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(COLORS[SURFACES[i]]); patch.set_alpha(0.7)
    np.random.seed(42)
    for i, d in enumerate(data_list):
        if not d: continue
        ax.scatter(np.random.normal(i, 0.04, size=len(d)), d, s=4, color='black', alpha=0.3, zorder=2)
        ax.scatter([i], [min(d)], s=140, color='gold', edgecolor='black', linewidth=2, zorder=5, marker='*')
    ax.set_xticks(positions); ax.set_xticklabels(SURFACES)
    ax.axhline(0, ls='--', color='gray', alpha=0.6, lw=1)
    ax.text(0.04, 0.96, f'({letter}) {ads_label}', transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
    ax.grid(True, alpha=0.3)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Top row: thermo ref (CH3OH - 1/2 H2) — previous
d_CO  = [[r['E']-E_SLAB[s]-E_CO for r in load_uniq(s,'CO')] for s in SURFACES]
d_CH_t= [[r['E']-E_SLAB[s]-E_CH3O_thermo for r in load_uniq(s,'CH3O')] for s in SURFACES]
d_co_t= [[r['E']-E_SLAB[s]-E_CO-E_CH3O_thermo for r in load_coads(s)] for s in SURFACES]
plot_panel(axes[0,0], d_CO,   'a', r'CO$^*$')
plot_panel(axes[0,1], d_CH_t, 'b', r'CH$_3$O$^*$')
plot_panel(axes[0,2], d_co_t, 'c', r'co-ads')
axes[0,0].set_ylabel(r'$E_{\mathrm{bind}}$(CO$^*$) / eV')
axes[0,1].set_ylabel(r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV')
axes[0,2].set_ylabel(r'$E_{\mathrm{bind}}$(CO+CH$_3$O$^*$) / eV')

# Bottom row: radical ref for CH3O (current advisor scheme)
d_CH_r= [[r['E']-E_SLAB[s]-E_CH3O_radical for r in load_uniq(s,'CH3O')] for s in SURFACES]
d_co_r= [[r['E']-E_SLAB[s]-E_CO-E_CH3O_radical for r in load_coads(s)] for s in SURFACES]
plot_panel(axes[1,0], d_CO,   'd', r'CO$^*$')
plot_panel(axes[1,1], d_CH_r, 'e', r'CH$_3$O$^*$')
plot_panel(axes[1,2], d_co_r, 'f', r'co-ads')
axes[1,0].set_ylabel(r'$E_{\mathrm{bind}}$(CO$^*$) / eV')
axes[1,1].set_ylabel(r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV')
axes[1,2].set_ylabel(r'$E_{\mathrm{bind}}$(CO+CH$_3$O$^*$) / eV')

plt.tight_layout()
plt.savefig(OUT/'F14_binding_2refs.png', dpi=160, bbox_inches='tight')
plt.close()
print(f'Saved: {OUT/"F14_binding_2refs.png"}')
