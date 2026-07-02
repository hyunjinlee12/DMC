"""Build F13_binding_energy.png — proper MACE+D3 binding energy distribution.

Y-axis: E_bind (eV)
  E_bind(CO*)    = E(slab+CO)    - E(slab)        - E(CO_gas)
  E_bind(CH3O*)  = E(slab+CH3O)  - E(slab)        - [E(CH3OH) - 1/2 E(H2)]
  E_bind(co-ads) = E(slab+co)    - E(slab)        - E(CO) - [E(CH3OH) - 1/2 E(H2)]

All E values from MACE-MH+D3+cueq (same calculator). Negative E_bind = bound.
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
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
OUT = ROOT / 'reports/predft_advisor_figures'
OUT.mkdir(parents=True, exist_ok=True)

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
COLORS = {'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs = json.load(open(G3 / 'mace_d3_references.json'))
E_SLAB = refs['slab']
E_CO_gas = refs['gas']['CO']
E_CH3O_ref = refs['gas']['CH3O_radical']
E_CH3O_thermo = refs['gas']['CH3O_ref']   # = CH3OH - 1/2 H2 (for co-ads per advisor)


def load_unique(sid, ads):
    f = G3 / SDIRS[sid] / 'MLIP_phase1' / f'unique_{ads}.json'
    return json.load(open(f)) if f.exists() else []


def load_unique_coads(sid):
    f = G3 / SDIRS[sid] / 'MLIP_phase2_filtered' / 'unique_SetA.json'
    return json.load(open(f)) if f.exists() else []


def binding_co(recs, sid):
    return [r['E'] - E_SLAB[sid] - E_CO_gas for r in recs]


def binding_ch3o(recs, sid):
    return [r['E'] - E_SLAB[sid] - E_CH3O_ref for r in recs]


def binding_coads(recs, sid):
    return [r['E'] - E_SLAB[sid] - E_CO_gas - E_CH3O_thermo for r in recs]


fig, axes = plt.subplots(1, 3, figsize=(18, 6.5))


def plot_panel(ax, data_list, label_letter, ads_label):
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
        x_jit = np.random.normal(i, 0.04, size=len(d))
        ax.scatter(x_jit, d, s=4, color='black', alpha=0.3, zorder=2)
        # mark top-1 (most negative = most bound)
        e_min = min(d)
        ax.scatter([i], [e_min], s=140, color='gold', edgecolor='black', linewidth=2, zorder=5,
                   marker='*')
    ax.set_xticks(positions); ax.set_xticklabels(SURFACES)
    ax.axhline(0, ls='--', color='gray', alpha=0.6, lw=1)
    ax.text(0.04, 0.96, f'({label_letter}) {ads_label}', transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
    ax.grid(True, alpha=0.3)


# === (a) CO* ===
data_CO = [binding_co(load_unique(sid, 'CO'), sid) for sid in SURFACES]
plot_panel(axes[0], data_CO, 'a', r'CO$^*$')
axes[0].set_ylabel(r'$E_{\mathrm{bind}}$(CO$^*$) / eV')

# === (b) CH3O* ===
data_CH = [binding_ch3o(load_unique(sid, 'CH3O'), sid) for sid in SURFACES]
plot_panel(axes[1], data_CH, 'b', r'CH$_3$O$^*$')
axes[1].set_ylabel(r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV')

# === (c) co-ads SetA ===
data_coads = [binding_coads(load_unique_coads(sid), sid) for sid in SURFACES]
plot_panel(axes[2], data_coads, 'c', r'co-ads (SetA)')
axes[2].set_ylabel(r'$E_{\mathrm{bind}}$(CO$^*$+CH$_3$O$^*$) / eV')

plt.tight_layout()
plt.savefig(OUT/'F13_binding_energy.png', dpi=160, bbox_inches='tight')
plt.close()
print(f'Saved: {OUT/"F13_binding_energy.png"}')

# Summary table
print('\n=== Binding energy (MACE+D3) summary per surface × adsorbate ===')
print(f"{'Sur':<5} {'ads':<8} {'n':<5} {'min/top-1':<11} {'median':<10} {'max':<10}")
for sid, dCO, dCH, dco in zip(SURFACES, data_CO, data_CH, data_coads):
    for name, d in [('CO', dCO), ('CH3O', dCH), ('coads', dco)]:
        if not d: continue
        d = np.array(d)
        print(f"  {sid:<5} {name:<8} {len(d):<5} {min(d):<11.3f} {np.median(d):<10.3f} {max(d):<10.3f}")

print('\nReferences used (MACE+D3):')
print(f"  E_CO_gas    = {E_CO_gas:.4f} eV")
print(f"  E_CH3O_ref  = {E_CH3O_ref:.4f} eV  [= E(CH3OH) - 1/2 E(H2)]")
for sid in SURFACES:
    print(f"  E_slab({sid}) = {E_SLAB[sid]:.4f} eV")
