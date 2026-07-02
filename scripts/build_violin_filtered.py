"""Filtered violin plots — exclude CO2*-formation, intramol-broken, and decomposed candidates."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size'] = 22
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.linewidth'] = 1.8

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT/'calculations/G3_adsorption'
G2 = ROOT/'calculations/G2_slab'
OUT = ROOT/'reports/predft_advisor_figures/violin_filtered'
OUT.mkdir(parents=True, exist_ok=True)

SURFACES = ['S1','S2','S3','S3b','S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
SURF_LABEL = {'S1':'Pd(100)','S2':'1 ML PdO(101)/Pd(100)',
              'S3':'O-rich PdO(100)','S3b':'Pd-rich PdO(100)',
              'S4':r'PdO$_2$(110)'}
COLORS = {'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs = json.load(open(G3/'mace_d3_references.json'))
E_SLAB = refs['slab']
E_CO = refs['gas']['CO']
E_CH3O = refs['gas']['CH3O_radical']

CO2_BOND_MAX = 1.50      # if d(C - other lattice O) < this → CO2 formation
INTRAMOL = {'CO':(1.05, 1.30), 'CH3O_OC':(1.30, 1.55), 'CH3O_CH':(0.90, 1.25)}
SUB_O_DISP_MAX = 0.30    # max allowed substrate O displacement from clean slab (Å)


def slab_o_intact(atoms, slab, n_ads):
    """Check substrate O atoms moved < SUB_O_DISP_MAX vs clean slab."""
    sub = atoms[:len(slab)]
    syms = sub.get_chemical_symbols()
    disp = np.linalg.norm(sub.positions - slab.positions, axis=1)
    o_disp = [d for s,d in zip(syms, disp) if s=='O']
    if not o_disp: return True
    return max(o_disp) < SUB_O_DISP_MAX


def filter_co(sid, recs):
    """For CO: pass intramol + no CO2* formation."""
    sdir = SDIRS[sid]
    traj = list(read(G3/sdir/'MLIP_phase1/relaxed_CO.traj', index=':'))
    slab = read(G2/sdir/'CONTCAR')
    out = []
    for r in recs:
        a = traj[r['idx']]
        if len(a) != len(slab)+2:
            ads = a[-2:]; a = slab.copy(); a += ads
        syms = a.get_chemical_symbols()
        c_idx = [i for i,s in enumerate(syms) if s=='C']
        if not c_idx: continue
        c = c_idx[-1]
        o_all = [i for i,s in enumerate(syms) if s=='O']
        d_c_o = sorted([(a.get_distance(c,oi,mic=True), oi) for oi in o_all])
        # d_c_o[0] = own CO oxygen, d_c_o[1] = nearest lattice O
        if not (INTRAMOL['CO'][0] <= d_c_o[0][0] <= INTRAMOL['CO'][1]):
            continue
        if len(d_c_o) > 1 and d_c_o[1][0] < CO2_BOND_MAX:
            continue   # CO2* formation
        if not slab_o_intact(a, slab, 2):
            continue   # slab O rearrangement artifact
        out.append(r)
    return out


def filter_ch3o(sid, recs):
    """For CH3O: pass intramol + no extra C-O_lat bond (avoid CH3 detached/CO2-like)."""
    sdir = SDIRS[sid]
    traj = list(read(G3/sdir/'MLIP_phase1/relaxed_CH3O.traj', index=':'))
    slab = read(G2/sdir/'CONTCAR')
    out = []
    for r in recs:
        a = traj[r['idx']]
        if len(a) != len(slab)+5:
            ads = a[-5:]; a = slab.copy(); a += ads
        syms = a.get_chemical_symbols()
        c_idx = [i for i,s in enumerate(syms) if s=='C']
        h_idx = [i for i,s in enumerate(syms) if s=='H']
        o_all = [i for i,s in enumerate(syms) if s=='O']
        if not c_idx or len(h_idx) != 3: continue
        c = c_idx[-1]
        d_c_h = [a.get_distance(c,h,mic=True) for h in h_idx]
        d_c_o = sorted([(a.get_distance(c,oi,mic=True), oi) for oi in o_all])
        if not all(INTRAMOL['CH3O_CH'][0] <= d <= INTRAMOL['CH3O_CH'][1] for d in d_c_h):
            continue
        if not (INTRAMOL['CH3O_OC'][0] <= d_c_o[0][0] <= INTRAMOL['CH3O_OC'][1]):
            continue
        # exclude additional C-O_lattice bond (only 1 O attached to C = methoxy)
        if len(d_c_o) > 1 and d_c_o[1][0] < CO2_BOND_MAX:
            continue
        if not slab_o_intact(a, slab, 5):
            continue   # slab O rearrangement artifact
        out.append(r)
    return out


def load_uniq(sid, ads):
    f = G3/SDIRS[sid]/'MLIP_phase1'/f'unique_{ads}.json'
    return json.load(open(f)) if f.exists() else []


def load_coads(sid):
    f = G3/SDIRS[sid]/'MLIP_phase2_filtered/unique_SetA.json'
    return json.load(open(f)) if f.exists() else []


def violin_panel(data_list, letter, ads_label, ylabel, fname):
    fig, ax = plt.subplots(figsize=(11, 8))
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
    ax.set_xticks(positions)
    ax.set_xticklabels([SURF_LABEL[s] for s in SURFACES], rotation=30, ha='right')
    ax.set_ylabel(ylabel)
    ax.axhline(0, ls='--', color='gray', alpha=0.6, lw=1)
    ax.text(0.04, 0.96, f'({letter}) {ads_label}', transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT/fname, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{fname} ✓')


# CO: apply filter
d_CO = []
for s in SURFACES:
    recs = load_uniq(s,'CO')
    n_pre = len(recs)
    recs_f = filter_co(s, recs)
    print(f'  {s} CO: {n_pre} → {len(recs_f)} after filter')
    d_CO.append([r['E']-E_SLAB[s]-E_CO for r in recs_f])

# CH3O: apply filter
d_CH = []
for s in SURFACES:
    recs = load_uniq(s,'CH3O')
    n_pre = len(recs)
    recs_f = filter_ch3o(s, recs)
    print(f'  {s} CH3O: {n_pre} → {len(recs_f)} after filter')
    d_CH.append([r['E']-E_SLAB[s]-E_CH3O for r in recs_f])

# co-ads: leave as-is (SetA already had intramol + band filter)
d_co = [[r['E']-E_SLAB[s]-E_CO-E_CH3O for r in load_coads(s)] for s in SURFACES]

violin_panel(d_CO, 'a', r'CO$^*$  (filtered)', r'$E_{\mathrm{bind}}$(CO$^*$) / eV', 'F13a_CO_filtered.png')
violin_panel(d_CH, 'b', r'CH$_3$O$^*$  (filtered)', r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV', 'F13b_CH3O_filtered.png')
violin_panel(d_co, 'c', r'co-ads  (SetA)', r'$E_{\mathrm{bind}}$(CO$^*$+CH$_3$O$^*$) / eV', 'F13c_coads_filtered.png')
