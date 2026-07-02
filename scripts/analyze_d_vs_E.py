"""Per-surface: d(Pd-C/O_Me) distribution and E_bind distribution + scatter."""
import json
import numpy as np
from pathlib import Path
from ase.io import read
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size']=16

ROOT=Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3=ROOT/'calculations/G3_adsorption'; G2=ROOT/'calculations/G2_slab'
OUT=ROOT/'reports/predft_advisor_figures/d_vs_E'; OUT.mkdir(parents=True,exist_ok=True)

SURFACES=['S1','S2','S3','S3b','S4']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100','S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs=json.load(open(G3/'mace_d3_references.json'))
E_SLAB=refs['slab']; E_CO=refs['gas']['CO']; E_CH3O=refs['gas']['CH3O_radical']

def analyze(ads):
    n_ads = 2 if ads=='CO' else 5
    by_surf = {}
    for sid in SURFACES:
        sdir = SDIRS[sid]
        slab = read(G2/sdir/'CONTCAR')
        unique = json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
        traj = list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
        rows = []
        for r in unique:
            a = traj[r['idx']]
            if len(a) != len(slab)+n_ads:
                ads_at = a[-n_ads:]; a = slab.copy(); a += ads_at
            syms = a.get_chemical_symbols()
            c_idx = [i for i,s in enumerate(syms) if s=='C']
            if not c_idx: continue
            c = c_idx[-1]
            d_co_sort = sorted([(a.get_distance(c,oi,mic=True),oi) for oi in range(len(a)) if syms[oi]=='O'])
            # Filter
            if ads=='CO':
                if not (1.05<=d_co_sort[0][0]<=1.30): continue
                if len(d_co_sort)>=2 and d_co_sort[1][0]<1.5: continue
                # use d(Pd-C)
                d_pd = min(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='Pd')
            else:  # CH3O
                h_idx = [i for i,s in enumerate(syms) if s=='H']
                if len(h_idx) != 3: continue
                d_ch = [a.get_distance(c,h,mic=True) for h in h_idx]
                if not all(0.90<=d<=1.25 for d in d_ch): continue
                if not (1.30<=d_co_sort[0][0]<=1.55): continue
                if len(d_co_sort)>=2 and d_co_sort[1][0]<1.5: continue
                # use d(Pd-O_methoxy)
                o_me = d_co_sort[0][1]
                d_pd = min(a.get_distance(o_me,i,mic=True) for i in range(len(a)) if syms[i]=='Pd')
            ref = E_CO if ads=='CO' else E_CH3O
            E_bind = r['E'] - E_SLAB[sid] - ref
            rows.append({'idx':r['idx'],'d_pd':d_pd,'E':E_bind})
        by_surf[sid] = rows
    return by_surf

for ads, lab, fname in [('CO','Pd–C', 'CO_d_vs_E.png'),
                         ('CH3O','Pd–O$_\\mathrm{Me}$', 'CH3O_d_vs_E.png')]:
    data = analyze(ads)
    # Scatter
    fig, axes = plt.subplots(1, 2, figsize=(18,7))
    # left: scatter d vs E
    for sid in SURFACES:
        rows = data[sid]
        d = [r['d_pd'] for r in rows]
        e = [r['E'] for r in rows]
        axes[0].scatter(d, e, s=30, alpha=0.5, color=CMAP[sid], label=sid)
    axes[0].set_xlabel(f'd({lab}) / Å'); axes[0].set_ylabel(r'$E_{\mathrm{bind}}$ / eV')
    axes[0].axhline(0, ls='--', color='gray', alpha=0.5)
    axes[0].axvline(2.5, ls='--', color='red', alpha=0.4, label='chem cutoff (2.5 Å)')
    axes[0].grid(True, alpha=0.3); axes[0].legend()
    # right: medians per surface
    pos = np.arange(len(SURFACES))
    d_med = [np.median([r['d_pd'] for r in data[s]]) for s in SURFACES]
    e_med = [np.median([r['E'] for r in data[s]]) for s in SURFACES]
    d_min = [np.min([r['d_pd'] for r in data[s]]) for s in SURFACES]
    ax2 = axes[1]
    ax2b = ax2.twinx()
    ax2.bar(pos-0.2, d_med, width=0.4, color='steelblue', alpha=0.7, label='median d')
    ax2.bar(pos-0.2, d_min, width=0.4, color='steelblue', alpha=0.4, label='min d (closest)')
    ax2b.bar(pos+0.2, e_med, width=0.4, color='orange', alpha=0.7, label='median E_bind')
    ax2.set_xticks(pos); ax2.set_xticklabels(SURFACES)
    ax2.set_ylabel(f'd({lab}) / Å', color='steelblue')
    ax2b.set_ylabel(r'$E_{\mathrm{bind}}$ / eV', color='orange')
    ax2.axhline(2.5, ls='--', color='red', alpha=0.4)
    for i,(d_,e_) in enumerate(zip(d_med, e_med)):
        ax2.text(i-0.2, d_+0.05, f'{d_:.2f}', ha='center', fontsize=11, color='steelblue')
        ax2b.text(i+0.2, e_-0.1, f'{e_:.2f}', ha='center', fontsize=11, color='darkorange')
    fig.suptitle(f'{ads}: distance vs binding energy (filtered, all candidates)', fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUT/fname, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'{fname} ✓')
    # Print numerical summary
    print(f'\n=== {ads} median statistics ===')
    print(f"{'Sur':<5} {'n':<5} {'d_min':<8} {'d_median':<10} {'E_min':<8} {'E_median':<10}")
    for sid in SURFACES:
        rows = data[sid]
        d = [r['d_pd'] for r in rows]
        e = [r['E'] for r in rows]
        print(f'{sid:<5} {len(rows):<5} {min(d):<8.2f} {np.median(d):<10.2f} {min(e):<8.2f} {np.median(e):<10.2f}')
