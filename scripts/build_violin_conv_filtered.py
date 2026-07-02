"""Conv A vs Conv B violin, CO2*+intramol filter applied."""
import json, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size']=18; plt.rcParams['axes.labelsize']=22
plt.rcParams['xtick.labelsize']=18; plt.rcParams['ytick.labelsize']=18
plt.rcParams['axes.linewidth']=1.6

ROOT=Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3=ROOT/'calculations/G3_adsorption'; G2=ROOT/'calculations/G2_slab'
OUT=ROOT/'reports/predft_advisor_figures/violin_conv_filtered'
OUT.mkdir(parents=True,exist_ok=True)

SURFACES=['S1','S2','S3','S3b','S4']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100','S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LBL={'S1':'Pd(100)','S2':'1 ML\nPdO(101)/Pd(100)','S3':'O-rich\nPdO(100)','S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}


def valid_idx(sid, ads):
    """Return set of idx that pass intramol+CO2 filter."""
    sdir = SDIRS[sid]
    slab = read(G2/sdir/'CONTCAR')
    traj = list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
    n_ads = 2 if ads=='CO' else 5
    unique = json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
    out = set()
    for r in unique:
        a = traj[r['idx']]
        if len(a) != len(slab)+n_ads:
            ads_atoms = a[-n_ads:]; a = slab.copy(); a += ads_atoms
        syms = a.get_chemical_symbols()
        c_idx = [i for i,s in enumerate(syms) if s=='C']
        if not c_idx: continue
        c = c_idx[-1]
        d_c_o = sorted(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='O')
        if ads == 'CO':
            if not (1.05 <= d_c_o[0] <= 1.30): continue
            if len(d_c_o) >= 2 and d_c_o[1] < 1.5: continue
        else:  # CH3O
            h_idx = [i for i,s in enumerate(syms) if s=='H']
            if len(h_idx) != 3: continue
            d_ch = [a.get_distance(c,h,mic=True) for h in h_idx]
            if not all(0.90<=d<=1.25 for d in d_ch): continue
            if not (1.30 <= d_c_o[0] <= 1.55): continue
            if len(d_c_o) >= 2 and d_c_o[1] < 1.5: continue
        out.add(r['idx'])
    return out


def violin(ax, data, title):
    pos=np.arange(len(SURFACES))
    parts=ax.violinplot(data,positions=pos,widths=0.7,showmeans=False,showmedians=False,showextrema=False)
    for i,pc in enumerate(parts['bodies']):
        pc.set_facecolor(CMAP[SURFACES[i]]); pc.set_alpha(0.4); pc.set_edgecolor('black')
    bp=ax.boxplot(data,positions=pos,widths=0.3,patch_artist=True,showfliers=False,medianprops={'color':'red','lw':2})
    for i,patch in enumerate(bp['boxes']):
        patch.set_facecolor(CMAP[SURFACES[i]]); patch.set_alpha(0.7)
    np.random.seed(42)
    for i,d in enumerate(data):
        if not d: continue
        ax.scatter(np.random.normal(i,0.04,size=len(d)),d,s=3,color='black',alpha=0.3,zorder=2)
        ax.scatter([i],[min(d)],s=120,color='gold',edgecolor='black',linewidth=2,zorder=5,marker='*')
    ax.set_xticks(pos); ax.set_xticklabels([LBL[s] for s in SURFACES],rotation=30,ha='right')
    ax.axhline(0,ls='--',color='gray',alpha=0.6,lw=1)
    ax.set_title(title,fontsize=22,fontweight='bold',pad=10)
    ax.grid(True,alpha=0.3)


for ads, ylab, fname in [('CO',r'$E_{\mathrm{bind}}$(CO$^*$) / eV','CO_convAB_filtered.png'),
                          ('CH3O',r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV','CH3O_convAB_filtered.png')]:
    A_list, B_list = [], []
    for s in SURFACES:
        valid = valid_idx(s, ads)
        recs = json.load(open(G3/'convB_binding'/f'{s}_{ads}_convB.json'))
        recs_f = [r for r in recs if r['idx'] in valid]
        print(f'  {s} {ads}: {len(recs)} → {len(recs_f)} after filter')
        A_list.append([r['E_bind_A_clean'] for r in recs_f])
        B_list.append([r['E_bind_B_frozen'] for r in recs_f])
    all_v = [v for d in A_list+B_list for v in d if d]
    ymin, ymax = min(all_v)-0.2, max(all_v)+0.2
    fig, axes = plt.subplots(1,2,figsize=(20,8), sharey=True)
    violin(axes[0], A_list, 'Convention A\n(clean slab)')
    violin(axes[1], B_list, 'Convention B\n(frozen post-relax slab)')
    axes[0].set_ylabel(ylab)
    for ax in axes: ax.set_ylim(ymin,ymax)
    plt.tight_layout()
    plt.savefig(OUT/fname,dpi=300,bbox_inches='tight'); plt.close()
    print(f'{fname} ✓')

print('\n=== S4 stats (filtered) ===')
for ads in ['CO','CH3O']:
    valid = valid_idx('S4', ads)
    recs = json.load(open(G3/'convB_binding'/f'S4_{ads}_convB.json'))
    recs_f = [r for r in recs if r['idx'] in valid]
    a = [r['E_bind_A_clean'] for r in recs_f]
    b = [r['E_bind_B_frozen'] for r in recs_f]
    if a:
        print(f'{ads}  A min={min(a):.2f} med={np.median(a):.2f}  |  B min={min(b):.2f} med={np.median(b):.2f}')
