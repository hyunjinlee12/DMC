"""Individual high-DPI: Conv B violin + d-vs-E scatter."""
import json, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size']=20; plt.rcParams['axes.labelsize']=24
plt.rcParams['xtick.labelsize']=20; plt.rcParams['ytick.labelsize']=20
plt.rcParams['axes.linewidth']=1.8

ROOT=Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3=ROOT/'calculations/G3_adsorption'; G2=ROOT/'calculations/G2_slab'
OUT=ROOT/'reports/predft_advisor_figures/final_individual'
OUT.mkdir(parents=True,exist_ok=True)

SURFACES=['S1','S2','S3','S3b','S4']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100','S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LBL={'S1':'Pd(100)','S2':'1 ML\nPdO(101)/Pd(100)','S3':'O-rich\nPdO(100)','S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs=json.load(open(G3/'mace_d3_references.json'))
E_SLAB=refs['slab']; E_CO=refs['gas']['CO']; E_CH3O=refs['gas']['CH3O_radical']


def valid_idx(sid, ads):
    sdir=SDIRS[sid]; slab=read(G2/sdir/'CONTCAR')
    traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj',index=':'))
    n=2 if ads=='CO' else 5
    unique=json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
    out=set()
    for r in unique:
        a=traj[r['idx']]
        if len(a)!=len(slab)+n: ads_a=a[-n:]; a=slab.copy(); a+=ads_a
        syms=a.get_chemical_symbols()
        c_idx=[i for i,s in enumerate(syms) if s=='C']
        if not c_idx: continue
        c=c_idx[-1]
        d_co=sorted(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='O')
        if ads=='CO':
            if not (1.05<=d_co[0]<=1.30): continue
            if len(d_co)>=2 and d_co[1]<1.5: continue
        else:
            h_idx=[i for i,s in enumerate(syms) if s=='H']
            if len(h_idx)!=3: continue
            d_ch=[a.get_distance(c,h,mic=True) for h in h_idx]
            if not all(0.90<=d<=1.25 for d in d_ch): continue
            if not (1.30<=d_co[0]<=1.55): continue
            if len(d_co)>=2 and d_co[1]<1.5: continue
        out.add(r['idx'])
    return out


def violin_panel(data_list, letter, ads_label, ylabel, fname):
    fig, ax = plt.subplots(figsize=(11,8))
    pos=np.arange(len(SURFACES))
    parts=ax.violinplot(data_list,positions=pos,widths=0.7,showmeans=False,showmedians=False,showextrema=False)
    for i,pc in enumerate(parts['bodies']):
        pc.set_facecolor(CMAP[SURFACES[i]]); pc.set_alpha(0.4); pc.set_edgecolor('black')
    bp=ax.boxplot(data_list,positions=pos,widths=0.3,patch_artist=True,showfliers=False,medianprops={'color':'red','lw':2})
    for i,patch in enumerate(bp['boxes']):
        patch.set_facecolor(CMAP[SURFACES[i]]); patch.set_alpha(0.7)
    np.random.seed(42)
    for i,d in enumerate(data_list):
        if not d: continue
        ax.scatter(np.random.normal(i,0.04,size=len(d)),d,s=4,color='black',alpha=0.3,zorder=2)
        ax.scatter([i],[min(d)],s=140,color='gold',edgecolor='black',linewidth=2,zorder=5,marker='*')
    ax.set_xticks(pos); ax.set_xticklabels([LBL[s] for s in SURFACES],rotation=30,ha='right')
    ax.set_ylabel(ylabel); ax.axhline(0,ls='--',color='gray',alpha=0.6,lw=1)
    ax.text(0.04,0.96,f'({letter}) {ads_label}  (Conv B, filtered)',transform=ax.transAxes,
            fontsize=20,fontweight='bold',va='top',
            bbox=dict(boxstyle='round',facecolor='white',alpha=0.85))
    ax.grid(True,alpha=0.3); plt.tight_layout()
    plt.savefig(OUT/fname,dpi=300,bbox_inches='tight'); plt.close()
    print(f'{fname} ✓')


# Conv B individual violins
for ads, letter, lab, ylab, fname in [
    ('CO','a',r'CO$^*$',r'$E_{\mathrm{bind}}^{B}$(CO$^*$) / eV','convB_CO.png'),
    ('CH3O','b',r'CH$_3$O$^*$',r'$E_{\mathrm{bind}}^{B}$(CH$_3$O$^*$) / eV','convB_CH3O.png'),
]:
    data=[]
    for s in SURFACES:
        valid=valid_idx(s,ads)
        recs=json.load(open(G3/'convB_binding'/f'{s}_{ads}_convB.json'))
        recs_f=[r for r in recs if r['idx'] in valid]
        data.append([r['E_bind_B_frozen'] for r in recs_f])
    violin_panel(data, letter, lab, ylab, fname)


# d vs E scatter individual
def d_vs_E_one(ads, fname, xlab):
    n_ads=2 if ads=='CO' else 5
    fig, ax = plt.subplots(figsize=(11,8))
    for sid in SURFACES:
        sdir=SDIRS[sid]
        slab=read(G2/sdir/'CONTCAR')
        unique=json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
        traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj',index=':'))
        d_list, e_list = [], []
        for r in unique:
            a=traj[r['idx']]
            if len(a)!=len(slab)+n_ads: ads_a=a[-n_ads:]; a=slab.copy(); a+=ads_a
            syms=a.get_chemical_symbols()
            c_idx=[i for i,s in enumerate(syms) if s=='C']
            if not c_idx: continue
            c=c_idx[-1]
            d_co_sort=sorted([(a.get_distance(c,oi,mic=True),oi) for oi in range(len(a)) if syms[oi]=='O'])
            if ads=='CO':
                if not (1.05<=d_co_sort[0][0]<=1.30): continue
                if len(d_co_sort)>=2 and d_co_sort[1][0]<1.5: continue
                d=min(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='Pd')
            else:
                h=[i for i,s in enumerate(syms) if s=='H']
                if len(h)!=3: continue
                if not all(0.90<=a.get_distance(c,hi,mic=True)<=1.25 for hi in h): continue
                if not (1.30<=d_co_sort[0][0]<=1.55): continue
                if len(d_co_sort)>=2 and d_co_sort[1][0]<1.5: continue
                o_me=d_co_sort[0][1]
                d=min(a.get_distance(o_me,i,mic=True) for i in range(len(a)) if syms[i]=='Pd')
            ref=E_CO if ads=='CO' else E_CH3O
            e=r['E']-E_SLAB[sid]-ref
            d_list.append(d); e_list.append(e)
        ax.scatter(d_list, e_list, s=60, alpha=0.6, color=CMAP[sid], label=LBL[sid].replace('\n',' '),
                   edgecolor='black', linewidth=0.5)
    ax.set_xlabel(f'd({xlab}) / Å')
    ax.set_ylabel(r'$E_{\mathrm{bind}}$ / eV')
    ax.axhline(0, ls='--', color='gray', alpha=0.5, lw=1)
    ax.axvline(2.5, ls='--', color='red', alpha=0.4, lw=1.5, label='chem cutoff (2.5 Å)')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=14, loc='best', frameon=True)
    plt.tight_layout()
    plt.savefig(OUT/fname, dpi=300, bbox_inches='tight'); plt.close()
    print(f'{fname} ✓')


d_vs_E_one('CO', 'd_vs_E_CO.png', 'Pd–C')
d_vs_E_one('CH3O', 'd_vs_E_CH3O.png', r'Pd–O$_\mathrm{Me}$')
