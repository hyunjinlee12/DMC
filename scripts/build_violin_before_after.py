"""Before vs after CO2 filter — same ylim for fair comparison."""
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
OUT=ROOT/'reports/predft_advisor_figures/violin_before_after'; OUT.mkdir(parents=True,exist_ok=True)

SURFACES=['S1','S2','S3','S3b','S4']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100','S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LBL={'S1':'Pd(100)','S2':'1 ML\nPdO(101)/Pd(100)','S3':'O-rich\nPdO(100)','S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs=json.load(open(G3/'mace_d3_references.json'))
E_SLAB=refs['slab']; E_CO=refs['gas']['CO']; E_CH3O=refs['gas']['CH3O_radical']

def filt_co(sid, recs):
    slab=read(G2/SDIRS[sid]/'CONTCAR')
    traj=list(read(G3/SDIRS[sid]/'MLIP_phase1/relaxed_CO.traj',index=':'))
    out=[]
    for r in recs:
        a=traj[r['idx']]
        if len(a)!=len(slab)+2: ads=a[-2:]; a=slab.copy(); a+=ads
        syms=a.get_chemical_symbols()
        c=[i for i,s in enumerate(syms) if s=='C'][-1]
        d=sorted(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='O')
        if not (1.05<=d[0]<=1.30): continue
        if len(d)>=2 and d[1]<1.5: continue
        out.append(r)
    return out

def filt_ch3o(sid, recs):
    slab=read(G2/SDIRS[sid]/'CONTCAR')
    traj=list(read(G3/SDIRS[sid]/'MLIP_phase1/relaxed_CH3O.traj',index=':'))
    out=[]
    for r in recs:
        a=traj[r['idx']]
        if len(a)!=len(slab)+5: ads=a[-5:]; a=slab.copy(); a+=ads
        syms=a.get_chemical_symbols()
        c=[i for i,s in enumerate(syms) if s=='C']; h=[i for i,s in enumerate(syms) if s=='H']
        if not c or len(h)!=3: continue
        c=c[-1]
        d_ch=[a.get_distance(c,hi,mic=True) for hi in h]
        d_co=sorted(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='O')
        if not all(0.90<=x<=1.25 for x in d_ch): continue
        if not (1.30<=d_co[0]<=1.55): continue
        if len(d_co)>=2 and d_co[1]<1.5: continue
        out.append(r)
    return out

def load_uniq(sid,ads):
    f=G3/SDIRS[sid]/'MLIP_phase1'/f'unique_{ads}.json'
    return json.load(open(f)) if f.exists() else []

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

def cmp(ads, ylab, fname):
    if ads=='CO': filt=filt_co
    else: filt=filt_ch3o
    data_b=[]; data_a=[]
    for s in SURFACES:
        raw=load_uniq(s,ads); f=filt(s,raw)
        ref = E_CO if ads=='CO' else E_CH3O
        data_b.append([r['E']-E_SLAB[s]-ref for r in raw])
        data_a.append([r['E']-E_SLAB[s]-ref for r in f])
    # Same ylim
    all_v = [v for d in data_b+data_a for v in d]
    ymin,ymax = min(all_v)-0.2, max(all_v)+0.2

    fig, axes = plt.subplots(1,2,figsize=(20,8), sharey=True)
    violin(axes[0], data_b, 'Before filter')
    violin(axes[1], data_a, 'After CO$_2$ filter')
    axes[0].set_ylabel(ylab)
    for ax in axes: ax.set_ylim(ymin,ymax)
    plt.tight_layout()
    plt.savefig(OUT/fname,dpi=300,bbox_inches='tight'); plt.close()
    print(f'{fname} ✓')

cmp('CO',  r'$E_{\mathrm{bind}}$(CO$^*$) / eV',  'CO_before_after.png')
cmp('CH3O',r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV','CH3O_before_after.png')
