"""Violin: only CO2* + intramol filter (no slab O rearrangement filter)."""
import json, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size']=22; plt.rcParams['axes.labelsize']=26
plt.rcParams['xtick.labelsize']=22; plt.rcParams['ytick.labelsize']=22
plt.rcParams['axes.linewidth']=1.8

ROOT=Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3=ROOT/'calculations/G3_adsorption'; G2=ROOT/'calculations/G2_slab'
OUT=ROOT/'reports/predft_advisor_figures/violin_co2only'; OUT.mkdir(parents=True,exist_ok=True)

SURFACES=['S1','S2','S3','S3b','S4']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100','S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LBL={'S1':'Pd(100)','S2':'1 ML PdO(101)/Pd(100)','S3':'O-rich PdO(100)','S3b':'Pd-rich PdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs=json.load(open(G3/'mace_d3_references.json'))
E_SLAB=refs['slab']; E_CO=refs['gas']['CO']; E_CH3O=refs['gas']['CH3O_radical']

def filt_co(sid, recs):
    slab=read(G2/SDIRS[sid]/'CONTCAR')
    traj=list(read(G3/SDIRS[sid]/'MLIP_phase1/relaxed_CO.traj',index=':'))
    out=[]
    for r in recs:
        a=traj[r['idx']]
        if len(a)!=len(slab)+2:
            ads=a[-2:]; a=slab.copy(); a+=ads
        syms=a.get_chemical_symbols()
        c=[i for i,s in enumerate(syms) if s=='C'][-1]
        d=sorted(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='O')
        if not (1.05<=d[0]<=1.30): continue
        if len(d)>=2 and d[1]<1.5: continue   # CO2*
        out.append(r)
    return out

def filt_ch3o(sid, recs):
    slab=read(G2/SDIRS[sid]/'CONTCAR')
    traj=list(read(G3/SDIRS[sid]/'MLIP_phase1/relaxed_CH3O.traj',index=':'))
    out=[]
    for r in recs:
        a=traj[r['idx']]
        if len(a)!=len(slab)+5:
            ads=a[-5:]; a=slab.copy(); a+=ads
        syms=a.get_chemical_symbols()
        c=[i for i,s in enumerate(syms) if s=='C']
        h=[i for i,s in enumerate(syms) if s=='H']
        if not c or len(h)!=3: continue
        c=c[-1]
        d_ch=[a.get_distance(c,hi,mic=True) for hi in h]
        d_co=sorted(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='O')
        if not all(0.90<=x<=1.25 for x in d_ch): continue
        if not (1.30<=d_co[0]<=1.55): continue
        if len(d_co)>=2 and d_co[1]<1.5: continue   # CO2-like
        out.append(r)
    return out

def load_uniq(sid,ads):
    f=G3/SDIRS[sid]/'MLIP_phase1'/f'unique_{ads}.json'
    return json.load(open(f)) if f.exists() else []
def load_coads(sid):
    f=G3/SDIRS[sid]/'MLIP_phase2_filtered/unique_SetA.json'
    return json.load(open(f)) if f.exists() else []

def panel(data, letter, lab, ylab, fname):
    fig,ax=plt.subplots(figsize=(11,8))
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
        ax.scatter(np.random.normal(i,0.04,size=len(d)),d,s=4,color='black',alpha=0.3,zorder=2)
        ax.scatter([i],[min(d)],s=140,color='gold',edgecolor='black',linewidth=2,zorder=5,marker='*')
    ax.set_xticks(pos); ax.set_xticklabels([LBL[s] for s in SURFACES],rotation=30,ha='right')
    ax.set_ylabel(ylab); ax.axhline(0,ls='--',color='gray',alpha=0.6,lw=1)
    ax.text(0.04,0.96,f'({letter}) {lab}',transform=ax.transAxes,fontsize=22,fontweight='bold',
            va='top',bbox=dict(boxstyle='round',facecolor='white',alpha=0.85))
    ax.grid(True,alpha=0.3); plt.tight_layout()
    plt.savefig(OUT/fname,dpi=300,bbox_inches='tight'); plt.close()
    print(f'{fname} ✓')

d_CO=[]
for s in SURFACES:
    raw=load_uniq(s,'CO'); f=filt_co(s,raw)
    print(f'  {s} CO: {len(raw)} → {len(f)}')
    d_CO.append([r['E']-E_SLAB[s]-E_CO for r in f])
d_CH=[]
for s in SURFACES:
    raw=load_uniq(s,'CH3O'); f=filt_ch3o(s,raw)
    print(f'  {s} CH3O: {len(raw)} → {len(f)}')
    d_CH.append([r['E']-E_SLAB[s]-E_CH3O for r in f])
d_co=[[r['E']-E_SLAB[s]-E_CO-E_CH3O for r in load_coads(s)] for s in SURFACES]

panel(d_CO,'a',r'CO$^*$',r'$E_{\mathrm{bind}}$(CO$^*$) / eV','F13a_CO_co2only.png')
panel(d_CH,'b',r'CH$_3$O$^*$',r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV','F13b_CH3O_co2only.png')
panel(d_co,'c',r'co-ads',r'$E_{\mathrm{bind}}$(CO$^*$+CH$_3$O$^*$) / eV','F13c_coads_co2only.png')
