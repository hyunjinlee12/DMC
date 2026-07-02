"""Individual violins (no titles) — before/after × CO/CH3O = 4 files."""
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
OUT=ROOT/'reports/predft_advisor_figures/final_comparison'

SURFACES=['S1','S2','S3','S3b','S4']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100','S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LBL={'S1':'Pd(100)','S2':'1 ML\nPdO(101)/Pd(100)','S3':'O-rich\nPdO(100)','S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}


def filt(sid, ads, strict):
    sdir=SDIRS[sid]; slab=read(G2/sdir/'CONTCAR')
    traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj',index=':'))
    n=2 if ads=='CO' else 5
    unique=json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
    convB_recs={x['idx']:x for x in json.load(open(G3/'convB_binding'/f'{sid}_{ads}_convB.json'))}
    out=set(); co2_cut = 2.0 if strict else 1.5
    for r in unique:
        if not r.get('converged',True): continue
        a=traj[r['idx']]
        if len(a)!=len(slab)+n: ads_a=a[-n:]; a=slab.copy(); a+=ads_a
        syms=a.get_chemical_symbols()
        c_idx=[i for i,s in enumerate(syms) if s=='C']
        if not c_idx: continue
        c=c_idx[-1]
        d_co=sorted(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='O')
        if ads=='CO':
            if not (1.05<=d_co[0]<=1.30): continue
            if len(d_co)>=2 and d_co[1]<co2_cut: continue
            if strict:
                d_pdc=min(a.get_distance(c,i,mic=True) for i in range(len(a)) if syms[i]=='Pd')
                E_B=convB_recs[r['idx']]['E_bind_B_frozen']
                if d_pdc>3.0 and E_B<-1.0: continue
        else:
            h_idx=[i for i,s in enumerate(syms) if s=='H']
            if len(h_idx)!=3: continue
            d_ch=[a.get_distance(c,h,mic=True) for h in h_idx]
            if not all(0.90<=d<=1.25 for d in d_ch): continue
            if not (1.30<=d_co[0]<=1.55): continue
            if len(d_co)>=2 and d_co[1]<co2_cut: continue
            if strict:
                o_me=[i for i in range(len(a)) if syms[i]=='O' and a.get_distance(i,c,mic=True)<1.6][0]
                d_pdo=min(a.get_distance(o_me,i,mic=True) for i in range(len(a)) if syms[i]=='Pd')
                E_B=convB_recs[r['idx']]['E_bind_B_frozen']
                if d_pdo>3.0 and E_B<-1.0: continue
        out.add(r['idx'])
    return out


def violin_one(data, ylab, fname, ylim=None):
    fig, ax = plt.subplots(figsize=(11,8))
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
    ax.set_xticks(pos); ax.set_xticklabels(['']*len(SURFACES))
    ax.set_yticklabels([])
    ax.axhline(0,ls='--',color='gray',alpha=0.6,lw=1)
    if ylim: ax.set_ylim(ylim)
    ax.grid(True,alpha=0.3); plt.tight_layout()
    plt.savefig(OUT/fname,dpi=300,bbox_inches='tight'); plt.close()
    print(f'{fname} ✓')


# Build CO before/after with shared ylim
co_before=[]; co_after=[]
for s in SURFACES:
    recs=json.load(open(G3/'convB_binding'/f'{s}_CO_convB.json'))
    vb=filt(s,'CO',False); va=filt(s,'CO',True)
    co_before.append([r['E_bind_B_frozen'] for r in recs if r['idx'] in vb])
    co_after.append([r['E_bind_B_frozen'] for r in recs if r['idx'] in va])
all_v=[v for d in co_before+co_after for v in d]
ylim_co = (min(all_v)-0.2, max(all_v)+0.2)
violin_one(co_before, r'$E_{\mathrm{bind}}^{B}$(CO$^*$) / eV', 'CO_before_only.png', ylim_co)
violin_one(co_after,  r'$E_{\mathrm{bind}}^{B}$(CO$^*$) / eV', 'CO_after_only.png',  ylim_co)

# CH3O
ch_before=[]; ch_after=[]
for s in SURFACES:
    recs=json.load(open(G3/'convB_binding'/f'{s}_CH3O_convB.json'))
    vb=filt(s,'CH3O',False); va=filt(s,'CH3O',True)
    ch_before.append([r['E_bind_B_frozen'] for r in recs if r['idx'] in vb])
    ch_after.append([r['E_bind_B_frozen'] for r in recs if r['idx'] in va])
all_v=[v for d in ch_before+ch_after for v in d]
ylim_ch = (min(all_v)-0.2, max(all_v)+0.2)
violin_one(ch_before, r'$E_{\mathrm{bind}}^{B}$(CH$_3$O$^*$) / eV', 'CH3O_before_only.png', ylim_ch)
violin_one(ch_after,  r'$E_{\mathrm{bind}}^{B}$(CH$_3$O$^*$) / eV', 'CH3O_after_only.png',  ylim_ch)
