"""Before vs After excluding 4 S4 outliers (idx=33, 80, 76, 61).
   Both panels: converged + intramol valid + CO2 (1.5) + Conv B.
   After adds: CO2 (2.0) + geometric-energetic mismatch.
Plus: render 4 excluded structures."""
import json, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read
import sys
sys.path.insert(0, str(Path(__file__).parent))
from render_structure import render

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size']=18; plt.rcParams['axes.labelsize']=22
plt.rcParams['xtick.labelsize']=18; plt.rcParams['ytick.labelsize']=18
plt.rcParams['axes.linewidth']=1.6

ROOT=Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3=ROOT/'calculations/G3_adsorption'; G2=ROOT/'calculations/G2_slab'
OUT=ROOT/'reports/predft_advisor_figures/final_comparison'
OUT.mkdir(parents=True,exist_ok=True)

SURFACES=['S1','S2','S3','S3b','S4']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100','S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LBL={'S1':'Pd(100)','S2':'1 ML\nPdO(101)/Pd(100)','S3':'O-rich\nPdO(100)','S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}


def filt_co_strict(sid, ads, strict):
    """strict=False: CO2 cutoff 1.5, no geom-energetic filter
       strict=True:  CO2 cutoff 2.0 + geom-energetic mismatch filter"""
    sdir=SDIRS[sid]; slab=read(G2/sdir/'CONTCAR')
    traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj',index=':'))
    n=2 if ads=='CO' else 5
    unique=json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
    convB_recs={x['idx']:x for x in json.load(open(G3/'convB_binding'/f'{sid}_{ads}_convB.json'))}
    out=set()
    co2_cut = 2.0 if strict else 1.5
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
    ax.set_title(title,fontsize=20,fontweight='bold',pad=10)
    ax.grid(True,alpha=0.3)


# Build data + medians
print('\n=== MEDIAN COMPARISON ===')
print(f"{'Sur':<5} {'ads':<6} {'before n':<9} {'before med':<11} {'after n':<9} {'after med':<11}")
for ads, ylab, fname in [
    ('CO',  r'$E_{\mathrm{bind}}^{B}$(CO$^*$) / eV',  'CO_before_after.png'),
    ('CH3O',r'$E_{\mathrm{bind}}^{B}$(CH$_3$O$^*$) / eV','CH3O_before_after.png'),
]:
    before=[]; after=[]
    for s in SURFACES:
        recs=json.load(open(G3/'convB_binding'/f'{s}_{ads}_convB.json'))
        valid_before=filt_co_strict(s, ads, strict=False)
        valid_after =filt_co_strict(s, ads, strict=True)
        b=[r['E_bind_B_frozen'] for r in recs if r['idx'] in valid_before]
        a=[r['E_bind_B_frozen'] for r in recs if r['idx'] in valid_after]
        before.append(b); after.append(a)
        if b and a:
            print(f"{s:<5} {ads:<6} {len(b):<9} {np.median(b):<+11.2f} {len(a):<9} {np.median(a):<+11.2f}")
    all_v=[v for d in before+after for v in d]
    ymin,ymax=min(all_v)-0.2, max(all_v)+0.2
    fig, axes=plt.subplots(1,2,figsize=(20,8), sharey=True)
    violin(axes[0], before, 'Before excluding 4 outliers')
    violin(axes[1], after,  'After excluding 4 outliers')
    axes[0].set_ylabel(ylab)
    for ax in axes: ax.set_ylim(ymin,ymax)
    plt.tight_layout()
    plt.savefig(OUT/fname,dpi=300,bbox_inches='tight'); plt.close()
    print(f'{fname} ✓')


# Render 4 excluded structures (S4 CO idx=33, 80, 76, 61)
print('\n=== Rendering 4 excluded structures ===')
slab=read(G2/'S4_PdO2_110/CONTCAR')
traj=list(read(G3/'S4_PdO2_110/MLIP_phase1/relaxed_CO.traj',index=':'))
fig, axes = plt.subplots(1,4,figsize=(20,7))
TITLES = {
    33: 'idx=33  (Pd-C 2.05+CO$_2$ boundary)',
    80: 'idx=80  (floating 3.90 Å)',
    76: 'idx=76  (floating 3.60 Å)',
    61: 'idx=61  (floating 3.52 Å)',
}
import matplotlib.image as mpimg
for ax, idx in zip(axes, [33, 80, 76, 61]):
    a=traj[idx]
    if len(a)!=len(slab)+2: ads=a[-2:]; a=slab.copy(); a+=ads
    a=a*(2,2,1)
    p = OUT / f'tmp_idx{idx}.png'
    render(a, p, rotation='-90x,-90y,0z', width=1500, show_cell=False)
    ax.imshow(mpimg.imread(p)); ax.set_axis_off()
    ax.set_title(TITLES[idx], fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(OUT/'4_excluded_structures.png', dpi=300, bbox_inches='tight')
plt.close()
print('4 excluded ✓')
