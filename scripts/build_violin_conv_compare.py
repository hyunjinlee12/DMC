"""Convention A (clean slab) vs Convention B (frozen post-relax slab) — violin compare."""
import json, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size']=18; plt.rcParams['axes.labelsize']=22
plt.rcParams['xtick.labelsize']=18; plt.rcParams['ytick.labelsize']=18
plt.rcParams['axes.linewidth']=1.6

ROOT=Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3=ROOT/'calculations/G3_adsorption'
OUT=ROOT/'reports/predft_advisor_figures/violin_conv_compare'
OUT.mkdir(parents=True,exist_ok=True)

SURFACES=['S1','S2','S3','S3b','S4']
LBL={'S1':'Pd(100)','S2':'1 ML\nPdO(101)/Pd(100)','S3':'O-rich\nPdO(100)','S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

def load_conv(sid, ads):
    recs = json.load(open(G3/'convB_binding'/f'{sid}_{ads}_convB.json'))
    return [r['E_bind_A_clean'] for r in recs], [r['E_bind_B_frozen'] for r in recs]

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

for ads, ylab, fname in [('CO',r'$E_{\mathrm{bind}}$(CO$^*$) / eV','CO_convA_vs_B.png'),
                          ('CH3O',r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV','CH3O_convA_vs_B.png')]:
    A_list, B_list = [], []
    for s in SURFACES:
        a,b = load_conv(s, ads)
        A_list.append(a); B_list.append(b)
    all_v = [v for d in A_list+B_list for v in d]
    ymin, ymax = min(all_v)-0.2, max(all_v)+0.2
    fig, axes = plt.subplots(1,2,figsize=(20,8), sharey=True)
    violin(axes[0], A_list, 'Convention A\n(clean slab)')
    violin(axes[1], B_list, 'Convention B\n(frozen post-relax slab)')
    axes[0].set_ylabel(ylab)
    for ax in axes: ax.set_ylim(ymin,ymax)
    plt.tight_layout()
    plt.savefig(OUT/fname,dpi=300,bbox_inches='tight'); plt.close()
    print(f'{fname} ✓')

# Print summary stats for S4
print('\n=== S4 stats ===')
for ads in ['CO','CH3O']:
    a,b = load_conv('S4', ads)
    print(f'{ads}  A(clean) min={min(a):.2f} median={np.median(a):.2f}  |  B(frozen) min={min(b):.2f} median={np.median(b):.2f}')
