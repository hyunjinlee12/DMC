"""Co-ads violin plot with Conv B correction. NO outlier deletion (per advisor)."""
import json, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size']=20; plt.rcParams['axes.labelsize']=24
plt.rcParams['xtick.labelsize']=20; plt.rcParams['ytick.labelsize']=20
plt.rcParams['axes.linewidth']=1.8

ROOT=Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3=ROOT/'calculations/G3_adsorption'
OUT=ROOT/'reports/predft_advisor_figures/coads_convB'
OUT.mkdir(parents=True, exist_ok=True)

SURFACES=['S1','S2','S3','S3b','S4']
LBL={'S1':'Pd(100)','S2':'1 ML\nPdO(101)/Pd(100)','S3':'O-rich\nPdO(100)',
     'S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}


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
        ax.scatter(np.random.normal(i,0.04,size=len(d)),d,s=3,color='black',alpha=0.2,zorder=2)
        ax.scatter([i],[min(d)],s=140,color='gold',edgecolor='black',linewidth=2,zorder=5,marker='*')
    ax.set_xticks(pos); ax.set_xticklabels([LBL[s] for s in SURFACES],rotation=30,ha='right')
    ax.set_ylabel(r'$E_{\mathrm{bind}}$ (CO$^*$+CH$_3$O$^*$) / eV')
    ax.axhline(0,ls='--',color='gray',alpha=0.6,lw=1)
    ax.grid(True,alpha=0.3)
    if title:
        ax.set_title(title, fontsize=22, fontweight='bold', pad=10)


# Recompute with CH3O radical ref (consistent with single CH3O)
refs = json.load(open(G3/'mace_d3_references.json'))
E_SLAB = refs['slab']; E_CO = refs['gas']['CO']; E_CH3O_rad = refs['gas']['CH3O_radical']
A_list=[]; B_list=[]
for s in SURFACES:
    recs = json.load(open(G3/'convB_binding'/f'{s}_coads_convB.json'))
    A_list.append([r['E_total'] - E_SLAB[s] - E_CO - E_CH3O_rad for r in recs])
    B_list.append([r['E_total'] - r['E_slab_post'] - E_CO - E_CH3O_rad for r in recs])

# Print stats
print(f"{'Sur':<5} {'n':<6} {'A median':<10} {'B median':<10} {'ΔE_recon':<9}")
for i,s in enumerate(SURFACES):
    print(f"{s:<5} {len(A_list[i]):<6} {np.median(A_list[i]):<+10.2f} {np.median(B_list[i]):<+10.2f} {np.median(A_list[i])-np.median(B_list[i]):<+9.2f}")

# 3 individual plots (no title, keep axes)
for title, data, fname in [
    ('Conv A (clean slab)', A_list, 'coads_convA.png'),
    ('Conv B (frozen post-relax slab)', B_list, 'coads_convB.png'),
]:
    fig, ax = plt.subplots(figsize=(11,8))
    violin(ax, data, None)
    plt.tight_layout()
    plt.savefig(OUT/fname, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{fname} ✓')

# Combined A vs B
fig, axes = plt.subplots(1,2,figsize=(20,8), sharey=True)
all_v=[v for d in A_list+B_list for v in d]
ymin, ymax = min(all_v)-0.2, max(all_v)+0.2
violin(axes[0], A_list, None)
violin(axes[1], B_list, None)
for ax in axes: ax.set_ylim(ymin,ymax)
axes[1].set_ylabel('')
plt.tight_layout()
plt.savefig(OUT/'coads_A_vs_B.png', dpi=300, bbox_inches='tight')
plt.close()
print('coads_A_vs_B ✓')
