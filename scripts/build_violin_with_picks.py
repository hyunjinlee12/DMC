"""Violin plots with top-5 picks overlaid as stars, per surface × ads."""
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
OUT=ROOT/'reports/predft_advisor_figures/violin_with_picks'
OUT.mkdir(parents=True,exist_ok=True)

SURFACES_SINGLE=['S1','S2','S3','S3b','S4']
SURFACES_COADS=['S1','S2','S3b','S3']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
       'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LBL={'S1':'Pd(100)','S2':'1 ML\nPdO(101)/Pd(100)','S3':'O-rich\nPdO(100)',
     'S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs=json.load(open(G3/'mace_d3_references.json'))
E_SLAB=refs['slab']; E_CO=refs['gas']['CO']; E_CH3O=refs['gas']['CH3O_radical']

# Load picks
picks_summary=json.load(open(G3/'DFT_shortlist_v3/summary.json'))
picks_by_key={}  # (sid, ads) -> [rank ordered idx,E]
for p in picks_summary:
    key=(p['sid'], p['ads'])
    picks_by_key.setdefault(key, []).append(p)


def compute_E_all(sid, ads):
    n_ads = 2 if ads=='CO' else (5 if ads=='CH3O' else 7)
    if ads=='coads':
        recs = json.load(open(G3/sdir_map[sid]/'MLIP_phase2_filtered/unique_SetA.json'))
    else:
        recs = json.load(open(G3/sdir_map[sid]/f'MLIP_phase1/unique_{ads}.json'))
    E_ads_ref = E_CO + (E_CH3O if ads=='coads' else 0) if ads=='coads' else (E_CO if ads=='CO' else E_CH3O)
    return {r['idx']: r['E'] - E_SLAB[sid] - E_ads_ref for r in recs}

sdir_map=SDIRS
STAR_COLORS=['#ffd700','#ff8c00','#dc143c','#8b0000','#4b0082']  # 5 ranks


def make_plot(ads_key, surfaces, ylabel, fname):
    fig, ax = plt.subplots(figsize=(11, 8))
    pos = np.arange(len(surfaces))
    all_data=[]
    for sid in surfaces:
        n_ads = 2 if ads_key=='CO' else (5 if ads_key=='CH3O' else 7)
        if ads_key == 'coads':
            recs = json.load(open(G3/sdir_map[sid]/'MLIP_phase2_filtered/unique_SetA.json'))
            E_ads_sum = E_CO + E_CH3O
        else:
            recs = json.load(open(G3/sdir_map[sid]/f'MLIP_phase1/unique_{ads_key}.json'))
            E_ads_sum = E_CO if ads_key=='CO' else E_CH3O
        d=[r['E']-E_SLAB[sid]-E_ads_sum for r in recs]
        all_data.append(d)
    # violin
    parts=ax.violinplot(all_data, positions=pos, widths=0.7,
                       showmeans=False, showmedians=False, showextrema=False)
    for i,pc in enumerate(parts['bodies']):
        pc.set_facecolor(CMAP[surfaces[i]]); pc.set_alpha(0.4); pc.set_edgecolor('black')
    bp=ax.boxplot(all_data, positions=pos, widths=0.3, patch_artist=True,
                  showfliers=False, medianprops={'color':'red','lw':2})
    for i,patch in enumerate(bp['boxes']):
        patch.set_facecolor(CMAP[surfaces[i]]); patch.set_alpha(0.6)
    # scatter all
    np.random.seed(42)
    for i,d in enumerate(all_data):
        ax.scatter(np.random.normal(i,0.04,size=len(d)), d, s=4, color='black', alpha=0.2, zorder=2)
    # Overlay picks (5 stars per surface)
    for i, sid in enumerate(surfaces):
        picks = picks_by_key.get((sid, ads_key), [])
        # Compute E_bind of each pick
        if ads_key == 'coads':
            E_ads_sum = E_CO + E_CH3O
        else:
            E_ads_sum = E_CO if ads_key=='CO' else E_CH3O
        for j, p in enumerate(picks):
            E_bind = p['E'] - E_SLAB[sid] - E_ads_sum
            ax.scatter(i, E_bind, s=180, color=STAR_COLORS[j%5], edgecolor='black',
                       linewidth=1.8, zorder=6, marker='*',
                       label=f'rank {j}' if i==0 else None)
    ax.set_xticks(pos); ax.set_xticklabels([LBL[s] for s in surfaces], rotation=30, ha='right')
    ax.set_ylabel(ylabel)
    ax.axhline(0, ls='--', color='gray', alpha=0.6, lw=1)
    ax.grid(True, alpha=0.3)
    # Legend for rank colors
    from matplotlib.lines import Line2D
    legend_els=[Line2D([0],[0], marker='*', color='w', markerfacecolor=STAR_COLORS[j],
                        markeredgecolor='black', markersize=15, label=f'rank {j}')
                for j in range(5)]
    ax.legend(handles=legend_els, loc='best', frameon=True, fontsize=15)
    plt.tight_layout()
    plt.savefig(OUT/fname, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{fname} ✓')


make_plot('CO', SURFACES_SINGLE, r'$E_{\mathrm{bind}}$(CO$^*$) / eV', 'CO_with_picks.png')
make_plot('CH3O', SURFACES_SINGLE, r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV', 'CH3O_with_picks.png')
make_plot('coads', SURFACES_COADS, r'$E_{\mathrm{bind}}$(CO$^*$+CH$_3$O$^*$) / eV', 'coads_with_picks.png')
