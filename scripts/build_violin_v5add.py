"""Violin plots with v4 (70) + v5-new (16) picks overlaid distinctly."""
import json, csv, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['font.size']=18
plt.rcParams['axes.labelsize']=22
plt.rcParams['xtick.labelsize']=18
plt.rcParams['ytick.labelsize']=18
plt.rcParams['axes.linewidth']=1.5

ROOT=Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3=ROOT/'calculations/G3_adsorption'
OUT=ROOT/'reports/predft_advisor_figures/violin_v5add'
OUT.mkdir(parents=True, exist_ok=True)

SURFACES_SINGLE=['S1','S2','S3','S3b','S4']
SURFACES_COADS=['S1','S2','S3','S3b']
SDIRS={'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
       'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LBL={'S1':'Pd(100)','S2':'PdO(101)\n/Pd(100)','S3':'O-rich\nPdO(100)',
     'S3b':'Pd-rich\nPdO(100)','S4':r'PdO$_2$(110)'}
CMAP={'S1':'#1f4e79','S2':'#2a9d8f','S3':'#e76f51','S3b':'#f4a261','S4':'#7b2cbf'}

refs=json.load(open(G3/'mace_d3_references.json'))
E_SLAB=refs['slab']; E_CO=refs['gas']['CO']; E_CH3O=refs['gas']['CH3O_radical']

# ---------- load picks: v4 (70) + v5-new (16) ----------
v4=json.load(open(G3/'DFT_shortlist_v3/summary.json'))
manifest=list(csv.DictReader(open(ROOT/'calculations/T1_16_DFT_L2_v5add16/manifest.csv')))
v5new=[r for r in manifest if r['origin']=='v5-new']

v4_by=dict(); v5_by=dict()
for p in v4:
    v4_by.setdefault((p['sid'],p['ads']), []).append((p['idx'], p['E']))
for r in v5new:
    v5_by.setdefault((r['sid'],r['ads']), []).append((int(r['idx']), float(r['E_MLIP']), r['priority']))

def make_plot(ads_key, surfaces, ylabel, fname, title):
    fig, ax = plt.subplots(figsize=(11, 8))
    pos = np.arange(len(surfaces))
    all_data=[]
    for sid in surfaces:
        if ads_key=='coads':
            recs=json.load(open(G3/SDIRS[sid]/'MLIP_phase2_filtered/unique_SetA.json'))
            E_ads_sum=E_CO+E_CH3O
        else:
            recs=json.load(open(G3/SDIRS[sid]/f'MLIP_phase1/unique_{ads_key}.json'))
            E_ads_sum=E_CO if ads_key=='CO' else E_CH3O
        d=[r['E']-E_SLAB[sid]-E_ads_sum for r in recs]
        all_data.append(d)
    # violin + box + scatter background
    parts=ax.violinplot(all_data, positions=pos, widths=0.7,
                        showmeans=False, showmedians=False, showextrema=False)
    for i,pc in enumerate(parts['bodies']):
        pc.set_facecolor(CMAP[surfaces[i]]); pc.set_alpha(0.4); pc.set_edgecolor('black')
    bp=ax.boxplot(all_data, positions=pos, widths=0.28, patch_artist=True,
                  showfliers=False, medianprops={'color':'red','lw':2})
    for i,patch in enumerate(bp['boxes']):
        patch.set_facecolor(CMAP[surfaces[i]]); patch.set_alpha(0.5)
    np.random.seed(42)
    for i,d in enumerate(all_data):
        ax.scatter(np.random.normal(i,0.04,size=len(d)), d, s=3,
                   color='black', alpha=0.15, zorder=2)
    # v4 picks (gold star)
    for i,sid in enumerate(surfaces):
        picks=v4_by.get((sid,ads_key),[])
        if ads_key=='coads': E_ads_sum=E_CO+E_CH3O
        else: E_ads_sum=E_CO if ads_key=='CO' else E_CH3O
        for idx,E in picks:
            Eb=E-E_SLAB[sid]-E_ads_sum
            ax.scatter(i+0.15, Eb, s=170, color='#ffd700', edgecolor='black',
                       linewidth=1.5, zorder=5, marker='*')
    # v5-new picks (red triangle for MUST, orange square for REC, cyan circle for OPTIONAL)
    style={'MUST':('#dc143c','^',200,'MUST'),
           'MUST-diagnostic':('#8b0000','^',200,'MUST-diagnostic'),
           'RECOMMENDED':('#ff8c00','s',170,'RECOMMENDED'),
           'OPTIONAL':('#00ced1','o',150,'OPTIONAL')}
    for i,sid in enumerate(surfaces):
        news=v5_by.get((sid,ads_key),[])
        if ads_key=='coads': E_ads_sum=E_CO+E_CH3O
        else: E_ads_sum=E_CO if ads_key=='CO' else E_CH3O
        for idx,E,prio in news:
            Eb=E-E_SLAB[sid]-E_ads_sum
            c,m,s,_ = style.get(prio, style['OPTIONAL'])
            ax.scatter(i-0.15, Eb, s=s, color=c, edgecolor='black',
                       linewidth=1.5, zorder=6, marker=m)
    ax.set_xticks(pos); ax.set_xticklabels([LBL[s] for s in surfaces], rotation=25, ha='right')
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=18)
    ax.axhline(0, ls='--', color='gray', alpha=0.6, lw=1)
    ax.grid(True, alpha=0.3)
    handles=[
        Line2D([0],[0], marker='*', color='w', markerfacecolor='#ffd700',
               markeredgecolor='black', markersize=15, label='v4 pick (existing, 70)'),
        Line2D([0],[0], marker='^', color='w', markerfacecolor='#dc143c',
               markeredgecolor='black', markersize=12, label='v5-new MUST'),
        Line2D([0],[0], marker='^', color='w', markerfacecolor='#8b0000',
               markeredgecolor='black', markersize=12, label='v5-new MUST-diagnostic'),
        Line2D([0],[0], marker='s', color='w', markerfacecolor='#ff8c00',
               markeredgecolor='black', markersize=12, label='v5-new RECOMMENDED'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#00ced1',
               markeredgecolor='black', markersize=12, label='v5-new OPTIONAL'),
    ]
    ax.legend(handles=handles, loc='best', frameon=True, fontsize=13)
    plt.tight_layout()
    plt.savefig(OUT/fname, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{fname} ✓')

make_plot('CO', SURFACES_SINGLE, r'$E_{\mathrm{bind}}$(CO$^*$) / eV', 'CO_v5add.png',
          'CO* — v4 picks (gold) + v5-new (red/orange/cyan)')
make_plot('CH3O', SURFACES_SINGLE, r'$E_{\mathrm{bind}}$(CH$_3$O$^*$) / eV', 'CH3O_v5add.png',
          'CH₃O* — v4 picks (gold) + v5-new (red/orange/cyan)')
make_plot('coads', SURFACES_COADS, r'$E_{\mathrm{bind}}$(CO$^*$+CH$_3$O$^*$) / eV', 'coads_v5add.png',
          'coads — v4 picks (gold) + v5-new (red/orange/cyan)')
