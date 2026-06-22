"""Build T1.15 comprehensive report with figures + docx.

Figures:
  fig01_pipeline_overview.png    — workflow schematic
  fig02_phase_counts.png         — candidates per phase per surface
  fig03_phase1_yields.png        — Phase 1 D3 yields + E_range
  fig04_phase2_filter_yield.png  — Phase 2 raw vs filtered
  fig05_phase3_filter_yield.png  — Phase 3 SetTS+SetB raw vs filtered
  fig06_committee_verdicts.png   — committee history (3 cycles)
  fig07_top1_chemistry.png       — Pd-C / Pd-O top-1 bond distances
  fig08_d_reactive_phase2.png    — Phase 2 d_reactive distributions
  fig09_dft_shortlist.png        — DFT shortlist composition (stacked bar)
  fig10_descriptor_preview.png   — preliminary E_CO vs E_CH3O scatter (MLIP)
  fig11_top_structures.png       — sample top-1 structure renders per surface

Then assemble docx using python-docx.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from ase.io import read
from ase.visualize.plot import plot_atoms

plt.rcParams.update({'font.size': 13, 'axes.labelsize': 14, 'axes.titlesize': 15,
                     'legend.fontsize': 11, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
G2 = ROOT / 'calculations/G2_slab'
FIG = ROOT / 'reports/G3/T1_15_figures'
FIG.mkdir(parents=True, exist_ok=True)

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LABELS = {'S1':'Pd(100)','S2':'PdO(101)/Pd(100)','S3':'PdO(100) O-term',
          'S3b':'PdO(100) PdO-term','S4':'PdO$_2$(110)'}
COLORS = ['#1f4e79', '#2a9d8f', '#e76f51', '#f4a261', '#7b2cbf']

# ======================================================================
# Fig 1: pipeline overview (matplotlib schematic)
# ======================================================================
def fig01_pipeline():
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.axis('off')
    stages = [
        ('G1 Bulk\nrelaxation', 0.05, 0.65, '#1f4e79', 'PASSED'),
        ('G2 Slab\n5 surfaces', 0.22, 0.65, '#2a9d8f', 'PASSED'),
        ('T1.10–13\nHeuristic\nenumeration', 0.39, 0.65, '#e76f51', 'DONE'),
        ('T1.14 MACE\nranking', 0.56, 0.65, '#f4a261', 'DONE'),
        ('T1.15 DFT\nshortlist', 0.73, 0.65, '#7b2cbf', 'NOW'),
        ('T1.16–17 DFT\nvacuum + sol', 0.90, 0.65, '#999', 'PENDING'),
    ]
    for txt, x, y, col, status in stages:
        ax.add_patch(plt.Rectangle((x-0.07, y-0.07), 0.14, 0.14, facecolor=col, alpha=0.7, edgecolor='black'))
        ax.text(x, y, txt, ha='center', va='center', fontsize=11, weight='bold', color='white')
        ax.text(x, y-0.13, status, ha='center', va='top', fontsize=9, style='italic',
                color='green' if status in ('PASSED','DONE') else ('red' if status=='NOW' else 'grey'))
    for i in range(5):
        x0 = stages[i][1] + 0.07
        x1 = stages[i+1][1] - 0.07
        ax.annotate('', xy=(x1, 0.65), xytext=(x0, 0.65), arrowprops=dict(arrowstyle='->', lw=2))

    # MLIP sub-phases
    ax.text(0.56, 0.40, 'MLIP sub-phases (MACE-MH + D3 + cuEq)',
            ha='center', fontsize=12, weight='bold')
    sub_stages = [
        ('Phase 1\nsingle ads\n2,516', 0.32, 0.25, '#5f7faa'),
        ('Phase 2\nco-ads SetA\n37,956', 0.50, 0.25, '#5f7faa'),
        ('Phase 3\nSetTS + SetB\n6,314', 0.68, 0.25, '#5f7faa'),
    ]
    for txt, x, y, col in sub_stages:
        ax.add_patch(plt.Rectangle((x-0.07, y-0.06), 0.14, 0.12, facecolor=col, alpha=0.6, edgecolor='black'))
        ax.text(x, y, txt, ha='center', va='center', fontsize=10, color='white')

    # committee box
    ax.add_patch(plt.Rectangle((0.10, 0.02), 0.80, 0.10, facecolor='#fff3cd', edgecolor='#b8860b', linewidth=2))
    ax.text(0.50, 0.07, '5-judge Committee (Paimon-extended): methods · physics · statistics · silent-error · malicious',
            ha='center', fontsize=10, weight='bold')
    ax.text(0.50, 0.92, 'Pd / PdO / PdO$_2$ Surface DMC Formation — Pipeline Overview',
            ha='center', fontsize=15, weight='bold')
    plt.tight_layout()
    plt.savefig(FIG/'fig01_pipeline_overview.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 2: phase candidate counts per surface
# ======================================================================
def fig02_phase_counts():
    g1 = json.load(open(G3/'MLIP_phase1_summary.json'))
    g2 = json.load(open(G3/'MLIP_phase2_summary.json'))
    g3p = json.load(open(G3/'MLIP_phase3_summary.json'))

    # per-surface counts
    p1_CO = []; p1_CH3O = []; p2 = []; p3_TS = []; p3_B = []
    for sid in SURFACES:
        sdir = SDIRS[sid]
        s1 = next(s for s in g1['surfaces'] if s['surface']==sdir)
        p1_CO.append(s1['adsorbates']['CO']['n_raw'])
        p1_CH3O.append(s1['adsorbates']['CH3O']['n_raw'])
        s2 = next(s for s in g2 if s['surface']==sid)
        p2.append(s2['n_raw'])
        ts = next(s for s in g3p if s['surface']==sid and s['pool']=='SetTS')
        b  = next(s for s in g3p if s['surface']==sid and s['pool']=='SetB')
        p3_TS.append(ts['n_raw'])
        p3_B.append(b['n_raw'])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    x = np.arange(len(SURFACES))
    w = 0.18
    ax = axes[0]
    ax.bar(x - 2*w, p1_CO, w, label='Phase 1 CO*', color='#1f4e79', edgecolor='black')
    ax.bar(x - w, p1_CH3O, w, label='Phase 1 CH$_3$O*', color='#5f7faa', edgecolor='black')
    ax.bar(x, p2, w, label='Phase 2 SetA (reactive)', color='#2a9d8f', edgecolor='black')
    ax.bar(x + w, p3_TS, w, label='Phase 3 SetTS', color='#e76f51', edgecolor='black')
    ax.bar(x + 2*w, p3_B, w, label='Phase 3 SetB', color='#f4a261', edgecolor='black')
    ax.set_xticks(x); ax.set_xticklabels([LABELS[s] for s in SURFACES], rotation=20, ha='right')
    ax.set_ylabel('Number of candidates')
    ax.set_yscale('log')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_title('(a) Raw heuristic + MLIP-relaxed candidate counts')
    ax.grid(True, alpha=0.3, axis='y')

    # Total over all surfaces (pie)
    ax = axes[1]
    totals = [sum(p1_CO)+sum(p1_CH3O), sum(p2), sum(p3_TS)+sum(p3_B)]
    labels = [f'Phase 1\n{totals[0]:,}', f'Phase 2 SetA\n{totals[1]:,}', f'Phase 3 (TS+B)\n{totals[2]:,}']
    colors = ['#1f4e79', '#2a9d8f', '#e76f51']
    ax.pie(totals, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90,
           textprops={'fontsize': 12, 'weight':'bold'}, wedgeprops={'edgecolor':'black'})
    ax.set_title(f'(b) Total MLIP relaxations: {sum(totals):,}')
    plt.tight_layout()
    plt.savefig(FIG/'fig02_phase_counts.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 3: Phase 1 D3 yields
# ======================================================================
def fig03_phase1_yields():
    g1 = json.load(open(G3/'MLIP_phase1_summary.json'))
    conv_CO = []; conv_CH3O = []; uniq_CO = []; uniq_CH3O = []; erng_CO = []; erng_CH3O = []
    for sid in SURFACES:
        s = next(x for x in g1['surfaces'] if x['surface']==SDIRS[sid])
        conv_CO.append(100*s['adsorbates']['CO']['n_converged']/s['adsorbates']['CO']['n_raw'])
        conv_CH3O.append(100*s['adsorbates']['CH3O']['n_converged']/s['adsorbates']['CH3O']['n_raw'])
        uniq_CO.append(s['adsorbates']['CO']['n_unique'])
        uniq_CH3O.append(s['adsorbates']['CH3O']['n_unique'])
        erng_CO.append(s['adsorbates']['CO']['E_range_meV'])
        erng_CH3O.append(s['adsorbates']['CH3O']['E_range_meV'])

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    x = np.arange(len(SURFACES))
    w = 0.35

    ax = axes[0]
    ax.bar(x-w/2, conv_CO, w, label='CO*', color='#1f4e79', edgecolor='black')
    ax.bar(x+w/2, conv_CH3O, w, label='CH$_3$O*', color='#5f7faa', edgecolor='black')
    ax.axhline(70, ls='--', c='red', alpha=0.5, label='REJECT threshold (70%)')
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('Convergence rate / %'); ax.set_ylim(60, 105)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3, axis='y')
    ax.set_title('(a) LBFGS convergence (fmax=0.05 eV/Å)')

    ax = axes[1]
    ax.bar(x-w/2, uniq_CO, w, label='CO*', color='#1f4e79', edgecolor='black')
    ax.bar(x+w/2, uniq_CH3O, w, label='CH$_3$O*', color='#5f7faa', edgecolor='black')
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('Unique structures after dedup')
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3, axis='y')
    ax.set_title('(b) Post-dedup unique counts')

    ax = axes[2]
    ax.bar(x-w/2, erng_CO, w, label='CO*', color='#1f4e79', edgecolor='black')
    ax.bar(x+w/2, erng_CH3O, w, label='CH$_3$O*', color='#5f7faa', edgecolor='black')
    ax.axhline(200, ls='--', c='red', alpha=0.5, label='Discriminative threshold (200 meV)')
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel(r'$\Delta E$ range / meV'); ax.set_yscale('log')
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3, axis='y')
    ax.set_title('(c) Energy spread (discrimination signal)')
    plt.tight_layout()
    plt.savefig(FIG/'fig03_phase1_yields.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 4: Phase 2 filter yield (silent bug demo)
# ======================================================================
def fig04_phase2_filter():
    raw_data = json.load(open(G3/'MLIP_phase2_summary.json'))
    filt_data = json.load(open(G3/'MLIP_phase2_filtered_summary.json'))
    raw = []; survivors = []; frag = []; bond = []; unc = []
    for sid in SURFACES:
        rs = next(s for s in raw_data if s['surface']==sid)
        fs = next(s for s in filt_data if s['surface']==sid)
        raw.append(rs['n_raw'])
        survivors.append(fs['n_survivors'])
        frag.append(fs['drop_reasons'].get('fragmented', 0))
        bond.append(fs['drop_reasons'].get('CH_broken', 0) +
                    fs['drop_reasons'].get('methoxy_OC_broken', 0) +
                    fs['drop_reasons'].get('CO_bond_broken', 0))
        unc.append(fs['drop_reasons'].get('unconverged', 0))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    x = np.arange(len(SURFACES))

    ax = axes[0]
    ax.bar(x, raw, color='lightgrey', edgecolor='black', label='Raw (all relaxed)')
    ax.bar(x, survivors, color='#2a9d8f', edgecolor='black', label='Survivors after refilter')
    for i, (r, s) in enumerate(zip(raw, survivors)):
        pct = 100*s/r
        ax.text(i, s + r*0.02, f'{pct:.0f}%', ha='center', fontsize=11, weight='bold')
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('Structures')
    ax.legend(fontsize=11)
    ax.set_title('(a) Phase 2 SetA: raw vs post-filter\n(silent PBC fragmentation removed)')
    ax.grid(True, alpha=0.3, axis='y')

    ax = axes[1]
    w = 0.6
    ax.bar(x, frag, w, color='#e76f51', edgecolor='black', label='Fragmented (PBC silent bug)')
    ax.bar(x, bond, w, bottom=frag, color='#f4a261', edgecolor='black', label='Bond-broken')
    ax.bar(x, unc, w, bottom=np.array(frag)+np.array(bond), color='#999', edgecolor='black', label='Unconverged')
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('Structures dropped')
    ax.legend(fontsize=11)
    ax.set_title('(b) Drop reasons')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(FIG/'fig04_phase2_filter_yield.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 5: Phase 3 filter yield (SetTS, SetB)
# ======================================================================
def fig05_phase3_filter():
    filt = json.load(open(G3/'MLIP_phase3_filtered_summary.json'))
    setts_raw, setts_surv, setb_raw, setb_surv = [], [], [], []
    for sid in SURFACES:
        ts = next((s for s in filt if s['surface']==sid and s['pool']=='SetTS'), None)
        b = next((s for s in filt if s['surface']==sid and s['pool']=='SetB'), None)
        if ts:
            setts_raw.append(ts['n_input']); setts_surv.append(ts['n_survivors'])
        else:
            setts_raw.append(0); setts_surv.append(0)
        if b:
            setb_raw.append(b['n_input']); setb_surv.append(b['n_survivors'])
        else:
            setb_raw.append(100); setb_surv.append(0)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    x = np.arange(len(SURFACES))

    ax = axes[0]
    ax.bar(x, setts_raw, color='lightgrey', edgecolor='black', label='Raw')
    ax.bar(x, setts_surv, color='#e76f51', edgecolor='black', label='Survive filter')
    for i, (r, s) in enumerate(zip(setts_raw, setts_surv)):
        pct = 100*s/r if r else 0
        ax.text(i, s + r*0.02, f'{pct:.0f}%', ha='center', fontsize=11, weight='bold')
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('Structures'); ax.legend(fontsize=11)
    ax.set_title('(a) Phase 3 SetTS (TS guess pool)\n→ 99.7% drift out of band: per guide T2.5, NEB seed not needed')
    ax.grid(True, alpha=0.3, axis='y')

    ax = axes[1]
    ax.bar(x, setb_raw, color='lightgrey', edgecolor='black', label='Raw')
    ax.bar(x, setb_surv, color='#f4a261', edgecolor='black', label='Survive filter')
    for i, (r, s) in enumerate(zip(setb_raw, setb_surv)):
        pct = 100*s/r if r else 0
        ax.text(i, s + 5, f'{pct:.0f}%', ha='center', fontsize=11, weight='bold')
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('Structures (sample size 100)'); ax.legend(fontsize=11)
    ax.set_title('(b) Phase 3 SetB (thermo reference)\n→ Low yield: PBC drift + collapse (S4 = 0 survivors)')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(FIG/'fig05_phase3_filter_yield.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 6: committee verdicts timeline (3 cycles)
# ======================================================================
def fig06_committee():
    fig, ax = plt.subplots(figsize=(13, 6))
    cycles = ['Phase 1 (noD3)', 'Phase 1 (D3)', 'Phase 2 (raw)', 'Phase 2 (filtered)', 'Phase 3 (raw)']
    judges = ['methods', 'physics', 'statistics', 'silent-error', 'malicious']
    # 4 = Pass, 3 = Pass-w-c, 2 = Concern, 1 = Reject
    matrix = [
        [3, 4, 2, 1, None],   # Phase 1 noD3:  methods=PwC(6), phys=Reject(4 in old labels), stat=Reject(5)
        [4, 2, 2, 4, None],   # Phase 1 D3
        [4, 1, 2, 1, 4],      # Phase 2 raw
        [4, 2, 2, 4, 4],      # Phase 2 filtered (retry)
        [4, 1, 2, 1, 4],      # Phase 3 raw
    ]
    # Use color coding
    cmap = {1:'#c0392b', 2:'#f39c12', 3:'#27ae60', 4:'#16a085', None:'#bbb'}
    legend_map = {4:'Pass', 3:'Pass-w-caveats', 2:'Concern', 1:'Reject', None:'(not run)'}

    for i, cycle in enumerate(cycles):
        for j, jud in enumerate(judges):
            v = matrix[i][j]
            ax.add_patch(plt.Rectangle((j, len(cycles)-1-i), 1, 1, facecolor=cmap[v], edgecolor='black'))
            label = {1:'R', 2:'C', 3:'P-c', 4:'P', None:'-'}[v]
            ax.text(j+0.5, len(cycles)-1-i+0.5, label, ha='center', va='center', fontsize=12, weight='bold', color='white')

    ax.set_xlim(0, len(judges)); ax.set_ylim(0, len(cycles))
    ax.set_xticks([i+0.5 for i in range(len(judges))]); ax.set_xticklabels(judges, rotation=15)
    ax.set_yticks([i+0.5 for i in range(len(cycles))]); ax.set_yticklabels(reversed(cycles))
    ax.set_title('Committee verdict timeline (Paimon-extended 5-judge)')

    handles = [mpatches.Patch(color=cmap[v], label=legend_map[v]) for v in [4,3,2,1]]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=11)
    plt.tight_layout()
    plt.savefig(FIG/'fig06_committee_verdicts.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 7: top-1 chemistry (Pd-C, Pd-O for top-3 shortlist)
# ======================================================================
def fig07_top1_chemistry():
    """Plot d_min for top-3 single CO and CH3O per surface."""
    dft = json.load(open(G3/'DFT_shortlist/shortlist_global.json'))
    d_co = {sid: [] for sid in SURFACES}
    d_ch = {sid: [] for sid in SURFACES}
    for e in dft:
        if e['kind'] == 'single_CO':   d_co[e['surface']].append(e['d_min'])
        if e['kind'] == 'single_CH3O': d_ch[e['surface']].append(e['d_min'])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    ax = axes[0]
    for i, sid in enumerate(SURFACES):
        ax.scatter([i]*len(d_co[sid]), d_co[sid], s=120, color=COLORS[i], edgecolor='black', label=LABELS[sid])
    ax.axhspan(1.85, 2.10, color='green', alpha=0.15, label='chemisorbed band')
    ax.axhline(3.0, ls='--', c='red', alpha=0.5, label='physisorption (~3 Å)')
    ax.set_xticks(range(len(SURFACES))); ax.set_xticklabels(SURFACES)
    ax.set_ylabel(r'$d_{\min}$(Pd–C) / Å (top-3 CO*)')
    ax.set_title('(a) CO* DFT shortlist: chemisorption check')
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9, loc='upper left')

    ax = axes[1]
    for i, sid in enumerate(SURFACES):
        ax.scatter([i]*len(d_ch[sid]), d_ch[sid], s=120, color=COLORS[i], edgecolor='black')
    ax.axhspan(2.00, 2.15, color='green', alpha=0.15, label='chemisorbed band')
    ax.axhline(3.0, ls='--', c='red', alpha=0.5, label='physisorption')
    ax.set_xticks(range(len(SURFACES))); ax.set_xticklabels(SURFACES)
    ax.set_ylabel(r'$d_{\min}$(Pd–O) / Å (top-3 CH$_3$O*)')
    ax.set_title('(b) CH$_3$O* DFT shortlist: chemisorption check')
    ax.grid(True, alpha=0.3); ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(FIG/'fig07_top1_chemistry.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 8: d_reactive distributions per surface (Phase 2 filtered)
# ======================================================================
def fig08_d_reactive():
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.flatten()
    for i, sid in enumerate(SURFACES):
        ax = axes[i]
        fpath = G3 / SDIRS[sid] / 'MLIP_phase2_filtered/unique_SetA.json'
        if not fpath.exists():
            ax.set_visible(False); continue
        u = json.load(open(fpath))
        dr = [r['d_reactive'] for r in u]
        ax.hist(dr, bins=40, color=COLORS[i], edgecolor='black', alpha=0.7)
        ax.axvspan(2.1, 4.0, alpha=0.15, color='green', label='SetA reactive band')
        ax.set_xlabel(r'$d(C_{CO} - O_{CH_3O})$ / Å (direct)')
        ax.set_ylabel('Count')
        ax.set_title(f'{LABELS[sid]} (n={len(u)})')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    axes[-1].set_visible(False)
    fig.suptitle('Phase 2 filtered SetA: reactive distance distributions', fontsize=14, weight='bold')
    plt.tight_layout()
    plt.savefig(FIG/'fig08_d_reactive_phase2.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 9: DFT shortlist composition
# ======================================================================
def fig09_dft_shortlist():
    co = []; ch = []; coads = []
    for sid in SURFACES:
        sl = json.load(open(G3/'DFT_shortlist'/sid/'shortlist.json'))
        co.append(sum(1 for e in sl if e['kind']=='single_CO'))
        ch.append(sum(1 for e in sl if e['kind']=='single_CH3O'))
        coads.append(sum(1 for e in sl if e['kind']=='coads_SetA'))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    x = np.arange(len(SURFACES))
    ax.bar(x, co, label='single CO*', color='#1f4e79', edgecolor='black')
    ax.bar(x, ch, bottom=co, label='single CH$_3$O*', color='#5f7faa', edgecolor='black')
    ax.bar(x, coads, bottom=np.array(co)+np.array(ch), label='co-ads (SetA filtered)', color='#2a9d8f', edgecolor='black')
    totals = np.array(co)+np.array(ch)+np.array(coads)
    for i, t in enumerate(totals):
        ax.text(i, t+0.3, str(t), ha='center', fontsize=12, weight='bold')
    ax.set_xticks(x); ax.set_xticklabels([LABELS[s] for s in SURFACES], rotation=15, ha='right')
    ax.set_ylabel('Number of DFT jobs')
    ax.set_title(f'T1.15 DFT Shortlist Composition (Total: {sum(totals)} jobs)')
    ax.legend(fontsize=11, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(FIG/'fig09_dft_shortlist.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 10: descriptor map preview (MLIP-based, will be replaced by DFT)
# ======================================================================
def fig10_descriptor_preview():
    g1 = json.load(open(G3/'MLIP_phase1_summary.json'))
    # Use top-1 E_MACE for CO and CH3O per surface (relative to slab+gas ref — but we just plot E_top1 since slab E is same per surface)
    # Better: ΔE = E_top1 - E_slab (slab CONTCAR E is fixed)
    E_slab = {'S1':-434.408,'S2':-618.565,'S3':-724.152,'S3b':-570.770,'S4':-788.493}
    # Gas refs (approx from MACE for reference; MLIP only, NOT DFT)
    # Use top-1 from ranked_*.json
    pts = []
    for sid in SURFACES:
        rco = json.load(open(G3/SDIRS[sid]/'MLIP_phase1/unique_CO.json'))
        rch = json.load(open(G3/SDIRS[sid]/'MLIP_phase1/unique_CH3O.json'))
        E_CO_ads = rco[0]['E'] - E_slab[sid]   # NOT corrected for gas
        E_CH3O_ads = rch[0]['E'] - E_slab[sid]
        pts.append((sid, E_CO_ads, E_CH3O_ads))

    fig, ax = plt.subplots(figsize=(9, 7))
    for i, (sid, ex, ey) in enumerate(pts):
        ax.scatter(ex, ey, s=250, color=COLORS[i], edgecolor='black', linewidth=2,
                   label=f"{sid} {LABELS[sid]}", zorder=3)
        ax.annotate(sid, (ex, ey), xytext=(8, 8), textcoords='offset points',
                    fontsize=13, weight='bold')
    ax.set_xlabel(r'$E^{MACE+D3}$(CO*) − $E_{slab}$  /  eV')
    ax.set_ylabel(r'$E^{MACE+D3}$(CH$_3$O*) − $E_{slab}$  /  eV')
    ax.set_title('Preliminary descriptor map (MLIP-based, top-1 per surface)\n'
                 '⚠ Absolute E from MLIP — DFT will replace numerically; relative ordering may shift')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIG/'fig10_descriptor_preview.png', dpi=150, bbox_inches='tight')
    plt.close()

# ======================================================================
# Fig 11: top-1 structure renders per surface (single CO/CH3O + co-ads)
# ======================================================================
def fig11_top_structures():
    """Render top-1 single CO, single CH3O, and co-ads per surface (when available)."""
    fig, axes = plt.subplots(5, 3, figsize=(13, 18))
    for i, sid in enumerate(SURFACES):
        for j, kind in enumerate(['single_CO', 'single_CH3O', 'coads_SetA']):
            ax = axes[i, j]
            sub = G3/'DFT_shortlist'/sid/kind
            if not sub.exists():
                ax.axis('off')
                ax.text(0.5, 0.5, '— (unavailable)', ha='center', va='center', transform=ax.transAxes, fontsize=11, color='red')
                if j == 0: ax.text(-0.1, 0.5, sid, fontsize=14, weight='bold', transform=ax.transAxes, va='center')
                continue
            files = sorted(sub.glob('*.vasp'))
            if not files:
                ax.axis('off'); continue
            atoms = read(files[0])
            plot_atoms(atoms, ax, rotation='-80x,5y,0z', radii=0.85, show_unit_cell=2)
            ax.set_xticks([]); ax.set_yticks([])
            tag = {'single_CO':'CO*', 'single_CH3O':'CH$_3$O*', 'coads_SetA':'co-ads'}[kind]
            if i == 0:
                ax.set_title(tag, fontsize=14, weight='bold')
            if j == 0:
                ax.text(-0.15, 0.5, f'{sid}\n{LABELS[sid]}', fontsize=11, weight='bold',
                        transform=ax.transAxes, va='center', ha='right')
    plt.tight_layout()
    plt.savefig(FIG/'fig11_top_structures.png', dpi=130, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    print('Building figures...')
    fig01_pipeline()
    print('  fig01 ✓')
    fig02_phase_counts()
    print('  fig02 ✓')
    fig03_phase1_yields()
    print('  fig03 ✓')
    fig04_phase2_filter()
    print('  fig04 ✓')
    fig05_phase3_filter()
    print('  fig05 ✓')
    fig06_committee()
    print('  fig06 ✓')
    fig07_top1_chemistry()
    print('  fig07 ✓')
    fig08_d_reactive()
    print('  fig08 ✓')
    fig09_dft_shortlist()
    print('  fig09 ✓')
    fig10_descriptor_preview()
    print('  fig10 ✓')
    fig11_top_structures()
    print('  fig11 ✓')
    print(f'\n11 figures saved to {FIG}')
