"""Comprehensive MLIP results analysis: Phase 1+2+3 per surface.

Outputs:
  reports/G3/MLIP_analysis/
    A1_phase1_E_ads.png          — relative E_ads CO* and CH3O* per surface
    A2_d_min_distribution.png    — d_min (chemisorption strength) per surface
    A3_site_preference.png       — atop/bridge/hollow share in top-N
    A4_phase2_d_react_E.png      — co-ads E vs d_reactive scatter
    A5_descriptor_preview.png    — E_CO vs E_CH3O scatter (5 surfaces)
    A6_oxidation_trend.png       — E_CO vs surface oxidation (S1→S2→S3/S3b→S4)
    A7_case_classification.png   — preliminary Case A-D placement
    MLIP_results_narrative.md    — written analysis
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ase.io import read

plt.rcParams.update({'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 14})

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
G2 = ROOT / 'calculations/G2_slab'
OUT = ROOT / 'reports/G3/MLIP_analysis'
OUT.mkdir(parents=True, exist_ok=True)

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
LABELS = {'S1':'Pd(100)\nPd⁰','S2':'PdO(101)/Pd\nPd⁰+Pd²⁺',
          'S3':'PdO(100) O-term\nPd²⁺ (O top)','S3b':'PdO(100) Pd-term\nPd²⁺ (Pd top)',
          'S4':'PdO₂(110)\nPd⁴⁺'}
COLORS = {'S1':'#1f4e79', 'S2':'#2a9d8f', 'S3':'#e76f51',
          'S3b':'#f4a261', 'S4':'#7b2cbf'}
# Slab reference E (E_sigma->0)
E_SLAB = {'S1':-434.380, 'S2':-618.565, 'S3':-724.103, 'S3b':-570.772, 'S4':-788.493}

# Load Phase 1 results
def load_phase1():
    data = {}
    for sid in SURFACES:
        d = G3 / SDIRS[sid] / 'MLIP_phase1'
        rco = json.load(open(d / 'unique_CO.json'))
        rch = json.load(open(d / 'unique_CH3O.json'))
        data[sid] = {'CO': rco, 'CH3O': rch}
    return data


# Load Phase 2 filtered
def load_phase2():
    data = {}
    for sid in SURFACES:
        f = G3 / SDIRS[sid] / 'MLIP_phase2_filtered/unique_SetA.json'
        if f.exists():
            data[sid] = json.load(open(f))
        else:
            data[sid] = []
    return data


p1 = load_phase1()
p2 = load_phase2()

# =====================================================================
# Figure A1: relative E_ads (E_top1 - E_slab) per surface
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# Per surface get top-10 relative E (relative to top-1)
for ax, ads, title in [(axes[0], 'CO', 'CO* relative E (MLIP+D3) — top 10 each'),
                       (axes[1], 'CH3O', 'CH₃O* relative E (MLIP+D3) — top 10 each')]:
    for i, sid in enumerate(SURFACES):
        recs = p1[sid][ads]
        if not recs: continue
        E0 = recs[0]['E']
        # All ΔE_rel for top 10
        dE = [(r['E'] - E0) * 1000 for r in recs[:10]]
        x_pos = [i + 0.05 * j for j in range(len(dE))]
        ax.scatter([i]*len(dE), dE, s=60, color=COLORS[sid], alpha=0.7, edgecolor='black')
        # mark top-1
        ax.scatter([i], [0], s=120, color=COLORS[sid], edgecolor='black', linewidth=2, marker='*')
    ax.set_xticks(range(len(SURFACES)))
    ax.set_xticklabels([LABELS[s] for s in SURFACES], fontsize=10)
    ax.set_ylabel(r'$\Delta E$ vs top-1 / meV')
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(OUT / 'A1_phase1_E_ads.png', dpi=140, bbox_inches='tight')
plt.close()

# =====================================================================
# Figure A2: d_min distributions per surface (chemisorption strength)
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
for ax, ads, title in [(axes[0], 'CO', 'CO* d_min (Pd–C) distribution — all unique'),
                       (axes[1], 'CH3O', 'CH₃O* d_min (Pd–O) distribution')]:
    for sid in SURFACES:
        recs = p1[sid][ads]
        d_mins = [r['d_min'] for r in recs]
        ax.hist(d_mins, bins=30, alpha=0.5, color=COLORS[sid], label=f'{sid} (n={len(d_mins)})', edgecolor='black')
    if ads == 'CO':
        ax.axvspan(1.85, 2.10, alpha=0.15, color='green', label='Pd-C chemisorbed band')
    else:
        ax.axvspan(2.00, 2.15, alpha=0.15, color='green', label='Pd-O chemisorbed band')
    ax.axvline(3.0, ls='--', c='red', alpha=0.4, label='physisorption (~3 Å)')
    ax.set_xlabel(r'$d_{\min}$ ads–substrate / Å')
    ax.set_ylabel('Count')
    ax.set_title(title)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT / 'A2_d_min_distribution.png', dpi=140, bbox_inches='tight')
plt.close()

# =====================================================================
# Figure A3: Phase 2 co-ads E vs d_reactive
# =====================================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()
for i, sid in enumerate(SURFACES):
    ax = axes[i]
    recs = p2.get(sid, [])
    if not recs:
        ax.set_visible(False); continue
    E0 = recs[0]['E']
    d_react = [r['d_reactive'] for r in recs]
    dE = [(r['E'] - E0) * 1000 for r in recs]
    ax.scatter(d_react, dE, s=15, alpha=0.4, color=COLORS[sid])
    ax.axvspan(2.1, 4.0, alpha=0.1, color='green', label='SetA band')
    ax.set_xlabel(r'$d(C_{CO} - O_{CH_3O})$ / Å (direct)')
    ax.set_ylabel(r'$\Delta E_{MACE}$ / meV')
    ax.set_title(f'{LABELS[sid]} (n={len(recs)})')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
axes[-1].set_visible(False)
fig.suptitle('Phase 2: co-adsorption E vs reactive distance', fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig(OUT / 'A3_phase2_d_react_E.png', dpi=140, bbox_inches='tight')
plt.close()

# =====================================================================
# Figure A4: Descriptor map preview (top-1 CO vs top-1 CH3O ads E)
# =====================================================================
fig, ax = plt.subplots(figsize=(9, 7))
pts = []
for sid in SURFACES:
    co_top = p1[sid]['CO'][0]['E']
    ch_top = p1[sid]['CH3O'][0]['E']
    E_slab = E_SLAB[sid]
    # ΔE_ads = E_slab+ads - E_slab - μ_ref (we plot just E_slab+ads - E_slab, no μ_ref)
    # This is "absolute binding strength" proxy. Real DFT will add μ_gas correction.
    dE_CO = co_top - E_slab
    dE_CH3O = ch_top - E_slab
    pts.append((sid, dE_CO, dE_CH3O))
    ax.scatter(dE_CO, dE_CH3O, s=300, color=COLORS[sid], edgecolor='black', linewidth=2, zorder=3)
    ax.annotate(f"  {sid}\n  {LABELS[sid].split(chr(10))[0]}", (dE_CO, dE_CH3O),
                xytext=(10, 5), textcoords='offset points', fontsize=10, weight='bold')

# Annotate Case framework regions (qualitative)
ax.axhline(np.mean([p[2] for p in pts]), color='grey', ls=':', alpha=0.5)
ax.axvline(np.mean([p[1] for p in pts]), color='grey', ls=':', alpha=0.5)
ax.text(0.02, 0.98, 'Case A\n(CO 강 + CH3O 강)\nDMC favorable', transform=ax.transAxes,
        fontsize=10, va='top', bbox=dict(boxstyle='round', facecolor='#e8f4f8', alpha=0.8))
ax.text(0.98, 0.02, 'Case D\n(둘 다 약)\nDMC inactive', transform=ax.transAxes,
        fontsize=10, va='bottom', ha='right', bbox=dict(boxstyle='round', facecolor='#fde8e8', alpha=0.8))
ax.text(0.98, 0.98, 'Case C\n(CO 약, CH3O 강)\nside-path', transform=ax.transAxes,
        fontsize=10, va='top', ha='right', bbox=dict(boxstyle='round', facecolor='#fff3cd', alpha=0.8))

ax.set_xlabel(r'$E^{MACE+D3}$(slab+CO*) − $E_{slab}$  /  eV    (more negative = CO* binds stronger)')
ax.set_ylabel(r'$E^{MACE+D3}$(slab+CH₃O*) − $E_{slab}$  /  eV')
ax.set_title('Preliminary descriptor map (MLIP top-1, slab subtracted)\n'
             'Note: μ_gas not subtracted (MLIP-only); axes invert sign convention vs literature')
ax.grid(True, alpha=0.3)
ax.invert_xaxis()
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(OUT / 'A4_descriptor_preview.png', dpi=140, bbox_inches='tight')
plt.close()

# =====================================================================
# Figure A5: oxidation trend — surface vs CO/CH3O d_min top-1
# =====================================================================
oxidation_order = ['S1', 'S2', 'S3b', 'S3', 'S4']  # Pd⁰ → composite → PdO Pd-side → PdO O-side → PdO₂
ox_labels = ['S1\nPd⁰', 'S2\nPd⁰+Pd²⁺', 'S3b\nPd²⁺ Pd-top', 'S3\nPd²⁺ O-top', 'S4\nPd⁴⁺']
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
ax = axes[0]
co_dmin = [p1[sid]['CO'][0]['d_min'] for sid in oxidation_order]
ch_dmin = [p1[sid]['CH3O'][0]['d_min'] for sid in oxidation_order]
x = range(len(oxidation_order))
ax.plot(x, co_dmin, 'o-', markersize=15, lw=2, color='#1f4e79', label='CO* (Pd–C)', markeredgecolor='black')
ax.plot(x, ch_dmin, 's-', markersize=15, lw=2, color='#e76f51', label='CH$_3$O* (Pd–O)', markeredgecolor='black')
for i, (xi, yi) in enumerate(zip(x, co_dmin)):
    ax.annotate(f'{yi:.2f}', (xi, yi), xytext=(5, 5), textcoords='offset points', fontsize=9)
for i, (xi, yi) in enumerate(zip(x, ch_dmin)):
    ax.annotate(f'{yi:.2f}', (xi, yi), xytext=(5, -12), textcoords='offset points', fontsize=9)
ax.axhspan(1.85, 2.15, alpha=0.15, color='green', label='chemisorbed band')
ax.axhline(3.0, ls='--', c='red', alpha=0.4, label='physisorption')
ax.set_xticks(x); ax.set_xticklabels(ox_labels)
ax.set_ylabel(r'top-1 $d_{\min}$ / Å')
ax.set_title('(a) Surface oxidation → chemisorption strength')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Right panel: E_range (discrimination strength) per surface
ax = axes[1]
co_erange = [(p1[sid]['CO'][-1]['E'] - p1[sid]['CO'][0]['E'])*1000 for sid in oxidation_order]
ch_erange = [(p1[sid]['CH3O'][-1]['E'] - p1[sid]['CH3O'][0]['E'])*1000 for sid in oxidation_order]
w = 0.35
ax.bar([i-w/2 for i in x], co_erange, w, color='#1f4e79', edgecolor='black', label='CO*')
ax.bar([i+w/2 for i in x], ch_erange, w, color='#e76f51', edgecolor='black', label='CH$_3$O*')
ax.axhline(200, ls='--', c='red', alpha=0.5, label='discriminative threshold')
ax.set_xticks(x); ax.set_xticklabels(ox_labels)
ax.set_ylabel(r'$E_{range}$ (top vs bottom unique) / meV')
ax.set_title('(b) Energy spread = ranking discrimination')
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(OUT / 'A5_oxidation_trend.png', dpi=140, bbox_inches='tight')
plt.close()

# =====================================================================
# Build summary stats for narrative
# =====================================================================
stats = {}
for sid in SURFACES:
    s = {}
    rco = p1[sid]['CO']; rch = p1[sid]['CH3O']
    s['n_CO_unique'] = len(rco)
    s['n_CH3O_unique'] = len(rch)
    s['CO_top1_E'] = rco[0]['E']
    s['CO_top1_d_min'] = rco[0]['d_min']
    s['CO_E_range_meV'] = (rco[-1]['E'] - rco[0]['E']) * 1000 if len(rco) > 1 else 0
    s['CH3O_top1_E'] = rch[0]['E']
    s['CH3O_top1_d_min'] = rch[0]['d_min']
    s['CH3O_E_range_meV'] = (rch[-1]['E'] - rch[0]['E']) * 1000 if len(rch) > 1 else 0
    coads = p2.get(sid, [])
    s['n_coads_unique'] = len(coads)
    if coads:
        s['coads_top1_E'] = coads[0]['E']
        s['coads_top1_d_react'] = coads[0]['d_reactive']
        s['coads_top1_d_min'] = coads[0]['d_min']
        s['coads_E_range_meV'] = (coads[-1]['E'] - coads[0]['E']) * 1000 if len(coads) > 1 else 0
    else:
        s['coads_top1_E'] = None
    # ΔE_ads vs slab (MLIP top-1)
    s['dE_CO_vs_slab'] = rco[0]['E'] - E_SLAB[sid]
    s['dE_CH3O_vs_slab'] = rch[0]['E'] - E_SLAB[sid]
    stats[sid] = s

json.dump(stats, open(OUT / 'stats_summary.json', 'w'), indent=2)
print('Stats saved.')

# Print summary table
print('\n' + '='*100)
print(f"{'Sur':<5} {'CO_d_min':>9} {'CO_Erng':>9} {'CH3O_d_min':>11} {'CH3O_Erng':>10} {'coads_d_react':>13} {'coads_d_min':>11}")
print('='*100)
for sid in SURFACES:
    s = stats[sid]
    coads_dr = f"{s.get('coads_top1_d_react', 0):.2f}" if s.get('coads_top1_d_react') else '—'
    coads_dm = f"{s.get('coads_top1_d_min', 0):.2f}" if s.get('coads_top1_d_min') else '—'
    print(f"{sid:<5} {s['CO_top1_d_min']:>9.3f} {s['CO_E_range_meV']:>9.0f} {s['CH3O_top1_d_min']:>11.3f} {s['CH3O_E_range_meV']:>10.0f} {coads_dr:>13} {coads_dm:>11}")
print()
print(f'Figures saved to {OUT}')
