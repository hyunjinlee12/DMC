"""v5 refinements per 2026-07-17 researcher spec:
- Correct v4 S1/CO characterization (1 atop_Pd + 4 br_PdPd, not 5 br_PdPd).
- Distinguish raw / converged / unique / valid pool counts and explain drops.
- Enrich proposed_additions_only.csv: E_MLIP, ΔE, site, region, coord,
  d_reactive, d_bin, xy dist to nearest v4 pick, selection_reason, priority.
- Attach anchor-neighbor distance tables for MUST + RECOMMENDED adds.
- Priorities per spec:
    MUST:   S3/coads 3481, 5161; S3b/coads 2051, 2754.
    RECOMMENDED: S1/coads 487, 4021; S2/coads 8079, 3633.
    OPTIONAL: single-ads v5-new.
    REVIEW-NEEDED: other coads v5-new that user did not name.

No structure files copied, no calc dirs created.
"""
import csv, json, math
from pathlib import Path
from collections import defaultdict
import numpy as np
from ase.io import read

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT/'calculations/G2_slab'
G3 = ROOT/'calculations/G3_adsorption'
V5 = G3/'DFT_shortlist_v5'
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

MUST = {('S3','coads',3481),('S3','coads',5161),
        ('S3b','coads',2051),('S3b','coads',2754)}
REC  = {('S1','coads',487),('S1','coads',4021),
        ('S2','coads',8079),('S2','coads',3633)}

CUTOFF = 2.6
NBR_RADIUS = 3.5   # for the descriptor table
CO2_CUT = 2.0

# ---------- helpers (reused from v5) ----------
def valid_CO(a):
    syms=a.get_chemical_symbols(); c_i=[i for i,s in enumerate(syms) if s=='C']
    o_i=[i for i,s in enumerate(syms) if s=='O']
    if len(c_i)!=1: return None
    c=c_i[0]
    d_o=sorted([(a.get_distance(c,oi,mic=True),oi) for oi in o_i])
    if not(1.05<=d_o[0][0]<=1.30): return None
    if len(d_o)>=2 and d_o[1][0]<CO2_CUT: return None
    return (c, d_o[0][1])

def valid_CH3O(a):
    syms=a.get_chemical_symbols()
    c_i=[i for i,s in enumerate(syms) if s=='C']
    h_i=[i for i,s in enumerate(syms) if s=='H']
    o_i=[i for i,s in enumerate(syms) if s=='O']
    if len(c_i)!=1 or len(h_i)!=3: return None
    c=c_i[0]
    if not all(0.90<=a.get_distance(c,h,mic=True)<=1.25 for h in h_i): return None
    d_o=sorted([(a.get_distance(c,oi,mic=True),oi) for oi in o_i])
    if not(1.30<=d_o[0][0]<=1.55): return None
    if len(d_o)>=2 and d_o[1][0]<CO2_CUT: return None
    return (c, d_o[0][1])

def valid_coads(a):
    syms=a.get_chemical_symbols()
    c_i=[i for i,s in enumerate(syms) if s=='C']
    h_i=[i for i,s in enumerate(syms) if s=='H']
    o_i=[i for i,s in enumerate(syms) if s=='O']
    if len(c_i)!=2 or len(h_i)!=3: return None
    me_c,co_c=None,None
    for c in c_i:
        nh=sum(1 for h in h_i if a.get_distance(c,h,mic=True)<1.3)
        if nh==3: me_c=c
        else: co_c=c
    if me_c is None or co_c is None: return None
    d_co=sorted([(a.get_distance(co_c,oi,mic=True),oi) for oi in o_i])
    if not(1.05<=d_co[0][0]<=1.30): return None
    if len(d_co)>=2 and d_co[1][0]<CO2_CUT: return None
    co_o=d_co[0][1]
    d_me=sorted([(a.get_distance(me_c,oi,mic=True),oi) for oi in o_i if oi!=co_o])
    if not d_me or not(1.30<=d_me[0][0]<=1.55): return None
    if len(d_me)>=2 and d_me[1][0]<CO2_CUT: return None
    return (co_c, co_o, me_c, d_me[0][1])

def xy_mic(a, p1, p2):
    cell=a.cell.array; inv=np.linalg.inv(cell)
    v=inv@np.array([p2[0]-p1[0], p2[1]-p1[1], 0.0])
    v[0]-=round(v[0]); v[1]-=round(v[1]); v[2]=0
    w=cell.T@v
    return math.hypot(w[0], w[1])

def site_label(a, anchor, slab_idx):
    syms=a.get_chemical_symbols()
    d=a.get_distances(anchor, slab_idx, mic=True)
    nbrs=[(slab_idx[i],d[i]) for i in range(len(slab_idx)) if d[i]<CUTOFF]
    n_pd=sum(1 for i,_ in nbrs if syms[i]=='Pd')
    n_o=sum(1 for i,_ in nbrs if syms[i]=='O')
    tot=n_pd+n_o
    if tot==0: return 'physi', n_pd, n_o
    if tot==1: return ('atop_Pd' if n_pd else 'atop_O'), n_pd, n_o
    if tot==2: return ('br_PdPd' if n_pd==2 else ('br_OO' if n_o==2 else 'br_PdO')), n_pd, n_o
    if tot==3: return ('h3Pd' if n_pd==3 else ('h3O' if n_o==3 else f'h3_{n_pd}Pd{n_o}O')), n_pd, n_o
    return f'{tot}f_{n_pd}Pd{n_o}O', n_pd, n_o

def surface_region(site, n_pd, n_o):
    if site=='physi': return 'physi'
    if n_pd>0 and n_o==0: return 'metal'
    if n_pd==0 and n_o>0: return 'oxide'
    return 'interface'

def neighbor_table(a, anchor, slab_idx, radius=NBR_RADIUS):
    """Return sorted list of (species, dist_A) for slab atoms within radius of anchor."""
    syms=a.get_chemical_symbols()
    d=a.get_distances(anchor, slab_idx, mic=True)
    tbl=[(syms[slab_idx[i]], float(d[i])) for i in range(len(slab_idx)) if d[i]<radius]
    tbl.sort(key=lambda t: t[1])
    return tbl

# ---------- pool count clarification ----------
def raw_pool_size(sid, ads):
    """Count raw AutoAdsorbate candidates (before MLIP). Read from candidates.traj if present."""
    sdir=SDIRS[sid]
    if ads=='coads':
        p=G3/sdir/'coads_guide/SetA.traj'
    else:
        p=G3/sdir/f'{ads}/candidates.traj'
    if not p.exists(): return None
    from ase.io import read as _r
    return len(list(_r(p, index=':')))

def pool_diagnostics():
    """For each group, count raw / MLIP-run / converged / unique / valid_after_filter."""
    rows=[]
    for sid in ['S1','S2','S3','S3b','S4']:
        for ads in ['CO','CH3O','coads']:
            if ads=='coads' and sid=='S4': continue
            sdir=SDIRS[sid]
            raw = raw_pool_size(sid, ads)
            # unique JSON
            if ads=='coads':
                uj = G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'
            else:
                uj = G3/sdir/f'MLIP_phase1/unique_{ads}.json'
            if not uj.exists():
                rows.append({'sid':sid,'ads':ads,'raw':raw or '?','unique_json':0,
                             'converged':0,'valid':0}); continue
            uni = json.load(open(uj))
            conv = [r for r in uni if r.get('converged',True)]
            # valid check: reuse from v5 characterize (load atoms and check)
            if ads=='coads':
                traj=list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
                vf = valid_coads
            else:
                traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
                vf = valid_CO if ads=='CO' else valid_CH3O
            slab_ref=read(G2/sdir/'CONTCAR'); n_ref=len(slab_ref)
            n_ads = 2 if ads=='CO' else 5 if ads=='CH3O' else 7
            valid_count=0
            for r in conv:
                a=traj[r['idx']]
                if len(a)!=n_ref+n_ads:
                    ads_at=a[-n_ads:]; a=slab_ref.copy(); a+=ads_at
                if vf(a) is not None: valid_count+=1
            rows.append({'sid':sid,'ads':ads,'raw':raw if raw is not None else '?',
                         'unique_json':len(uni),'converged':len(conv),'valid':valid_count})
    return rows

# ---------- v4 KEEP xy positions (for xy_dist_to_nearest_v4) ----------
def v4_pick_xy(sid, ads):
    """Return list of (anchor_xy_tuples) for v4 picks of this group."""
    v4=json.load(open(G3/'DFT_shortlist_v3/summary.json'))
    sdir=SDIRS[sid]
    if ads=='coads':
        traj=list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
    else:
        traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
    slab_ref=read(G2/sdir/'CONTCAR'); n_ref=len(slab_ref)
    n_ads = 2 if ads=='CO' else 5 if ads=='CH3O' else 7
    picks=[]
    for r in v4:
        if r['sid']!=sid or r['ads']!=ads: continue
        a=traj[r['idx']]
        if len(a)!=n_ref+n_ads:
            ads_at=a[-n_ads:]; a=slab_ref.copy(); a+=ads_at
        if ads=='CO':
            v=valid_CO(a); anchor=v[0]  # C
        elif ads=='CH3O':
            v=valid_CH3O(a); anchor=v[1]  # O
        else:
            v=valid_coads(a); anchor=v[0]  # co_C
        picks.append((r['idx'], a, anchor, tuple(a.positions[anchor,:2])))
    return picks

# ---------- addition characterization ----------
def build_addition_row(sid, ads, idx):
    sdir=SDIRS[sid]
    if ads=='coads':
        uj=json.load(open(G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'))
        traj=list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
    else:
        uj=json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
        traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
    ur = next(r for r in uj if r['idx']==idx)
    E = float(ur['E'])
    # global min in this group (lowest E from valid+converged pool)
    slab_ref=read(G2/sdir/'CONTCAR'); n_ref=len(slab_ref)
    n_ads = 2 if ads=='CO' else 5 if ads=='CH3O' else 7
    vf = valid_CO if ads=='CO' else valid_CH3O if ads=='CH3O' else valid_coads
    valid_pool=[]
    for r in uj:
        if not r.get('converged',True): continue
        a=traj[r['idx']]
        if len(a)!=n_ref+n_ads:
            ads_at=a[-n_ads:]; a=slab_ref.copy(); a+=ads_at
        if vf(a) is not None: valid_pool.append((r['idx'], float(r['E'])))
    E_min = min(e for _,e in valid_pool)
    dE = E - E_min

    a=traj[idx]
    if len(a)!=n_ref+n_ads:
        ads_at=a[-n_ads:]; a=slab_ref.copy(); a+=ads_at
    v=vf(a)
    if ads=='CO':
        c, o_ads = v
        ads_set={c,o_ads}
        anchor=c; anchor_role='C'
        slab_idx=[i for i in range(len(a)) if i not in ads_set]
        site,npd,no=site_label(a, anchor, slab_idx)
        d_min=float(min(a.get_distances(anchor, slab_idx, mic=True)))
        nbrs=neighbor_table(a, anchor, slab_idx)
        row={'sid':sid,'ads':ads,'idx':idx,'E_MLIP':round(E,4),'dE_from_global':round(dE,4),
             'site':site,'region':surface_region(site,npd,no),'n_Pd_nbr':npd,'n_O_nbr':no,
             'd_anchor_surf':round(d_min,3),
             'd_CO_bond':round(a.get_distance(c,o_ads,mic=True),3),
             'd_CH3O_bond':'','d_reactive':'','d_bin':'','site_combo':'',
             'neighbors_within_3.5A':';'.join(f'{s}:{d:.2f}' for s,d in nbrs)}
    elif ads=='CH3O':
        c, o_bonded = v
        h_i=[i for i,s in enumerate(a.get_chemical_symbols()) if s=='H']
        ads_set={c,o_bonded,*h_i}
        anchor=o_bonded; anchor_role='O'
        slab_idx=[i for i in range(len(a)) if i not in ads_set]
        site,npd,no=site_label(a, anchor, slab_idx)
        d_min=float(min(a.get_distances(anchor, slab_idx, mic=True)))
        nbrs=neighbor_table(a, anchor, slab_idx)
        row={'sid':sid,'ads':ads,'idx':idx,'E_MLIP':round(E,4),'dE_from_global':round(dE,4),
             'site':site,'region':surface_region(site,npd,no),'n_Pd_nbr':npd,'n_O_nbr':no,
             'd_anchor_surf':round(d_min,3),
             'd_CO_bond':'','d_CH3O_bond':round(a.get_distance(c,o_bonded,mic=True),3),
             'd_reactive':'','d_bin':'','site_combo':'',
             'neighbors_within_3.5A':';'.join(f'{s}:{d:.2f}' for s,d in nbrs)}
    else:
        co_c, co_o, me_c, me_o = v
        h_i=[i for i,s in enumerate(a.get_chemical_symbols()) if s=='H']
        ads_set={co_c,co_o,me_c,me_o,*h_i}
        slab_idx=[i for i in range(len(a)) if i not in ads_set]
        s_co,npd_co,no_co=site_label(a, co_c, slab_idx)
        s_me,npd_me,no_me=site_label(a, me_o, slab_idx)
        d_react=float(a.get_distance(co_c, me_o, mic=True))
        d_min=float(min(a.get_distances(co_c, slab_idx, mic=True)))
        nbrs_co=neighbor_table(a, co_c, slab_idx)
        nbrs_me=neighbor_table(a, me_o, slab_idx)
        def bin_(d):
            if d<2.1: return 'product-like'
            if d<3.0: return 'reactive-close'
            if d<4.0: return 'reactive-loose'
            if d<5.0: return 'separated'
            return 'thermodynamic'
        row={'sid':sid,'ads':ads,'idx':idx,'E_MLIP':round(E,4),'dE_from_global':round(dE,4),
             'site':'','region':'',
             'n_Pd_nbr':f'CO:{npd_co}/CH3O:{npd_me}',
             'n_O_nbr':f'CO:{no_co}/CH3O:{no_me}',
             'd_anchor_surf':round(d_min,3),
             'd_CO_bond':round(a.get_distance(co_c,co_o,mic=True),3),
             'd_CH3O_bond':round(a.get_distance(me_c,me_o,mic=True),3),
             'd_reactive':round(d_react,3),
             'd_bin':bin_(d_react),
             'site_combo':f'{s_co}+{s_me}',
             'neighbors_within_3.5A':
                 'C_CO:['+';'.join(f'{s}:{d:.2f}' for s,d in nbrs_co)+'] '+
                 'O_CH3:['+';'.join(f'{s}:{d:.2f}' for s,d in nbrs_me)+']'}
    # xy distance to nearest v4 pick (for the same anchor role)
    if ads=='coads':
        anchor=co_c
    v4=v4_pick_xy(sid, ads)
    if v4:
        xy_new=tuple(a.positions[anchor,:2])
        dists=[xy_mic(a, xy_new, xy_old) for _,_,_,xy_old in v4]
        row['xy_dist_to_nearest_v4']=round(min(dists),3)
        # which v4 idx it's closest to
        row['nearest_v4_idx']=v4[int(np.argmin(dists))][0]
    else:
        row['xy_dist_to_nearest_v4']=''; row['nearest_v4_idx']=''
    return row

# ---------- build refined additions ----------
def refine():
    v4v5=list(csv.DictReader(open(V5/'v4_vs_v5_comparison.csv')))
    new_keys=[(r['sid'],r['ads'],int(r['idx'])) for r in v4v5 if r['status'].startswith('v5-new')]
    rows=[]
    for k in new_keys:
        sid,ads,idx=k
        r=build_addition_row(sid, ads, idx)
        if k in MUST:
            r['priority']='MUST'
            r['selection_reason']='distance-bin gap (reactive) not covered by v4'
        elif k in REC:
            r['priority']='RECOMMENDED'
            if ads=='coads' and 'reactive' in r['d_bin']:
                r['selection_reason']=f'distinct-site {r["site_combo"]} + {r["d_bin"]}'
            else:
                r['selection_reason']=f'{r["d_bin"]} representative'
        elif ads in ('CO','CH3O'):
            r['priority']='OPTIONAL'
            r['selection_reason']=f'single-ads {r["site"]}/{r["region"]} diversity'
        else:
            # coads news not in MUST/REC list → REVIEW-NEEDED
            r['priority']='REVIEW-NEEDED'
            r['selection_reason']=(f'{r["d_bin"]} secondary pick / combo {r["site_combo"]} '
                                   f'— not a reactive rep')
        rows.append(r)

    # write refined proposed_additions_only.csv
    keys=['priority','sid','ads','idx','E_MLIP','dE_from_global','site','region',
          'n_Pd_nbr','n_O_nbr','d_anchor_surf','d_CO_bond','d_CH3O_bond',
          'd_reactive','d_bin','site_combo','xy_dist_to_nearest_v4','nearest_v4_idx',
          'selection_reason','neighbors_within_3.5A']
    order={'MUST':0,'RECOMMENDED':1,'REVIEW-NEEDED':2,'OPTIONAL':3}
    rows.sort(key=lambda r:(order[r['priority']], r['sid'], r['ads'], r['idx']))
    with open(V5/'proposed_additions_only.csv','w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    return rows

# ---------- pool_diagnostics.csv ----------
def write_pool_diag(rows):
    keys=['sid','ads','raw','unique_json','converged','valid']
    with open(V5/'pool_diagnostics.csv','w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)

# ---------- errata & report addendum ----------
def write_addendum(refined_rows, pool_rows):
    md=['# v5 addendum (2026-07-17 refinements)\n\n',
        '## 1. Errata\n\n',
        '- **v4 S1/CO composition** was mis-summarized as "5 br_PdPd" in the earlier text.\n',
        '  The v4 picks were **atop_Pd 1 (idx 64) + br_PdPd 4 (idx 3, 8, 28, 30)**\n',
        '  per `picks_analysis.csv`. The v5 selection log itself was correct; only the\n',
        '  narrative summary was wrong.\n\n',
        '## 2. Pool count clarification\n\n',
        '| sid | ads | raw (AutoAdsorbate) | unique_json (MLIP dedup) | converged | valid (final v5 pool) |\n',
        '|---|---|---|---|---|---|\n']
    for r in pool_rows:
        md.append(f"| {r['sid']} | {r['ads']} | {r['raw']} | {r['unique_json']} | {r['converged']} | {r['valid']} |\n")
    md.append('\n**Drop reasons**:\n\n')
    md.append('- `unique_json → converged`: MLIP relaxation did not reach EDIFFG stop → excluded.\n')
    md.append('- `converged → valid`: intramolecular geometry check (`valid_CO/CH3O/coads`) rejects\n')
    md.append('  broken adsorbates — C–O bond outside allowed range, CH₃ lost an H, or a 2nd O\n')
    md.append('  is within 2.0 Å of the C (CO₂-like collapse). This is the dominant loss on **S4/CO\n')
    md.append('  (52 → 26)** — half of the "converged" candidates collapsed to CO₂-like on the O-rich\n')
    md.append('  PdO₂(110) surface. `S2/CH3O 321→311`, `S3/CH3O 172→168` = a handful of dissociated methoxides.\n\n')
    md.append('## 3. Refined additions (proposed_additions_only.csv)\n\n')
    md.append('Priority sorted: MUST → RECOMMENDED → REVIEW-NEEDED → OPTIONAL.\n')
    md.append('Each row includes E_MLIP, ΔE(global), site+region, neighbor coordination,\n')
    md.append('distance bin (coads), xy MIC distance to nearest v4 pick, and neighbor distance table.\n\n')
    md.append('### MUST (4)\n')
    for r in refined_rows:
        if r['priority']!='MUST': continue
        md.append(f'- **{r["sid"]}/{r["ads"]} idx={r["idx"]}** — E={r["E_MLIP"]} eV, '
                  f'ΔE={r["dE_from_global"]} eV, d_reactive={r["d_reactive"]} Å '
                  f'({r["d_bin"]}), combo {r["site_combo"]}. '
                  f'xy dist to nearest v4 = {r["xy_dist_to_nearest_v4"]} Å.\n')
        md.append(f"    reason: {r['selection_reason']}\n")
    md.append('\n### RECOMMENDED (4)\n')
    for r in refined_rows:
        if r['priority']!='RECOMMENDED': continue
        md.append(f'- **{r["sid"]}/{r["ads"]} idx={r["idx"]}** — E={r["E_MLIP"]} eV, '
                  f'ΔE={r["dE_from_global"]} eV, d_reactive={r["d_reactive"]} Å '
                  f'({r["d_bin"]}), combo {r["site_combo"]}. '
                  f'xy dist to nearest v4 = {r["xy_dist_to_nearest_v4"]} Å.\n')
        md.append(f"    reason: {r['selection_reason']}\n")
    md.append('\n### REVIEW-NEEDED (auto-added coads picks NOT on researcher list)\n')
    for r in refined_rows:
        if r['priority']!='REVIEW-NEEDED': continue
        md.append(f'- {r["sid"]}/{r["ads"]} idx={r["idx"]} — {r["d_bin"]} '
                  f'combo {r["site_combo"]}, ΔE={r["dE_from_global"]} eV. '
                  f'*Not a reactive rep — likely redundant with global min or '
                  f'another separated/thermo pick. RECOMMENDED action: drop unless '
                  f'the site combo is genuinely novel.*\n')
    md.append('\n### OPTIONAL (single-ads v5-new — 다양성 이유이지만 ΔE 검토 필요)\n')
    for r in refined_rows:
        if r['priority']!='OPTIONAL': continue
        md.append(f'- {r["sid"]}/{r["ads"]} idx={r["idx"]} — site={r["site"]}/{r["region"]}, '
                  f'ΔE={r["dE_from_global"]} eV, d_anchor_surf={r["d_anchor_surf"]} Å. '
                  f'neighbors: {r["neighbors_within_3.5A"][:120]}...\n')
    md.append('\n## 4. Site verification (anchor-neighbor distance table)\n\n')
    md.append('Each MUST/RECOMMENDED candidate below shows all slab atoms within 3.5 Å of\n')
    md.append('the anchor, so the site label is grounded in actual coordination:\n\n')
    for r in refined_rows:
        if r['priority'] not in ('MUST','RECOMMENDED'): continue
        md.append(f'**{r["sid"]}/{r["ads"]} idx={r["idx"]}** (site_combo={r["site_combo"]})\n')
        md.append(f'  neighbors within 3.5 Å: {r["neighbors_within_3.5A"]}\n\n')
    md.append('## 5. S4 coadsorption\n\nStill intentionally excluded per researcher decision.\n\n')
    md.append('## 6. What was NOT done (guardrails)\n\n')
    md.append('- No file under `T1_16_DFT_L2/` or the pending65 bundle was touched.\n')
    md.append('- No new calc directories were created.\n')
    md.append('- No jobs submitted or cancelled.\n')
    md.append('- Existing S1/CO 5 completed candidates (idx 64, 3, 8, 28, 30) untouched;\n')
    md.append('  v4-only idx=30 is NOT replaced.\n')
    (V5/'addendum.md').write_text(''.join(md))

def main():
    print('Refining v5 additions per researcher spec...')
    refined = refine()
    pool_rows = pool_diagnostics()
    write_pool_diag(pool_rows)
    write_addendum(refined, pool_rows)
    # summary counts
    from collections import Counter
    c=Counter(r['priority'] for r in refined)
    print(f'\nPriority counts: {dict(c)}')
    for r in refined:
        if r['priority']=='MUST':
            print(f"  MUST         {r['sid']}/{r['ads']:<5} idx={r['idx']:>5}  d_reactive={r['d_reactive']}  ΔE={r['dE_from_global']}")
    for r in refined:
        if r['priority']=='RECOMMENDED':
            print(f"  RECOMMENDED  {r['sid']}/{r['ads']:<5} idx={r['idx']:>5}  d_reactive={r['d_reactive']}  ΔE={r['dE_from_global']}")

if __name__=='__main__': main()
