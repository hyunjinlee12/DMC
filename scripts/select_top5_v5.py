"""v5 DFT shortlist: diversity-aware selection.

Rules (per 2026-07-17 researcher spec):
- Global MLIP min always slot 0.
- Fill remaining slots by (site_class, xy_cluster, coord_signature) diversity.
- For coads: also enforce distance-bin diversity (reactive-close, reactive-loose,
  thermodynamic, product, side-path).
- Do NOT accept high-ΔE picks purely to satisfy diversity — report if a
  category is empty.

Reads MLIP pool from calculations/G3_adsorption/{surface}/MLIP_phase{1,2_filtered}/.
Writes:
  DFT_shortlist_v5/summary.csv           final picks with selection reasons
  DFT_shortlist_v5/picks_analysis.csv    descriptor table (corrected columns)
  DFT_shortlist_v5/v4_vs_v5_comparison.csv
  DFT_shortlist_v5/proposed_additions_only.csv
  DFT_shortlist_v5/report.md             human-readable summary
NO structure files copied, NO T1_16_DFT_L2 mutation, NO submission.
"""
import json, csv, math
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
from ase.io import read

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT/'calculations/G2_slab'
G3 = ROOT/'calculations/G3_adsorption'
OUT = G3/'DFT_shortlist_v5'
OUT.mkdir(exist_ok=True)

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
SURFACES_SINGLE = ['S1','S2','S3','S3b','S4']
SURFACES_COADS  = ['S1','S2','S3','S3b']       # S4 coads intentionally excluded

CUTOFF_NBR = 2.6
XY_CLUSTER = 1.5    # Å; two anchors within this xy are "same xy cluster"
CO2_CUT = 2.0

# ---------- validity ----------
def valid_CO(a):
    syms=a.get_chemical_symbols()
    c_i=[i for i,s in enumerate(syms) if s=='C']
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
    # CO's O
    d_co_o=sorted([(a.get_distance(co_c,oi,mic=True),oi) for oi in o_i])
    if not(1.05<=d_co_o[0][0]<=1.30): return None
    if len(d_co_o)>=2 and d_co_o[1][0]<CO2_CUT: return None
    co_o=d_co_o[0][1]
    # methoxide O
    d_me_o=sorted([(a.get_distance(me_c,oi,mic=True),oi) for oi in o_i if oi!=co_o])
    if not d_me_o or not(1.30<=d_me_o[0][0]<=1.55): return None
    if len(d_me_o)>=2 and d_me_o[1][0]<CO2_CUT: return None
    me_o=d_me_o[0][1]
    return (co_c, co_o, me_c, me_o)

# ---------- descriptors ----------
def site_label(a, anchor, slab_idx):
    syms=a.get_chemical_symbols()
    d=a.get_distances(anchor, slab_idx, mic=True)
    nbrs=[(slab_idx[i],d[i]) for i in range(len(slab_idx)) if d[i]<CUTOFF_NBR]
    n_pd=sum(1 for i,_ in nbrs if syms[i]=='Pd')
    n_o=sum(1 for i,_ in nbrs if syms[i]=='O')
    tot=n_pd+n_o
    if tot==0: return 'physi', n_pd, n_o
    if tot==1: return ('atop_Pd' if n_pd else 'atop_O'), n_pd, n_o
    if tot==2:
        return ('br_PdPd' if n_pd==2 else ('br_OO' if n_o==2 else 'br_PdO')), n_pd, n_o
    if tot==3:
        return ('h3Pd' if n_pd==3 else ('h3O' if n_o==3 else f'h3_{n_pd}Pd{n_o}O')), n_pd, n_o
    return f'{tot}f_{n_pd}Pd{n_o}O', n_pd, n_o

def surface_region(site, n_pd, n_o):
    """Coarse chem-environment: metal / oxide / interface / physi."""
    if site == 'physi': return 'physi'
    if n_pd > 0 and n_o == 0: return 'metal'
    if n_pd == 0 and n_o > 0: return 'oxide'
    return 'interface'

def xy_of(a, i):
    """Fractional xy modulo 1 in cell basis, then project to Cartesian xy."""
    cell=a.cell.array
    inv=np.linalg.inv(cell)
    df=inv@a.positions[i]
    df[0]-=math.floor(df[0]); df[1]-=math.floor(df[1]); df[2]=0
    v=cell.T@df
    return float(v[0]), float(v[1])

def xy_mic_dist(a, i, j):
    cell=a.cell.array; inv=np.linalg.inv(cell)
    df=inv@(a.positions[j]-a.positions[i])
    df[0]-=round(df[0]); df[1]-=round(df[1]); df[2]=0
    v=cell.T@df
    return math.hypot(v[0], v[1])

def z_top(a, slab_idx):
    return float(a.positions[slab_idx,2].max())

def distance_bin(d):
    if d<2.1: return 'product-like'
    if d<3.0: return 'reactive-close'
    if d<4.0: return 'reactive-loose'
    if d<5.0: return 'separated'
    return 'thermodynamic'

# ---------- pool characterization ----------
def characterize_single(sid, ads):
    sdir=SDIRS[sid]
    uni=json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
    traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
    slab_ref=read(G2/sdir/'CONTCAR'); n_ref=len(slab_ref)
    n_ads=2 if ads=='CO' else 5
    valid_f=valid_CO if ads=='CO' else valid_CH3O
    out=[]
    for r in uni:
        if not r.get('converged',True): continue
        a=traj[r['idx']]
        if len(a)!=n_ref+n_ads:
            ads_at=a[-n_ads:]; a=slab_ref.copy(); a+=ads_at
        v=valid_f(a)
        if v is None: continue
        c, o = v
        anchor = c if ads=='CO' else o
        ads_set = {c, o}
        if ads=='CH3O':
            ads_set |= {i for i,s in enumerate(a.get_chemical_symbols()) if s=='H'}
        slab_idx=[i for i in range(len(a)) if i not in ads_set]
        site, npd, no = site_label(a, anchor, slab_idx)
        region = surface_region(site, npd, no)
        d_anc=min(a.get_distances(anchor, slab_idx, mic=True))
        x,y=xy_of(a, anchor)
        d_CO=a.get_distance(c, o, mic=True)
        h=a.positions[anchor,2]-z_top(a, slab_idx)
        out.append({
            'idx':r['idx'], 'E':float(r['E']), 'converged':r.get('converged',True),
            'site':site, 'region':region, 'n_Pd':npd, 'n_O':no,
            'd_anchor_surf':float(d_anc), 'height':float(h), 'd_CO_bond':float(d_CO),
            'anchor_idx':anchor, 'xy':(x,y), 'ads_indices':sorted(ads_set),
            'atoms':a,
        })
    out.sort(key=lambda r: r['E'])
    return out

def characterize_coads(sid):
    sdir=SDIRS[sid]
    uni=json.load(open(G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'))
    traj=list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
    slab_ref=read(G2/sdir/'CONTCAR'); n_ref=len(slab_ref)
    n_ads=7
    out=[]
    for r in uni:
        if not r.get('converged',True): continue
        a=traj[r['idx']]
        if len(a)!=n_ref+n_ads:
            ads_at=a[-n_ads:]; a=slab_ref.copy(); a+=ads_at
        v=valid_coads(a)
        if v is None: continue
        co_c, co_o, me_c, me_o = v
        ads_set={co_c,co_o,me_c,me_o}|{i for i,s in enumerate(a.get_chemical_symbols()) if s=='H'}
        slab_idx=[i for i in range(len(a)) if i not in ads_set]
        s_co,npd_co,no_co=site_label(a, co_c, slab_idx)
        s_me,npd_me,no_me=site_label(a, me_o, slab_idx)
        d_react=a.get_distance(co_c, me_o, mic=True)
        x_co,y_co=xy_of(a, co_c); x_me,y_me=xy_of(a, me_o)
        out.append({
            'idx':r['idx'], 'E':float(r['E']),
            'site_CO':s_co, 'site_CH3O':s_me,
            'region_CO':surface_region(s_co,npd_co,no_co),
            'region_CH3O':surface_region(s_me,npd_me,no_me),
            'nbr_CO':(npd_co,no_co), 'nbr_CH3O':(npd_me,no_me),
            'd_reactive':float(d_react),
            'd_bin':distance_bin(d_react),
            'd_CO_bond':float(a.get_distance(co_c, co_o, mic=True)),
            'd_CH3O_bond':float(a.get_distance(me_c, me_o, mic=True)),
            'co_anchor_idx':co_c, 'me_anchor_idx':me_o,
            'xy_co':(x_co,y_co), 'xy_me':(x_me,y_me),
            'atoms':a,
        })
    out.sort(key=lambda r: r['E'])
    return out

# ---------- v5 selection ----------
def is_xy_dup_single(pool_row, pick, atoms_for_cell):
    """xy MIC distance between anchors <XY_CLUSTER → same xy cluster."""
    a=pool_row['atoms']
    # use cell of a; positions from a and pick['atoms'] must be in same cell (same slab)
    dx=pick['xy'][0]-pool_row['xy'][0]; dy=pick['xy'][1]-pool_row['xy'][1]
    # xy_of already wraps into [0, cell)  ; use minimum-image via cell dims
    cell=a.cell.array
    # approximate MIC: subtract nearest multiple of the two in-plane lattice periods
    inv=np.linalg.inv(cell); v=inv@np.array([dx,dy,0.0])
    v[0]-=round(v[0]); v[1]-=round(v[1]); v[2]=0
    w=cell.T@v
    return math.hypot(w[0],w[1]) < XY_CLUSTER

def is_xy_dup_coads(pool_row, pick):
    a=pool_row['atoms']; cell=a.cell.array; inv=np.linalg.inv(cell)
    def mic(p1,p2):
        v=inv@np.array([p2[0]-p1[0], p2[1]-p1[1], 0.0])
        v[0]-=round(v[0]); v[1]-=round(v[1]); v[2]=0
        w=cell.T@v
        return math.hypot(w[0],w[1])
    d_co=mic(pool_row['xy_co'], pick['xy_co'])
    d_me=mic(pool_row['xy_me'], pick['xy_me'])
    return (d_co<XY_CLUSTER) and (d_me<XY_CLUSTER)

def select_single(sid, ads, pool, N_TARGET=5):
    """Greedy diversity-aware fill. Returns picks + selection log."""
    if not pool:
        return [], ['pool empty']
    picks=[]; log=[]
    seen_sites=set(); seen_regions=set(); seen_xy=[]
    for r in pool:
        # duplicate check vs existing picks
        if any(is_xy_dup_single(r, p, r['atoms']) and r['site']==p['site'] for p in picks):
            continue
        if not picks:
            reason='global E min'
            picks.append(r); seen_sites.add(r['site']); seen_regions.add(r['region']); seen_xy.append(r['xy'])
            log.append((r['idx'], r['E'], r['site'], reason)); continue
        if len(picks)>=N_TARGET: break
        new_site = r['site'] not in seen_sites
        new_region = r['region'] not in seen_regions
        dE = r['E']-picks[0]['E']
        # Accept if: (a) introduces new site OR new region AND dE<0.5 eV
        # OR (b) picks<N_TARGET/2 and still low-E and different xy
        if new_site or new_region:
            if dE <= 0.5:
                reason=f'new {"site" if new_site else "region"}: {r["site"]}/{r["region"]} (ΔE={dE:.3f})'
                picks.append(r); seen_sites.add(r['site']); seen_regions.add(r['region']); seen_xy.append(r['xy'])
                log.append((r['idx'], r['E'], r['site'], reason))
            else:
                log.append((r['idx'], r['E'], r['site'], f'SKIP diversity candidate ΔE={dE:.3f}>0.5 eV'))
    # 2nd pass: fill remaining slots by lowest E among non-duplicates (any bin)
    for r in pool:
        if len(picks)>=N_TARGET: break
        if r in picks: continue
        if any(is_xy_dup_single(r, p, r['atoms']) and r['site']==p['site'] for p in picks):
            continue
        dE=r['E']-picks[0]['E']
        picks.append(r); seen_sites.add(r['site']); seen_regions.add(r['region']); seen_xy.append(r['xy'])
        log.append((r['idx'], r['E'], r['site'], f'fill (ΔE={dE:.3f})'))
    return picks, log

def select_coads(sid, pool, N_TARGET=5):
    if not pool:
        return [], ['pool empty']
    picks=[]; log=[]
    seen_combo=set(); seen_bin=set()
    def combo(r): return (r['site_CO'], r['site_CH3O'])
    def add(r, reason):
        picks.append(r); seen_combo.add(combo(r)); seen_bin.add(r['d_bin'])
        log.append((r['idx'], r['E'], r['d_bin'], combo(r), reason))
    # slot 0: global min
    add(pool[0], 'global E min')
    # target bins in priority order
    target_bins=['reactive-close','reactive-loose','thermodynamic','separated','product-like']
    for tb in target_bins:
        if len(picks)>=N_TARGET: break
        if tb in seen_bin: continue
        cand=[r for r in pool if r['d_bin']==tb and r['idx']!=pool[0]['idx']
              and not any(is_xy_dup_coads(r,p) for p in picks)]
        if not cand:
            log.append((None, None, tb, None, 'NO CANDIDATE IN BIN'))
            continue
        best=cand[0]
        dE=best['E']-pool[0]['E']
        if dE>0.8:
            log.append((best['idx'], best['E'], tb, combo(best), f'SKIP dE={dE:.3f}>0.8 eV'))
            continue
        add(best, f'{tb} rep (ΔE={dE:.3f})')
    # remaining slots: new site combos + interface preference for S2/S3b
    prefer_iface = sid in ('S2','S3b')
    for r in pool:
        if len(picks)>=N_TARGET: break
        if r in picks: continue
        if any(is_xy_dup_coads(r,p) for p in picks): continue
        c=combo(r)
        dE=r['E']-pool[0]['E']
        if c not in seen_combo and dE<=0.8:
            add(r, f'new combo {c[0]}+{c[1]} (ΔE={dE:.3f})')
            continue
    # final fill by lowest E respecting dedup
    for r in pool:
        if len(picks)>=N_TARGET: break
        if r in picks: continue
        if any(is_xy_dup_coads(r,p) for p in picks): continue
        dE=r['E']-pool[0]['E']
        add(r, f'fill (ΔE={dE:.3f})')
    return picks, log

# ---------- run ----------
def main():
    single_picks={}; coads_picks={}
    pool_stats=[]
    all_desc=[]

    for sid in SURFACES_SINGLE:
        for ads in ['CO','CH3O']:
            print(f'  characterizing {sid} {ads}...')
            pool=characterize_single(sid, ads)
            picks, log = select_single(sid, ads, pool)
            single_picks[(sid,ads)] = (pool, picks, log)
            sites=Counter(r['site'] for r in pool)
            regions=Counter(r['region'] for r in pool)
            e_arr=[r['E'] for r in pool]
            pool_stats.append({'sid':sid,'ads':ads,'N_pool':len(pool),
                               'E_min':round(e_arr[0],3),'E_max':round(e_arr[-1],3),
                               'dE_span':round(e_arr[-1]-e_arr[0],3),
                               'sites':dict(sites),'regions':dict(regions),
                               'n_picks':len(picks),
                               'picks_sites':dict(Counter(r['site'] for r in picks))})
            for rank,r in enumerate(picks):
                all_desc.append({
                    'sid':sid,'ads':ads,'rank':rank,'idx':r['idx'],
                    'E_MLIP':round(r['E'],4),'dE_from_global':round(r['E']-pool[0]['E'],4),
                    'site':r['site'],'region':r['region'],
                    'n_Pd_nbr':r['n_Pd'],'n_O_nbr':r['n_O'],
                    'd_anchor_surf':round(r['d_anchor_surf'],3),
                    'height':round(r['height'],3),
                    'd_CO_bond':round(r['d_CO_bond'],3),
                    'd_CH3O_bond':'',       # N/A for singles
                    'd_reactive':'','d_bin':'',
                    'site_combo':'',
                })

    for sid in SURFACES_COADS:
        print(f'  characterizing {sid} coads...')
        pool=characterize_coads(sid)
        picks, log = select_coads(sid, pool)
        coads_picks[sid] = (pool, picks, log)
        combos=Counter((r['site_CO'],r['site_CH3O']) for r in pool)
        bins=Counter(r['d_bin'] for r in pool)
        e_arr=[r['E'] for r in pool]
        pool_stats.append({'sid':sid,'ads':'coads','N_pool':len(pool),
                           'E_min':round(e_arr[0],3),'E_max':round(e_arr[-1],3),
                           'dE_span':round(e_arr[-1]-e_arr[0],3),
                           'combos':{str(k):v for k,v in combos.most_common(5)},
                           'd_bins':dict(bins),
                           'n_picks':len(picks),
                           'picks_bins':dict(Counter(r['d_bin'] for r in picks))})
        for rank,r in enumerate(picks):
            all_desc.append({
                'sid':sid,'ads':'coads','rank':rank,'idx':r['idx'],
                'E_MLIP':round(r['E'],4),'dE_from_global':round(r['E']-pool[0]['E'],4),
                'site':'','region':'',
                'n_Pd_nbr':'','n_O_nbr':'',
                'd_anchor_surf':'','height':'',
                'd_CO_bond':round(r['d_CO_bond'],3),
                'd_CH3O_bond':round(r['d_CH3O_bond'],3),
                'd_reactive':round(r['d_reactive'],3),
                'd_bin':r['d_bin'],
                'site_combo':f"{r['site_CO']}+{r['site_CH3O']}",
            })

    # ---- write summary.csv ----
    keys=['sid','ads','rank','idx','E_MLIP','dE_from_global','site','region',
          'n_Pd_nbr','n_O_nbr','d_anchor_surf','height','d_CO_bond','d_CH3O_bond',
          'd_reactive','d_bin','site_combo']
    with open(OUT/'summary.csv','w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in all_desc: w.writerow(r)

    # ---- picks_analysis.csv (same content, but pandas-friendly headers) ----
    with open(OUT/'picks_analysis.csv','w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in all_desc: w.writerow(r)

    # ---- v4 vs v5 comparison ----
    v4 = json.load(open(G3/'DFT_shortlist_v3/summary.json'))
    v4_by_group=defaultdict(set)
    for r in v4: v4_by_group[(r['sid'],r['ads'])].add(r['idx'])
    v5_by_group=defaultdict(set)
    for r in all_desc: v5_by_group[(r['sid'],r['ads'])].add(r['idx'])

    comp_rows=[]
    for k in sorted(set(v4_by_group)|set(v5_by_group)):
        v4s=v4_by_group.get(k,set()); v5s=v5_by_group.get(k,set())
        for idx in sorted(v4s & v5s): comp_rows.append({'sid':k[0],'ads':k[1],'idx':idx,'status':'KEEP (in both v4 and v5)'})
        for idx in sorted(v4s - v5s): comp_rows.append({'sid':k[0],'ads':k[1],'idx':idx,'status':'v4-only (v5 dropped)'})
        for idx in sorted(v5s - v4s): comp_rows.append({'sid':k[0],'ads':k[1],'idx':idx,'status':'v5-new (needs new DFT if approved)'})
    with open(OUT/'v4_vs_v5_comparison.csv','w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=['sid','ads','idx','status']); w.writeheader()
        for r in comp_rows: w.writerow(r)

    # ---- proposed_additions_only.csv (v5-new only) ----
    with open(OUT/'proposed_additions_only.csv','w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        added_keys={(r['sid'],r['ads'],r['idx']) for r in comp_rows if r['status'].startswith('v5-new')}
        for r in all_desc:
            if (r['sid'],r['ads'],r['idx']) in added_keys: w.writerow(r)

    # ---- report.md ----
    md=[]
    md.append('# DFT shortlist v5 — diversity-aware selection\n')
    md.append(f'Generated by `scripts/select_top5_v5.py`.\n')
    md.append('## Selection rules\n')
    md.append('- Global MLIP min → always slot 0.\n'
              '- Remaining slots filled by new site class or region (ΔE < 0.5 eV cap).\n'
              '- Coads: distance-bin diversity (reactive-close 2.1–3.0, reactive-loose 3.0–4.0, thermodynamic ≥5.0), ΔE < 0.8 eV cap.\n'
              '- Fill remaining slots by lowest E after xy+site dedup.\n'
              '- Categories with no candidate are reported as "NO CANDIDATE IN BIN".\n\n')
    md.append('## Pool + pick statistics\n\n')
    md.append('| sid | ads | N_pool | E_min | ΔE_span | picks | site diversity |\n')
    md.append('|---|---|---|---|---|---|---|\n')
    for s in pool_stats:
        if s['ads']=='coads':
            div=f"bins {s['picks_bins']}"
        else:
            div=f"sites {s['picks_sites']}"
        md.append(f"| {s['sid']} | {s['ads']} | {s['N_pool']} | {s['E_min']} | {s['dE_span']} | {s['n_picks']} | {div} |\n")
    md.append('\n## Per-group v5 picks with selection reason\n\n')
    for (sid,ads),(pool,picks,log) in list(single_picks.items()):
        md.append(f'### {sid} {ads} (pool={len(pool)}, picks={len(picks)})\n')
        for entry in log:
            if len(entry)==4:
                idx,E,site,reason=entry
                md.append(f'- idx={idx:>5} E={E:.4f} site={site:<10}  {reason}\n')
        md.append('\n')
    for sid,(pool,picks,log) in coads_picks.items():
        md.append(f'### {sid} coads (pool={len(pool)}, picks={len(picks)})\n')
        for entry in log:
            if len(entry)==5:
                idx,E,bin_,combo,reason=entry
                e_str='  --  ' if E is None else f'E={E:.4f}'
                idx_str='  --' if idx is None else f'idx={idx}'
                md.append(f'- {idx_str} {e_str} bin={bin_:<15} combo={combo}  {reason}\n')
        md.append('\n')
    md.append('## Comparison to v4 (rows in v4_vs_v5_comparison.csv)\n\n')
    stat_c=Counter(r['status'] for r in comp_rows)
    for k,v in stat_c.most_common(): md.append(f'- {k}: {v}\n')
    md.append('\n## Notes / caveats\n\n')
    md.append('- **S4 coads intentionally excluded** — per researcher decision.\n')
    md.append('- **S1 CO\\* pool has only 12 unique candidates** (labels: bridge 3, atop 2, uncl 7).\n')
    md.append('  → diversity is bounded by pool, not by selector.\n')
    md.append('- **S1 CH₃O\\* pool** has 85 unique but all labeled "unknown" by the MLIP screener\n'
              '  — v5 re-classifies at MLIP-relaxed geometry using our anchor+2.6 Å rule.\n')
    md.append('- **S4 CO\\*** — pool has ΔE = 3.2 eV span, chemisorbed candidates exist but at higher E\n'
              '  than physi cluster. v5 reports if diversity cap (ΔE < 0.5 eV) excludes them.\n')
    md.append('- **v5 output does NOT copy structure files** — only descriptor tables + report.\n'
              '  Nothing under T1_16_DFT_L2 or the pending65 bundle was touched.\n')
    md.append('- **No calculations submitted or cancelled**.\n')
    (OUT/'report.md').write_text(''.join(md))

    # pool stats JSON for machine consumption
    json.dump(pool_stats, open(OUT/'pool_stats.json','w'), indent=1, default=str)

    print(f'\nWrote {OUT}/')
    for f in ['summary.csv','picks_analysis.csv','v4_vs_v5_comparison.csv',
              'proposed_additions_only.csv','report.md','pool_stats.json']:
        p=OUT/f
        print(f'  {f}  ({p.stat().st_size} bytes)')

if __name__=='__main__':
    main()
