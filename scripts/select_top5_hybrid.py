"""Hybrid top-5 selection:
Pass 1: site-diverse (distinct site combos, ΔE > 0.02 eV) — as before
Pass 2: if <5 after Pass 1, fill remaining with duplicate-filter (ΔE<0.03 + xy<1.5 = duplicate)
"""
import json, shutil
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT/'calculations/G3_adsorption'
G2 = ROOT/'calculations/G2_slab'
OUT = G3/'DFT_shortlist_v3'
if OUT.exists(): shutil.rmtree(OUT)
OUT.mkdir()

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
SURFACES_SINGLE = ['S1','S2','S3','S3b','S4']
SURFACES_COADS  = ['S1','S2','S3b','S3']
N_TARGET = 5
DELTA_E_PASS1 = 0.02   # ΔE >= this to be considered new site rank
DELTA_E_DUP = 0.03     # duplicate detection (Pass 2)
XY_DUP = 1.5           # Å
CO2_CUTOFF = 2.0
E_CEILING = 0.5        # eV — top-1 대비 이내만 허용 (Pass 2 안 뽑음 넘으면)


def site_type(a, anchor, n_ads):
    syms=a.get_chemical_symbols()
    sub=list(range(len(a)-n_ads))
    d=a.get_distances(anchor,sub,mic=True)
    nbrs=[(sub[i],d[i]) for i in range(len(sub)) if d[i]<2.6]
    n_pd=sum(1 for i,_ in nbrs if syms[i]=='Pd')
    n_o=sum(1 for i,_ in nbrs if syms[i]=='O')
    tot=n_pd+n_o
    if tot==0: return 'physi'
    if tot==1: return 'atop_Pd' if n_pd else 'atop_O'
    if tot==2: return 'br_PdPd' if n_pd==2 else ('br_OO' if n_o==2 else 'br_PdO')
    if tot==3: return 'h3Pd' if n_pd==3 else ('h3O' if n_o==3 else f'h3({n_pd}Pd{n_o}O)')
    return f'{tot}f({n_pd}Pd{n_o}O)'


def valid_CO(a):
    syms=a.get_chemical_symbols()
    c_idx=[i for i,s in enumerate(syms) if s=='C']
    o_idx=[i for i,s in enumerate(syms) if s=='O']
    if len(c_idx)!=1: return None
    c=c_idx[0]
    d_o=sorted([(a.get_distance(c,oi,mic=True), oi) for oi in o_idx])
    if not (1.05<=d_o[0][0]<=1.30): return None
    if len(d_o)>=2 and d_o[1][0]<CO2_CUTOFF: return None
    return c

def valid_CH3O(a):
    syms=a.get_chemical_symbols()
    c_idx=[i for i,s in enumerate(syms) if s=='C']
    h_idx=[i for i,s in enumerate(syms) if s=='H']
    o_idx=[i for i,s in enumerate(syms) if s=='O']
    if len(c_idx)!=1 or len(h_idx)!=3: return None
    c=c_idx[0]
    if not all(0.90<=a.get_distance(c,h,mic=True)<=1.25 for h in h_idx): return None
    d_o=sorted([(a.get_distance(c,oi,mic=True), oi) for oi in o_idx])
    if not (1.30<=d_o[0][0]<=1.55): return None
    if len(d_o)>=2 and d_o[1][0]<CO2_CUTOFF: return None
    return d_o[0][1]

def valid_coads(a):
    syms=a.get_chemical_symbols()
    c_idx=[i for i,s in enumerate(syms) if s=='C']
    h_idx=[i for i,s in enumerate(syms) if s=='H']
    o_idx=[i for i,s in enumerate(syms) if s=='O']
    if len(c_idx)!=2 or len(h_idx)!=3: return None
    me_c, co_c = None, None
    for c in c_idx:
        n_h = sum(1 for h in h_idx if a.get_distance(c,h,mic=True)<1.3)
        if n_h==3: me_c=c
        else: co_c=c
    if me_c is None or co_c is None: return None
    d_o_co=sorted([(a.get_distance(co_c,oi,mic=True), oi) for oi in o_idx])
    if not (1.05<=d_o_co[0][0]<=1.30): return None
    if len(d_o_co)>=2 and d_o_co[1][0]<CO2_CUTOFF: return None
    if not all(0.90<=a.get_distance(me_c,h,mic=True)<=1.25 for h in h_idx): return None
    d_o_me=sorted([(a.get_distance(me_c,oi,mic=True), oi) for oi in o_idx])
    if not (1.30<=d_o_me[0][0]<=1.55): return None
    if len(d_o_me)>=2 and d_o_me[1][0]<CO2_CUTOFF: return None
    return (co_c, d_o_me[0][1])


def xy_mic(a1, i, a2, j):
    cell = a1.cell.array
    dx = a1.positions[i] - a2.positions[j]
    inv = np.linalg.inv(cell)
    df = inv @ dx
    df[0] -= np.round(df[0]); df[1] -= np.round(df[1]); df[2]=0
    v = cell.T @ df
    return float(np.sqrt(v[0]**2+v[1]**2))


def get_valid_candidates(sid, ads_key):
    sdir = SDIRS[sid]
    slab = read(G2/sdir/'CONTCAR'); n_sub = len(slab)
    if ads_key == 'coads':
        n_ads = 7
        unique = json.load(open(G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'))
        traj = list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
        valid_f = valid_coads
    else:
        n_ads = 2 if ads_key == 'CO' else 5
        unique = json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads_key}.json'))
        traj = list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads_key}.traj', index=':'))
        valid_f = valid_CO if ads_key == 'CO' else valid_CH3O
    out=[]
    for r in sorted(unique, key=lambda r: r['E']):
        if not r.get('converged', True): continue
        a = traj[r['idx']]
        if len(a) != n_sub + n_ads:
            ads_at = a[-n_ads:]; a = slab.copy(); a += ads_at
        v = valid_f(a)
        if v is None: continue
        # site combo
        if ads_key == 'coads':
            combo = (site_type(a, v[0], n_ads), site_type(a, v[1], n_ads))
        else:
            combo = (site_type(a, v, n_ads),)
        out.append({'idx':r['idx'], 'E':r['E'], 'atoms':a, 'anc':v, 'combo':combo})
    return out


def pick_hybrid(sid, ads_key):
    cands = get_valid_candidates(sid, ads_key)
    if not cands: return []
    E_top1 = cands[0]['E']
    # Median 계산 (valid candidate 만 대상)
    Es = [c['E'] for c in cands]
    E_median = np.median(Es)
    picks=[]
    combos_seen=set()
    # Pass 1: site-diverse — 단 E < median 인 후보만
    for c in cands:
        if c['E'] > E_median: break
        if c['combo'] in combos_seen: continue
        if picks and abs(c['E'] - picks[-1]['E']) < DELTA_E_PASS1:
            continue
        combos_seen.add(c['combo'])
        picks.append(c)
        if len(picks) >= N_TARGET: return picks
    # Pass 2: 나머지 채움 — E < median, top-1+ceiling 이내, duplicate 아니면
    for c in cands:
        if c['E'] > E_median: break
        if c['E'] > E_top1 + E_CEILING: break
        if any(c['idx']==p['idx'] for p in picks): continue
        dup=False
        for p in picks:
            if abs(c['E'] - p['E']) < DELTA_E_DUP:
                if ads_key == 'coads':
                    d1 = xy_mic(c['atoms'], c['anc'][0], p['atoms'], p['anc'][0])
                    d2 = xy_mic(c['atoms'], c['anc'][1], p['atoms'], p['anc'][1])
                    if d1 < XY_DUP and d2 < XY_DUP: dup=True; break
                else:
                    d = xy_mic(c['atoms'], c['anc'], p['atoms'], p['anc'])
                    if d < XY_DUP: dup=True; break
        if dup: continue
        picks.append(c)
        if len(picks) >= N_TARGET: return picks
    return picks


def fix_bottom_half(atoms, n_sub):
    z = atoms.positions[:n_sub, 2]; zm = np.median(z)
    fixed = [i for i in range(n_sub) if atoms.positions[i,2] < zm]
    atoms.set_constraint(FixAtoms(indices=fixed))


summary=[]
for sid in SURFACES_SINGLE:
    slab=read(G2/SDIRS[sid]/'CONTCAR'); n_sub=len(slab)
    for ads in ['CO','CH3O']:
        picks = pick_hybrid(sid, ads)
        out_dir = OUT/sid/ads; out_dir.mkdir(parents=True)
        for k, p in enumerate(picks):
            atoms = p['atoms'].copy()
            fix_bottom_half(atoms, n_sub)
            combo_str = '_'.join(p['combo'])
            fname = f"{k:02d}_{ads}_{combo_str}_idx{p['idx']:05d}.vasp"
            write(str(out_dir/fname), atoms, format='vasp', direct=True, sort=True, vasp5=True)
            summary.append({'sid':sid,'ads':ads,'rank':k,'idx':p['idx'],'E':p['E'],'combo':'_'.join(p['combo'])})
        print(f'  {sid} {ads:<6}: {len(picks)} picks')

for sid in SURFACES_COADS:
    slab=read(G2/SDIRS[sid]/'CONTCAR'); n_sub=len(slab)
    picks = pick_hybrid(sid, 'coads')
    out_dir = OUT/sid/'coads'; out_dir.mkdir(parents=True)
    for k, p in enumerate(picks):
        atoms = p['atoms'].copy()
        fix_bottom_half(atoms, n_sub)
        combo_str = f"CO-{p['combo'][0]}_OMe-{p['combo'][1]}"
        fname = f"{k:02d}_coads_{combo_str}_idx{p['idx']:05d}.vasp"[:80]
        write(str(out_dir/fname), atoms, format='vasp', direct=True, sort=True, vasp5=True)
        summary.append({'sid':sid,'ads':'coads','rank':k,'idx':p['idx'],'E':p['E'],'combo':combo_str})
    print(f'  {sid} coads : {len(picks)} picks')

json.dump(summary, open(OUT/'summary.json','w'), indent=2)
print(f'\nTotal picks: {len(summary)}')
