"""Top-5 per surface × ads. Exclusion = E 겹침 + xy 위치 겹침 (site label 사용 X).

기준:
  - 확실한 중복: ΔE < 0.03 eV AND anchor xy MIC 거리 < 1.5 Å
  - 위 조건 아니면 다른 구조로 인정 → pick
  - Filter: converged + intramol valid + no CO2 formation
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

DELTA_E_DUP = 0.03    # eV; below → same E cluster
XY_DUP = 1.5          # Å; below → same physical position
N_TARGET = 5
CO2_CUTOFF = 2.0


def valid_CO(atoms):
    syms = atoms.get_chemical_symbols()
    c_idx = [i for i,s in enumerate(syms) if s=='C']
    o_idx = [i for i,s in enumerate(syms) if s=='O']
    if len(c_idx) != 1: return None
    c = c_idx[0]
    d_o = sorted([(atoms.get_distance(c,oi,mic=True), oi) for oi in o_idx])
    if not (1.05 <= d_o[0][0] <= 1.30): return None
    if len(d_o) >= 2 and d_o[1][0] < CO2_CUTOFF: return None
    return c


def valid_CH3O(atoms):
    syms = atoms.get_chemical_symbols()
    c_idx = [i for i,s in enumerate(syms) if s=='C']
    h_idx = [i for i,s in enumerate(syms) if s=='H']
    o_idx = [i for i,s in enumerate(syms) if s=='O']
    if len(c_idx) != 1 or len(h_idx) != 3: return None
    c = c_idx[0]
    d_h = [atoms.get_distance(c,h,mic=True) for h in h_idx]
    if not all(0.90<=d<=1.25 for d in d_h): return None
    d_o = sorted([(atoms.get_distance(c,oi,mic=True), oi) for oi in o_idx])
    if not (1.30 <= d_o[0][0] <= 1.55): return None
    if len(d_o) >= 2 and d_o[1][0] < CO2_CUTOFF: return None
    return d_o[0][1]


def valid_coads(atoms):
    syms = atoms.get_chemical_symbols()
    c_idx = [i for i,s in enumerate(syms) if s=='C']
    h_idx = [i for i,s in enumerate(syms) if s=='H']
    o_idx = [i for i,s in enumerate(syms) if s=='O']
    if len(c_idx) != 2 or len(h_idx) != 3: return None
    me_c, co_c = None, None
    for c in c_idx:
        n_h = sum(1 for h in h_idx if atoms.get_distance(c,h,mic=True) < 1.3)
        if n_h == 3: me_c = c
        else: co_c = c
    if me_c is None or co_c is None: return None
    d_o_co = sorted([(atoms.get_distance(co_c,oi,mic=True), oi) for oi in o_idx])
    if not (1.05 <= d_o_co[0][0] <= 1.30): return None
    if len(d_o_co) >= 2 and d_o_co[1][0] < CO2_CUTOFF: return None
    d_ch = [atoms.get_distance(me_c,h,mic=True) for h in h_idx]
    if not all(0.90<=d<=1.25 for d in d_ch): return None
    d_o_me = sorted([(atoms.get_distance(me_c,oi,mic=True), oi) for oi in o_idx])
    if not (1.30 <= d_o_me[0][0] <= 1.55): return None
    if len(d_o_me) >= 2 and d_o_me[1][0] < CO2_CUTOFF: return None
    return (co_c, d_o_me[0][1])


def mic_xy_dist(a, i, j):
    """MIC distance in xy only."""
    cell = a.cell.array
    dx = a.positions[j] - a.positions[i]
    # convert to fractional, wrap
    inv = np.linalg.inv(cell)
    df = inv @ dx
    df[0] -= np.round(df[0])
    df[1] -= np.round(df[1])
    df[2] = 0  # ignore z difference
    v = cell.T @ df
    return float(np.sqrt(v[0]**2 + v[1]**2))


def is_duplicate(a1, anchor1, a2, anchor2):
    """Two structures duplicate if xy position of anchor is close."""
    # Both structures should be same shape (same slab). Use a1 for cell.
    dx = a1.positions[anchor1] - a2.positions[anchor2]
    # Wrap through cell1
    cell = a1.cell.array
    inv = np.linalg.inv(cell)
    df = inv @ dx
    df[0] -= np.round(df[0])
    df[1] -= np.round(df[1])
    df[2] = 0
    v = cell.T @ df
    return float(np.sqrt(v[0]**2 + v[1]**2))


def pick_top5(sid, ads_key):
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
    picks = []
    for r in sorted(unique, key=lambda r: r['E']):
        if not r.get('converged', True): continue
        a = traj[r['idx']]
        if len(a) != n_sub + n_ads:
            ads_at = a[-n_ads:]; a = slab.copy(); a += ads_at
        v = valid_f(a)
        if v is None: continue
        # Check duplicate vs prior picks
        dup = False
        for p in picks:
            if abs(r['E'] - p['E']) < DELTA_E_DUP:
                # Same E, now check xy
                if ads_key == 'coads':
                    d1 = is_duplicate(a, v[0], p['atoms'], p['anc'][0])
                    d2 = is_duplicate(a, v[1], p['atoms'], p['anc'][1])
                    if d1 < XY_DUP and d2 < XY_DUP:
                        dup = True; break
                else:
                    d = is_duplicate(a, v, p['atoms'], p['anc'])
                    if d < XY_DUP: dup = True; break
        if dup: continue
        picks.append({'idx': r['idx'], 'E': r['E'], 'atoms': a, 'anc': v})
        if len(picks) >= N_TARGET: break
    return picks


def fix_bottom_half(atoms, n_sub):
    z = atoms.positions[:n_sub, 2]; zm = np.median(z)
    fixed = [i for i in range(n_sub) if atoms.positions[i,2] < zm]
    atoms.set_constraint(FixAtoms(indices=fixed))


summary = []
for sid in SURFACES_SINGLE:
    slab = read(G2/SDIRS[sid]/'CONTCAR'); n_sub = len(slab)
    for ads in ['CO','CH3O']:
        picks = pick_top5(sid, ads)
        out_dir = OUT/sid/ads; out_dir.mkdir(parents=True)
        for k, p in enumerate(picks):
            atoms = p['atoms'].copy()
            fix_bottom_half(atoms, n_sub)
            fname = f"{k:02d}_{ads}_idx{p['idx']:05d}.vasp"
            write(str(out_dir/fname), atoms, format='vasp', direct=True, sort=True, vasp5=True)
            summary.append({'sid':sid,'ads':ads,'rank':k,'idx':p['idx'],'E':p['E']})
        print(f'  {sid} {ads:<6}: {len(picks)} picks')

for sid in SURFACES_COADS:
    slab = read(G2/SDIRS[sid]/'CONTCAR'); n_sub = len(slab)
    picks = pick_top5(sid, 'coads')
    out_dir = OUT/sid/'coads'; out_dir.mkdir(parents=True)
    for k, p in enumerate(picks):
        atoms = p['atoms'].copy()
        fix_bottom_half(atoms, n_sub)
        fname = f"{k:02d}_coads_idx{p['idx']:05d}.vasp"
        write(str(out_dir/fname), atoms, format='vasp', direct=True, sort=True, vasp5=True)
        summary.append({'sid':sid,'ads':'coads','rank':k,'idx':p['idx'],'E':p['E']})
    print(f'  {sid} coads : {len(picks)} picks')

json.dump(summary, open(OUT/'summary.json','w'), indent=2)
print(f'\nTotal picks: {len(summary)}')
