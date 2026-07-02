"""Top-5 site-diverse DFT shortlist per (surface × ads type).
Filters: converged + intramol valid + CO2* excluded + distinct site combo + ΔE >= threshold.
Output: shortlist_v3/<surface>/<ads>/00_..._idx?????.vasp
"""
import json, shutil
from pathlib import Path
from collections import defaultdict
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
SURFACES_SINGLE = ['S1','S2','S3','S3b','S4']    # CO*, CH3O*
SURFACES_COADS  = ['S1','S2','S3b','S3']         # per advisor (no S4/PdO2)

DELTA_E_MIN = 0.02  # eV — skip candidates within 20 meV of prior pick (near-degenerate)
CO2_CUTOFF = 2.0    # C-O_lattice < this = exclude (CO2 or boundary)

def site_type(atoms, anchor, n_ads):
    syms = atoms.get_chemical_symbols()
    n_sub_atoms = len(atoms) - n_ads
    sub = list(range(n_sub_atoms))
    d = atoms.get_distances(anchor, sub, mic=True)
    nbrs = [(sub[i], d[i]) for i in range(len(sub)) if d[i] < 2.6]
    n_pd = sum(1 for i,_ in nbrs if syms[i]=='Pd')
    n_o = sum(1 for i,_ in nbrs if syms[i]=='O')
    total = n_pd + n_o
    if total == 0: return 'physi'
    if total == 1: return 'atop_Pd' if n_pd else 'atop_O'
    if total == 2:
        if n_pd == 2: return 'br_PdPd'
        if n_o == 2: return 'br_OO'
        return 'br_PdO'
    if total == 3:
        if n_pd == 3: return 'h3Pd'
        if n_o == 3: return 'h3O'
        return f'h3({n_pd}Pd{n_o}O)'
    return f'{total}f({n_pd}Pd{n_o}O)'


def valid_single_CO(atoms):
    syms = atoms.get_chemical_symbols()
    c_idx = [i for i,s in enumerate(syms) if s=='C']
    o_idx = [i for i,s in enumerate(syms) if s=='O']
    if len(c_idx) != 1: return None
    c = c_idx[0]
    d_o = sorted([(atoms.get_distance(c,oi,mic=True), oi) for oi in o_idx])
    if not (1.05 <= d_o[0][0] <= 1.30): return None
    if len(d_o) >= 2 and d_o[1][0] < CO2_CUTOFF: return None
    return c   # anchor = C


def valid_single_CH3O(atoms):
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
    return d_o[0][1]   # anchor = methoxy O


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
    # CO validity
    d_o_co = sorted([(atoms.get_distance(co_c,oi,mic=True), oi) for oi in o_idx])
    if not (1.05 <= d_o_co[0][0] <= 1.30): return None
    if len(d_o_co) >= 2 and d_o_co[1][0] < CO2_CUTOFF: return None
    # CH3O validity
    d_ch = [atoms.get_distance(me_c,h,mic=True) for h in h_idx]
    if not all(0.90<=d<=1.25 for d in d_ch): return None
    d_o_me = sorted([(atoms.get_distance(me_c,oi,mic=True), oi) for oi in o_idx])
    if not (1.30 <= d_o_me[0][0] <= 1.55): return None
    if len(d_o_me) >= 2 and d_o_me[1][0] < CO2_CUTOFF: return None
    return (co_c, d_o_me[0][1])   # anchors (CO's C, methoxy's O)


def fix_bottom_half(atoms, n_sub):
    z = atoms.positions[:n_sub, 2]; zm = np.median(z)
    fixed = [i for i in range(n_sub) if atoms.positions[i,2] < zm]
    atoms.set_constraint(FixAtoms(indices=fixed))


def pick_top5(sid, ads_key, valid_func, get_site_labels):
    """ads_key: 'CO','CH3O','coads'. Returns picks list."""
    sdir = SDIRS[sid]
    slab = read(G2/sdir/'CONTCAR'); n_sub = len(slab)
    if ads_key == 'coads':
        n_ads = 7
        unique = json.load(open(G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'))
        traj = list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
    else:
        n_ads = 2 if ads_key == 'CO' else 5
        unique = json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads_key}.json'))
        traj = list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads_key}.traj', index=':'))
    picks = []
    combos_seen = set()
    for r in sorted(unique, key=lambda r: r['E']):
        if not r.get('converged', True): continue
        a = traj[r['idx']]
        if len(a) != n_sub + n_ads:
            ads_at = a[-n_ads:]; a = slab.copy(); a += ads_at
        v = valid_func(a)
        if v is None: continue
        combo = get_site_labels(a, v, n_ads)
        if combo in combos_seen: continue
        # ΔE from most recent pick
        if picks and abs(r['E'] - picks[-1]['E']) < DELTA_E_MIN:
            continue
        combos_seen.add(combo)
        picks.append({'idx': r['idx'], 'E': r['E'], 'combo': combo, 'atoms': a})
        if len(picks) >= 5: break
    return picks


def site_labels_CO(atoms, anchor, n_ads):
    return (site_type(atoms, anchor, n_ads),)

def site_labels_CH3O(atoms, anchor, n_ads):
    return (site_type(atoms, anchor, n_ads),)

def site_labels_coads(atoms, anchors, n_ads):
    co_c, me_o = anchors
    return (site_type(atoms, co_c, n_ads), site_type(atoms, me_o, n_ads))


# === Main ===
summary = []

for sid in SURFACES_SINGLE:
    slab = read(G2/SDIRS[sid]/'CONTCAR'); n_sub = len(slab)
    for ads, valid_f, sl_f in [
        ('CO', valid_single_CO, site_labels_CO),
        ('CH3O', valid_single_CH3O, site_labels_CH3O),
    ]:
        picks = pick_top5(sid, ads, valid_f, sl_f)
        out_dir = OUT/sid/ads; out_dir.mkdir(parents=True)
        for k, p in enumerate(picks):
            atoms = p['atoms'].copy()
            fix_bottom_half(atoms, n_sub)
            combo_str = '_'.join(p['combo'])
            fname = f"{k:02d}_{ads}_{combo_str}_idx{p['idx']:05d}.vasp"
            write(str(out_dir/fname), atoms, format='vasp', direct=True, sort=True, vasp5=True)
            summary.append({'sid':sid,'ads':ads,'rank':k,'idx':p['idx'],'E':p['E'],'combo':combo_str})
        print(f'  {sid} {ads}: {len(picks)} picks')

for sid in SURFACES_COADS:
    slab = read(G2/SDIRS[sid]/'CONTCAR'); n_sub = len(slab)
    picks = pick_top5(sid, 'coads', valid_coads, site_labels_coads)
    out_dir = OUT/sid/'coads'; out_dir.mkdir(parents=True)
    for k, p in enumerate(picks):
        atoms = p['atoms'].copy()
        fix_bottom_half(atoms, n_sub)
        combo_str = f"CO-{p['combo'][0]}_OMe-{p['combo'][1]}"
        fname = f"{k:02d}_coads_{combo_str}_idx{p['idx']:05d}.vasp"[:80]
        write(str(out_dir/fname), atoms, format='vasp', direct=True, sort=True, vasp5=True)
        summary.append({'sid':sid,'ads':'coads','rank':k,'idx':p['idx'],'E':p['E'],'combo':combo_str})
    print(f'  {sid} coads: {len(picks)} picks')

# Save summary
json.dump(summary, open(OUT/'summary.json','w'), indent=2)
print(f'\nTotal picks: {len(summary)}')
print(f'Saved to {OUT}')
