"""Convention B binding energy — frozen post-relax slab single-point.

For each MLIP-relaxed (slab+ads) candidate:
  E_slab^{post-relax} = MACE single-point on (atoms - ads)
  E_bind^B = E(slab+ads) − E_slab^{post-relax} − E(ads_gas)

Saves: convB_binding_<surface>_<ads>.json with per-candidate dict.
"""
import os, json, time
os.environ['CUDA_VISIBLE_DEVICES']='0'
import warnings; warnings.filterwarnings('ignore')

import numpy as np
from pathlib import Path
from ase.io import read
from mace.calculators import mace_mp

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT/'calculations/G3_adsorption'
G2 = ROOT/'calculations/G2_slab'
OUT = G3/'convB_binding'
OUT.mkdir(exist_ok=True)

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

refs = json.load(open(G3/'mace_d3_references.json'))
E_CO = refs['gas']['CO']
E_CH3O = refs['gas']['CH3O_radical']

print('Loading MACE...')
calc = mace_mp(model='mh-1', head='oc20_usemppbe',
               default_dtype='float64', enable_cueq=True, device='cuda',
               dispersion=True, damping='bj', dispersion_xc='pbe')

def single_point(atoms):
    a = atoms.copy()
    a.calc = calc
    return float(a.get_potential_energy())

t0 = time.time()
for sid, sdir in SDIRS.items():
    slab = read(G2/sdir/'CONTCAR')
    n_sub = len(slab)
    for ads, n_ads in [('CO',2), ('CH3O',5)]:
        unique = json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
        traj = list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
        out_recs = []
        for r in unique:
            a = traj[r['idx']]
            if len(a) != n_sub + n_ads:
                ads_atoms = a[-n_ads:]
                a = slab.copy(); a += ads_atoms
            # frozen slab: just take first n_sub atoms (substrate)
            slab_post = a[:n_sub]
            E_slab_post = single_point(slab_post)
            ref = E_CO if ads=='CO' else E_CH3O
            E_bind_A = r['E'] - (-389.5571 if sid=='S1' else -559.9475 if sid=='S2' else -673.2548 if sid=='S3' else -531.0268 if sid=='S3b' else -738.9077) - ref
            E_bind_B = r['E'] - E_slab_post - ref
            out_recs.append({
                'idx': r['idx'],
                'E_total': r['E'],
                'E_slab_post': E_slab_post,
                'E_bind_A_clean': E_bind_A,
                'E_bind_B_frozen': E_bind_B,
                'delta_recon': E_bind_A - E_bind_B,
            })
        out_f = OUT / f'{sid}_{ads}_convB.json'
        json.dump(out_recs, open(out_f, 'w'), indent=2)
        print(f'  {sid} {ads}: {len(out_recs)} cands, {time.time()-t0:.0f}s elapsed')

print(f'\nDONE in {time.time()-t0:.0f}s')
