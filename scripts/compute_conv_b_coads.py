"""Convention B binding energy for co-adsorption — frozen post-relax slab single-point.

E_bind^B(coads) = E(slab+CO+CH3O) − E(slab^{post-relax}, no ads) − E(CO_gas) − E(CH3OH_gas − 1/2 H2_gas)
Save: convB_binding/<surface>_coads_convB.json
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
E_CH3O_thermo = refs['gas']['CH3O_ref']   # CH3OH - 1/2 H2 (advisor's coads convention)
E_CLEAN = refs['slab']

print('Loading MACE...')
calc = mace_mp(model='mh-1', head='oc20_usemppbe',
               default_dtype='float64', enable_cueq=True, device='cuda',
               dispersion=True, damping='bj', dispersion_xc='pbe')

def single_point(atoms):
    a = atoms.copy(); a.calc = calc
    return float(a.get_potential_energy())

n_ads = 7   # 2 (CO) + 5 (CH3O)
t0 = time.time()
for sid, sdir in SDIRS.items():
    slab = read(G2/sdir/'CONTCAR'); n_sub = len(slab)
    unique = json.load(open(G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'))
    traj = list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
    out_recs = []
    for i, r in enumerate(unique):
        a = traj[r['idx']]
        if len(a) != n_sub + n_ads:
            ads_a = a[-n_ads:]; a = slab.copy(); a += ads_a
        slab_post = a[:n_sub]
        E_slab_post = single_point(slab_post)
        E_A = r['E'] - E_CLEAN[sid] - E_CO - E_CH3O_thermo
        E_B = r['E'] - E_slab_post - E_CO - E_CH3O_thermo
        out_recs.append({'idx': r['idx'], 'E_total': r['E'],
                          'E_slab_post': E_slab_post,
                          'E_bind_A_clean': E_A,
                          'E_bind_B_frozen': E_B,
                          'delta_recon': E_A - E_B})
        if (i+1) % 500 == 0:
            print(f'  {sid} {i+1}/{len(unique)}, {time.time()-t0:.0f}s')
    json.dump(out_recs, open(OUT/f'{sid}_coads_convB.json','w'), indent=2)
    print(f'{sid}: {len(out_recs)} cands done, elapsed {time.time()-t0:.0f}s')

print(f'\nDONE in {time.time()-t0:.0f}s')
