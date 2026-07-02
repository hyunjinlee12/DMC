"""Re-relax unconverged unique candidates with more steps."""
import os, json, time
os.environ['CUDA_VISIBLE_DEVICES']='0'
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from pathlib import Path
from ase.io import read, Trajectory
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from mace.calculators import mace_mp

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT/'calculations/G3_adsorption'
G2 = ROOT/'calculations/G2_slab'
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

print('Loading MACE...')
calc = mace_mp(model='mh-1', head='oc20_usemppbe',
               default_dtype='float64', enable_cueq=True, device='cuda',
               dispersion=True, damping='bj', dispersion_xc='pbe')

def fix_bottom(atoms, n_sub):
    z = atoms.positions[:n_sub, 2]; zm = np.median(z)
    fixed = [i for i in range(n_sub) if atoms.positions[i,2]<zm]
    atoms.set_constraint(FixAtoms(indices=fixed))

t_all = time.time()
for sid, sdir in SDIRS.items():
    slab = read(G2/sdir/'CONTCAR'); n_sub=len(slab)
    for ads in ['CO','CH3O']:
        n_ads = 2 if ads=='CO' else 5
        unique_f = G3/sdir/f'MLIP_phase1/unique_{ads}.json'
        traj_f = G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj'
        unique = json.load(open(unique_f))
        n_uncv = sum(1 for r in unique if not r.get('converged', True))
        if n_uncv == 0:
            print(f'  {sid} {ads}: all converged, skip')
            continue
        traj_list = list(read(traj_f, index=':'))
        updated = 0
        for i, r in enumerate(unique):
            if r.get('converged', True): continue
            a = traj_list[r['idx']]
            if len(a) != n_sub + n_ads:
                ads_a = a[-n_ads:]; a = slab.copy(); a += ads_a
            fix_bottom(a, n_sub)
            a.calc = calc
            opt = LBFGS(a, logfile=None)
            opt.run(fmax=0.05, steps=400)
            r['E'] = float(a.get_potential_energy())
            r['converged'] = bool(opt.converged())
            r['n_steps'] = int(opt.nsteps) + r.get('n_steps', 0)
            traj_list[r['idx']] = a   # update traj entry
            if opt.converged(): updated += 1
        # write back json
        json.dump(unique, open(unique_f, 'w'), indent=2)
        # write back traj (rebuild)
        w = Trajectory(str(traj_f), 'w')
        for at in traj_list: w.write(at)
        w.close()
        print(f'  {sid} {ads}: re-tried {n_uncv}, newly converged {updated}, elapsed {time.time()-t_all:.0f}s')

print(f'\nDONE in {time.time()-t_all:.0f}s')
