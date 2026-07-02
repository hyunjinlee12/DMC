"""Compute MACE+D3 reference energies for proper binding energy computation.

Outputs all energies with the SAME calculator (MACE-MH+D3+cueq) so that:
  E_bind(CO*)    = E(slab+CO)    - E(slab)        - E(CO_gas)
  E_bind(CH3O*)  = E(slab+CH3O)  - E(slab)        - E(CH3O_ref)
                                                    where E(CH3O_ref) = E(CH3OH_gas) - 1/2 E(H2_gas)
  E_bind(co-ads) = E(slab+CO+CH3O) - E(slab) - E(CO_gas) - E(CH3O_ref)

Saves to: calculations/G3_adsorption/mace_d3_references.json
"""
import os, json
os.environ['CUDA_VISIBLE_DEVICES'] = '0'   # use idle GPU
import warnings; warnings.filterwarnings('ignore')

from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from mace.calculators import mace_mp

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab'
OUT_F = ROOT / 'calculations/G3_adsorption/mace_d3_references.json'

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}


def make_calc():
    return mace_mp(model='mh-1', head='oc20_usemppbe',
                   default_dtype='float64', enable_cueq=True, device='cuda',
                   dispersion=True, damping='bj', dispersion_xc='pbe')


def fix_bottom_half(atoms):
    z = atoms.positions[:, 2]
    z_med = np.median(z)
    fixed = [i for i in range(len(atoms)) if atoms.positions[i, 2] < z_med]
    atoms.set_constraint(FixAtoms(indices=fixed))


def relax(atoms, calc, fmax=0.03, steps=200, name=''):
    atoms = atoms.copy()
    atoms.calc = calc
    e_init = atoms.get_potential_energy()
    opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=fmax, steps=steps)
    e_final = atoms.get_potential_energy()
    print(f"  {name:15s} E_init={e_init:.4f}  E_final={e_final:.4f}  nsteps={opt.nsteps}  conv={opt.converged()}")
    return float(e_final)


def main():
    calc = make_calc()
    refs = {'slab': {}, 'gas': {}}

    # ---- 1) Clean slab energies (MACE+D3, fixed bottom half = same as ads) ----
    print("=== Clean slabs (MACE+D3, fixed bottom-half) ===")
    for sid, sdir in SDIRS.items():
        slab = read(G2 / sdir / 'CONTCAR')
        fix_bottom_half(slab)
        e = relax(slab, calc, name=f'{sid}_slab')
        refs['slab'][sid] = e

    # ---- 2) Gas-phase references ----
    print("\n=== Gas-phase molecules (10 Å vacuum box) ===")
    # CO gas
    co = Atoms('CO', positions=[[0,0,0],[1.143,0,0]], cell=[15,15,15], pbc=True)
    co.center()
    e_co = relax(co, calc, name='CO_gas')
    refs['gas']['CO'] = e_co

    # CH3OH gas (for methoxy reference)
    ch3oh = Atoms('COH4', positions=[
        [0.0, 0.0, 0.0],          # C
        [0.0, 0.0, 1.42],          # O (bonded to C, OH side)
        [1.03, 0.0, -0.36],        # H1 on C
        [-0.51, 0.89, -0.36],      # H2 on C
        [-0.51,-0.89, -0.36],      # H3 on C
        [0.94, 0.0, 1.72],         # H on O (methanol H)
    ], cell=[15,15,15], pbc=True)
    ch3oh.center()
    e_ch3oh = relax(ch3oh, calc, name='CH3OH_gas')
    refs['gas']['CH3OH'] = e_ch3oh

    # H2 gas
    h2 = Atoms('H2', positions=[[0,0,0],[0.741,0,0]], cell=[15,15,15], pbc=True)
    h2.center()
    e_h2 = relax(h2, calc, name='H2_gas')
    refs['gas']['H2'] = e_h2

    # Derived methoxy radical reference
    refs['gas']['CH3O_ref'] = e_ch3oh - 0.5 * e_h2
    print(f"\n  CH3O_ref(thermo) = E(CH3OH) - 1/2 E(H2) = {refs['gas']['CH3O_ref']:.4f} eV")

    # Direct CH3O radical (methoxy)
    ch3o = Atoms('COH3', positions=[
        [0.0, 0.0, 0.0],          # C
        [0.0, 0.0, 1.36],          # O
        [1.03, 0.0, -0.36],        # H1
        [-0.51, 0.89, -0.36],      # H2
        [-0.51,-0.89, -0.36],      # H3
    ], cell=[15,15,15], pbc=True)
    ch3o.center()
    e_ch3o = relax(ch3o, calc, name='CH3O_radical')
    refs['gas']['CH3O_radical'] = e_ch3o

    # ---- Save ----
    OUT_F.parent.mkdir(parents=True, exist_ok=True)
    json.dump(refs, open(OUT_F, 'w'), indent=2)
    print(f"\nSaved: {OUT_F}")

    # ---- Print summary table ----
    print("\n=== Summary ===")
    print("Surface  E_slab (eV)")
    for sid in SDIRS:
        print(f"  {sid:4s}   {refs['slab'][sid]:.4f}")
    print("\nGas references:")
    for k, v in refs['gas'].items():
        print(f"  {k:10s} {v:.4f}")


if __name__ == '__main__':
    main()
