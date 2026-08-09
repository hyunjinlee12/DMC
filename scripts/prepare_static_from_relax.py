"""Populate paired-static gas dirs (CH3O_vacuum_{15,20}A_static and
CH3OH_vacuum_{15,20}A_static) with the SAME relaxed atomic positions,
but each dir keeps its own box size.

Reviewer round 3: copying CONTCAR wholesale (with its cell) into a
20 Å dir would revert the cell to 15 Å and invalidate the finite-cell
test. This script extracts positions only, re-centers them in the
target box, and writes POSCAR with the correct cell.

Idempotent. Only writes into dirs whose source _vacuum relaxation has
completed (CONTCAR present with recognisable energy).
"""
import subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
GAS = ROOT/'calculations/gas_references'

PAIRS = [
    ('CH3O',  [15.0, 20.0]),
    ('CH3OH', [15.0, 20.0]),
]

def is_relaxed(contcar):
    if not contcar.exists(): return False
    # A CONTCAR from a converged relax will have "reached required accuracy" in OUTCAR
    outcar = contcar.parent/'OUTCAR'
    if not outcar.exists(): return False
    return 'reached required accuracy' in outcar.read_text()

def prepare(mol):
    src_dir = GAS/f'{mol}_vacuum'
    src_contcar = src_dir/'CONTCAR'
    if not is_relaxed(src_contcar):
        return {'mol':mol,'status':'skipped','reason':'source relax not yet converged'}
    atoms = read(src_contcar)     # has 15 Å cell + relaxed coords
    results = []
    for box in [15.0, 20.0]:
        dest_dir = GAS/f'{mol}_vacuum_{int(box)}A_static'
        if not dest_dir.exists():
            results.append({'mol':mol,'box':box,'status':'skipped','reason':'dest dir missing'})
            continue
        a = atoms.copy()
        a.set_cell([box, box, box], scale_atoms=False)   # cell only, NOT positions
        a.set_pbc(True)
        a.center()
        # POTCAR order + selective dynamics: no selective dynamics for isolated
        # molecules; just species-grouped POSCAR via sort=True
        write(str(dest_dir/'POSCAR'), a, format='vasp',
              direct=True, vasp5=True, sort=True)
        # sanity: POSCAR cell diagonal is exactly `box`
        head = (dest_dir/'POSCAR').read_text().splitlines()
        cell_lines = head[2:5]
        vals = [float(cell_lines[i].split()[i]) for i in range(3)]
        if not all(abs(v - box) < 1e-4 for v in vals):
            results.append({'mol':mol,'box':box,'status':'FAIL',
                           'reason':f'cell diag {vals} != {box}'})
            continue
        results.append({'mol':mol,'box':box,'status':'written',
                       'natoms':len(a),
                       'cell_A': box})
    return results

def main():
    all_results = []
    for mol, _ in PAIRS:
        r = prepare(mol)
        all_results.extend(r if isinstance(r, list) else [r])
    for r in all_results:
        print(f'  {r}')
    n_ok = sum(1 for r in all_results if r.get('status')=='written')
    n_pending = sum(1 for r in all_results if r.get('status')=='skipped')
    n_fail = sum(1 for r in all_results if r.get('status')=='FAIL')
    print(f'\nPaired-static POSCARs: {n_ok} written, {n_pending} pending, {n_fail} FAIL')
    if n_fail:
        print('FAIL — do not submit until resolved.')

if __name__=='__main__':
    main()
