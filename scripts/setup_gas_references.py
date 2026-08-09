"""Gas-phase reference DFT setup for T1.18 adsorption energies.

Creates isolated-molecule single-point relaxations in a 15 Å cubic box for:
  CO, CH3O (radical, spin-polarized), CH3OH, H2

Two variants per molecule:
  {mol}_vacuum      PBE-D3, no VASPsol   ← matches vacuum-adsorbate E_bind convention
  {mol}_vaspsol     PBE-D3, VASPsol on   ← matches solvated E_bind for T1.17 output

Output: calculations/gas_references/<mol>_<variant>/
Files: POSCAR (isolated in box), INCAR, POTCAR, submit script, metadata.json
No jobs submitted.
"""
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
OUT = ROOT/'calculations/gas_references'
OUT.mkdir(exist_ok=True)
POT_LIB = Path('/home/hyunjin/POTENTIAL/potpaw_PBE')
POT_FOLDER = {'C':'C','H':'H','O':'O','Pd':'Pd_pv'}
BOX = 15.0   # Å

# ---------- molecules ----------
def make_CO():
    # linear CO along z, centered
    a = Atoms('CO', positions=[[0,0,0],[0,0,1.128]], cell=[BOX,BOX,BOX], pbc=False)
    a.center()
    return a

def make_CH3O():
    # methoxide radical (CH3-O)
    # C at origin, O along +z at 1.42 Å, 3 H tetrahedrally around C in -z hemisphere
    C = np.array([0.,0.,0.])
    O = np.array([0.,0.,1.42])
    # 3 H at tetrahedral positions; C-H = 1.09 Å; angle H-C-O ≈ 109.5°
    h_len = 1.09
    tilt = np.radians(180 - 109.5)   # from +z axis
    H1 = C + h_len*np.array([np.sin(tilt), 0, np.cos(tilt)])
    H2 = C + h_len*np.array([np.sin(tilt)*np.cos(2*np.pi/3), np.sin(tilt)*np.sin(2*np.pi/3), np.cos(tilt)])
    H3 = C + h_len*np.array([np.sin(tilt)*np.cos(4*np.pi/3), np.sin(tilt)*np.sin(4*np.pi/3), np.cos(tilt)])
    a = Atoms('COHHH', positions=[C, O, H1, H2, H3], cell=[BOX,BOX,BOX], pbc=False)
    a.center()
    return a

def make_CH3OH():
    # methanol: C-O-H bent, 3 H on C tetrahedrally
    C = np.array([0.,0.,0.])
    O = np.array([0.,0.,1.43])
    # O-H at 108.5° from C-O direction
    OH_len = 0.96
    OH_ang = np.radians(108.5)
    H_OH = O + OH_len*np.array([np.sin(OH_ang), 0, np.cos(OH_ang)])
    # methyl H's
    h_len = 1.09
    tilt = np.radians(180 - 109.5)
    H1 = C + h_len*np.array([np.sin(tilt), 0, np.cos(tilt)])
    H2 = C + h_len*np.array([np.sin(tilt)*np.cos(2*np.pi/3), np.sin(tilt)*np.sin(2*np.pi/3), np.cos(tilt)])
    H3 = C + h_len*np.array([np.sin(tilt)*np.cos(4*np.pi/3), np.sin(tilt)*np.sin(4*np.pi/3), np.cos(tilt)])
    a = Atoms('COHHHH', positions=[C, O, H_OH, H1, H2, H3], cell=[BOX,BOX,BOX], pbc=False)
    a.center()
    return a

def make_H2():
    a = Atoms('HH', positions=[[0,0,0],[0,0,0.74]], cell=[BOX,BOX,BOX], pbc=False)
    a.center()
    return a

MOLECULES = {
    'CO':   {'atoms': make_CO(),   'spin': False, 'magmom': None,
             'note': 'CO, linear'},
    # After sort=True the POSCAR order is C(1), H(3), O(1) — MAGMOM follows that order.
    'CH3O': {'atoms': make_CH3O(), 'spin': True,  'magmom': '0 0 0 0 1',
             'note': 'CH3O radical (open-shell, S=1/2), spin-polarised. '
                     'MAGMOM per POSCAR order [C H H H O]: unpaired electron on O.'},
    'CH3OH':{'atoms': make_CH3OH(),'spin': False, 'magmom': None,
             'note': 'methanol (closed-shell)'},
    'H2':   {'atoms': make_H2(),   'spin': False, 'magmom': None,
             'note': 'molecular hydrogen'},
}

# ---------- INCAR templates ----------
def incar_vacuum(spin, magmom_line):
    lines = [
        "SYSTEM = gas reference",
        "ENCUT = 520",
        "PREC = Accurate",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
        f"ISPIN = {'2' if spin else '1'}",
    ]
    if spin and magmom_line:
        lines.append(f"MAGMOM = {magmom_line}")
    lines += [
        "IVDW = 12",
        "EDIFF = 1e-06",
        "NELM = 500",
        "NELMIN = 5",
        "ALGO = Normal",
        "NCORE = 1",
        "LREAL = Auto",
        "LWAVE = .FALSE.",
        "LCHARG = .FALSE.",
        "ISYM = 0",
        "IBRION = 2",
        "NSW = 200",
        "ISIF = 2",
        "ISMEAR = 0",
        "SIGMA = 0.01",
        "EDIFFG = -0.01",
        "KSPACING = 1.0",    # effectively gamma-only for an isolated molecule
    ]
    return '\n'.join(lines) + '\n'

def incar_vaspsol(spin, magmom_line):
    base = incar_vacuum(spin, magmom_line).rstrip()
    return base + """

# --- VASPsol (methanol) ---
LSOL = .TRUE.
EB_K = 32.6
TAU = 0
LAMBDA_D_K = 3.0
"""

# ---------- submit scripts ----------
SUBMIT_VACUUM = """#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err

unset PYTHONPATH PYTHONHOME CONDA_PREFIX CONDA_DEFAULT_ENV
unset CONDA_SHLVL CONDA_PROMPT_MODIFIER
export LD_LIBRARY_PATH=""

export NVHPC=$HOME/nvhpc
export NVARCH=Linux_x86_64
export NVVERSION=25.9
export PATH=$NVHPC/$NVARCH/$NVVERSION/compilers/bin:$PATH
export PATH=$NVHPC/$NVARCH/$NVVERSION/comm_libs/mpi/bin:$PATH
export LD_LIBRARY_PATH=$NVHPC/$NVARCH/$NVVERSION/comm_libs/mpi/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVHPC/$NVARCH/$NVVERSION/compilers/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVHPC/$NVARCH/$NVVERSION/compilers/extras/qd/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/fftw/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/cuda-13.0/lib64:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1

VASP_BIN=${VASP_BIN:-/home/hyunjin/vasp.6.4.3/bin/vasp_std}
NPROCS=${SLURM_NTASKS:-1}
echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | VASP: ${VASP_BIN}"
echo "Start: $(date)"
mpirun --bind-to none -np ${NPROCS} ${VASP_BIN}
echo "End: $(date)"
"""

SUBMIT_VASPSOL = SUBMIT_VACUUM.replace(
    'VASP_BIN=${VASP_BIN:-/home/hyunjin/vasp.6.4.3/bin/vasp_std}',
    'VASP_SOL_BIN=${VASP_SOL_BIN:-/home/hyunjin/vasp.6.4.3_sol/bin/vasp_std}\nVASP_BIN=${VASP_SOL_BIN}',
)

# ---------- build ----------
def build_dir(mol_name, variant):
    m = MOLECULES[mol_name]
    a = m['atoms'].copy()
    dest = OUT/f'{mol_name}_{variant}'
    dest.mkdir(exist_ok=True)
    write(str(dest/'POSCAR'), a, format='vasp', direct=True, sort=True, vasp5=True)
    incar = incar_vacuum(m['spin'], m['magmom']) if variant=='vacuum' else incar_vaspsol(m['spin'], m['magmom'])
    (dest/'INCAR').write_text(incar)
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    potcar = ''.join((POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species)
    (dest/'POTCAR').write_text(potcar)
    submit = SUBMIT_VACUUM if variant=='vacuum' else SUBMIT_VASPSOL
    (dest/'submit_vasp.sh').write_text(submit)
    (dest/'metadata.json').write_text(json.dumps({
        'molecule': mol_name,
        'variant': variant,
        'natoms': len(a),
        'species': species,
        'box_A': BOX,
        'ISPIN': 2 if m['spin'] else 1,
        'MAGMOM': m['magmom'],
        'purpose': (f'Gas-phase {mol_name} reference for T1.18 adsorption energy '
                    f'({"solvated" if variant=="vaspsol" else "vacuum"} convention).'),
        'note': m['note'],
    }, indent=2))
    return dest

def main():
    dirs = []
    for mol in MOLECULES:
        for variant in ['vacuum', 'vaspsol']:
            d = build_dir(mol, variant)
            dirs.append(d)
            print(f'  {mol}_{variant} → {d.relative_to(ROOT)}')
    readme = """# gas_references — isolated-molecule DFT references for T1.18

Regenerate: `python scripts/setup_gas_references.py` (idempotent — overwrites files).

## Purpose

Provide DFT reference energies for adsorption energy formulae (T1.18):
```
ΔG_CO*        = G_slab+CO         − G_slab − μ_CO
ΔG_CH3O*_rad  = G_slab+CH3O       − G_slab − G_CH3O(radical)
ΔG_CH3O*_MeOH(U) = G_slab+CH3O + ½G_H2 − G_slab − G_CH3OH − eU
```

## Contents

8 dirs total (4 molecules × 2 variants):

| dir | ISPIN | MAGMOM | LSOL | note |
|---|---|---|---|---|
| CO_vacuum        | 1 | –                    | off | vacuum reference for vacuum E_bind |
| CO_vaspsol       | 1 | –                    | on  | solvated reference for T1.17 E_bind |
| CH3O_vacuum      | 2 | O 1 (radical)        | off | open-shell doublet |
| CH3O_vaspsol     | 2 | O 1 (radical)        | on  | " |
| CH3OH_vacuum     | 1 | –                    | off | closed-shell methanol |
| CH3OH_vaspsol    | 1 | –                    | on  | " |
| H2_vacuum        | 1 | –                    | off | for CHE reference |
| H2_vaspsol       | 1 | –                    | on  | " |

## Setup

- 15 Å cubic box (non-periodic PBC off — VASP still uses periodic images but
  15 Å is large enough for isolated behavior).
- KSPACING = 1.0 → effectively gamma-only for such a large cell.
- ENCUT 520, PREC=Accurate, LASPH, ADDGRID, IVDW=12 — identical to slab settings.
- ISMEAR=0, SIGMA=0.01 for isolated molecules (very sharp).
- EDIFFG = -0.01 (tighter for gas molecule geometry).

## Submitting

- **Vacuum variants**: use standard vasp_std (`VASP_BIN` in submit script).
- **VASPsol variants**: require VASPsol-enabled build (`VASP_SOL_BIN`).

## Post-processing (once complete)

Extract from OUTCAR (`grep "energy(sigma->0)" OUTCAR | tail -1`) and record in
paper_data/03_mace_references.csv as *DFT* references (currently only MACE
references are in that table).
"""
    (OUT/'README.md').write_text(readme)
    print(f'\n{len(dirs)} gas reference dirs created under {OUT.relative_to(ROOT)}')

if __name__=='__main__':
    main()
