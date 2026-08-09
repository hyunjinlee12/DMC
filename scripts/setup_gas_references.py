"""Gas-phase reference DFT setup for T1.18 adsorption energies.

Post-review revision (2026-07-17):
- Removed LAMBDA_D_K (electrolyte screening, not needed for neutral solvent).
- CH3O radical: added NUPDOWN=1 (enforce total moment = 1 μ_B).
- Explicit Γ-only KPOINTS file for isolated molecules (unambiguous).
- Renamed submit script env var from VASP_SOL_BIN to VASP_BIN with note
  that VASP ≥5.4.1 standard builds usually support LSOL.
- MAGMOM for CH3O verified against ASE sort=True POSCAR order (C H H H O).
- Optional 20 Å box variants for polar molecules (CH3OH, CH3O) to check
  finite-cell error — build_dir() called with box_A parameter.

Creates 8 primary dirs (4 molecules × {vacuum, vaspsol}) + 2 sanity-check
20 Å variants for polar molecules (CH3OH_vacuum_20A, CH3O_vacuum_20A).

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

def make_CO():
    a = Atoms('CO', positions=[[0,0,0],[0,0,1.128]], pbc=False)
    return a
def make_CH3O():
    C = np.array([0.,0.,0.]); O = np.array([0.,0.,1.42])
    h_len = 1.09; tilt = np.radians(180 - 109.5)
    H1 = C + h_len*np.array([np.sin(tilt), 0, np.cos(tilt)])
    H2 = C + h_len*np.array([np.sin(tilt)*np.cos(2*np.pi/3), np.sin(tilt)*np.sin(2*np.pi/3), np.cos(tilt)])
    H3 = C + h_len*np.array([np.sin(tilt)*np.cos(4*np.pi/3), np.sin(tilt)*np.sin(4*np.pi/3), np.cos(tilt)])
    return Atoms('COHHH', positions=[C, O, H1, H2, H3], pbc=False)
def make_CH3OH():
    C = np.array([0.,0.,0.]); O = np.array([0.,0.,1.43])
    OH_len = 0.96; OH_ang = np.radians(108.5)
    H_OH = O + OH_len*np.array([np.sin(OH_ang), 0, np.cos(OH_ang)])
    h_len = 1.09; tilt = np.radians(180 - 109.5)
    H1 = C + h_len*np.array([np.sin(tilt), 0, np.cos(tilt)])
    H2 = C + h_len*np.array([np.sin(tilt)*np.cos(2*np.pi/3), np.sin(tilt)*np.sin(2*np.pi/3), np.cos(tilt)])
    H3 = C + h_len*np.array([np.sin(tilt)*np.cos(4*np.pi/3), np.sin(tilt)*np.sin(4*np.pi/3), np.cos(tilt)])
    return Atoms('COHHHH', positions=[C, O, H_OH, H1, H2, H3], pbc=False)
def make_H2():
    return Atoms('HH', positions=[[0,0,0],[0,0,0.74]], pbc=False)

# spin + magmom: MAGMOM is per POSCAR-order after sort=True (ASE sorts by symbol
# alphabetically: C→H→O). For CH3O the order is C(1) H(3) O(1), so the last
# entry corresponds to O — the atom that carries the unpaired electron.
MOLECULES = {
    'CO':    {'atoms_fn': make_CO,    'spin': False, 'nupdown': None, 'magmom': None,
              'note': 'CO, linear (closed-shell singlet)'},
    'CH3O':  {'atoms_fn': make_CH3O,  'spin': True,  'nupdown': 1,     'magmom': '0 0 0 0 1',
              'note': ('CH3O radical (open-shell doublet, S=1/2). '
                       'POSCAR sort=True order = C H H H O; MAGMOM places unpaired μ on O. '
                       'NUPDOWN=1 enforces total moment = 1 μ_B.')},
    'CH3OH': {'atoms_fn': make_CH3OH, 'spin': False, 'nupdown': None, 'magmom': None,
              'note': ('methanol, closed-shell singlet. '
                       'Convention note: G(CH3OH) here is the ISOLATED-GAS reference; '
                       'if the T1.18 formula requires a liquid/solution standard state, '
                       'a Δμ_solvation correction must be added post-hoc.')},
    'H2':    {'atoms_fn': make_H2,    'spin': False, 'nupdown': None, 'magmom': None,
              'note': 'molecular hydrogen (for ½ E(H2) in CHE reference)'},
}

def incar_common(spin, nupdown, magmom, solvate):
    lines = [
        "SYSTEM = gas reference (T1.18 μ_gas)",
        "ENCUT = 520",
        "PREC = Accurate",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
        f"ISPIN = {2 if spin else 1}",
    ]
    if spin:
        if nupdown is not None:
            lines.append(f"NUPDOWN = {nupdown}")
        if magmom:
            lines.append(f"MAGMOM = {magmom}")
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
        # KPOINTS handled via explicit file (Γ-only); do NOT set KSPACING here
        # to avoid mesh ambiguity.
    ]
    if solvate:
        lines += ["",
                  "# ---- VASPsol (methanol implicit solvent) ----",
                  "LSOL = .TRUE.",
                  "EB_K = 32.6",
                  "TAU = 0"]
    return '\n'.join(lines) + '\n'

KPOINTS_GAMMA = """Gamma-only
0
Gamma
1 1 1
0 0 0
"""

SUBMIT = """#!/bin/bash
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

# VASP ≥5.4.1 standard builds typically include LSOL support. If this run's
# OUTCAR shows "unknown INCAR tag" for LSOL/EB_K/TAU, rebuild with solvation.
VASP_BIN=${VASP_BIN:-/home/hyunjin/vasp.6.4.3/bin/vasp_std}
NPROCS=${SLURM_NTASKS:-1}
echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | VASP: ${VASP_BIN}"
echo "Start: $(date)"
mpirun --bind-to none -np ${NPROCS} ${VASP_BIN}
echo "End: $(date)"
"""

def build_dir(mol_name, variant, box_A=15.0, suffix=''):
    m = MOLECULES[mol_name]
    a = m['atoms_fn']()
    a.set_cell([box_A, box_A, box_A]); a.center()
    solvate = (variant == 'vaspsol')
    dest = OUT/f'{mol_name}_{variant}{suffix}'
    dest.mkdir(exist_ok=True)
    write(str(dest/'POSCAR'), a, format='vasp', direct=True, sort=True, vasp5=True)
    (dest/'INCAR').write_text(incar_common(m['spin'], m['nupdown'], m['magmom'], solvate))
    (dest/'KPOINTS').write_text(KPOINTS_GAMMA)
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    potcar = ''.join((POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species)
    (dest/'POTCAR').write_text(potcar)
    (dest/'submit_vasp.sh').write_text(SUBMIT)
    (dest/'metadata.json').write_text(json.dumps({
        'molecule': mol_name,
        'variant': variant,
        'box_A': box_A,
        'natoms': len(a),
        'poscar_species_order': species,
        'ISPIN': 2 if m['spin'] else 1,
        'NUPDOWN': m['nupdown'],
        'MAGMOM': m['magmom'],
        'MAGMOM_atom_by_atom_note':
            f'MAGMOM entries correspond to POSCAR sort=True order = {species}.'
            if m['magmom'] else None,
        'LSOL': solvate,
        'KPOINTS': 'Gamma-only (explicit file)',
        'purpose': (f'Gas-phase {mol_name} reference for T1.18 adsorption energy '
                    f'({"solvated" if solvate else "vacuum"} convention).'),
        'note': m['note'],
    }, indent=2))
    return dest

def main():
    dirs = []
    # 8 primary
    for mol in MOLECULES:
        for variant in ['vacuum', 'vaspsol']:
            d = build_dir(mol, variant, box_A=15.0)
            dirs.append(d)
            print(f'  {d.name}')
    # 2 box-size sanity checks for polar molecules (Δμ_finite_cell < 5 meV target)
    for mol in ['CH3OH','CH3O']:
        d = build_dir(mol, 'vacuum', box_A=20.0, suffix='_20A')
        dirs.append(d)
        print(f'  {d.name}  (20 Å sanity check)')

    readme = """# gas_references — isolated-molecule DFT (T1.18 references)

Regenerate: `python scripts/setup_gas_references.py` (idempotent).

## Purpose

Provides reference energies for T1.18 adsorption-energy formulae:

```
ΔG_CO*(vac)         = G(slab+CO)_vac    − G(slab)_vac    − μ_CO(gas, vacuum)
ΔG_CO*(sol)         = G(slab+CO)_sol    − G(slab)_sol    − μ_CO(gas, solvated)

ΔG_CH3O*_rad(vac)   = G(slab+CH3O)_vac  − G(slab)_vac    − G_CH3O(radical, vacuum)
ΔG_CH3O*_rad(sol)   = G(slab+CH3O)_sol  − G(slab)_sol    − G_CH3O(radical, vaspsol)

ΔG_CH3O*_MeOH(U) = G(slab+CH3O) + ½G(H2) − G(slab) − G(CH3OH) − eU
```

**Never mix vacuum-slab and solvated-adsorbate energies in one formula.**
Always use references from the same phase (vacuum ↔ vacuum, solvated ↔ solvated).

## Contents

10 dirs (4 molecules × 2 variants + 2 sanity checks):

| dir | ISPIN | NUPDOWN | MAGMOM (POSCAR order) | LSOL | note |
|---|---|---|---|---|---|
| CO_vacuum        | 1 | – | – | off | vacuum μ_CO |
| CO_vaspsol       | 1 | – | – | on  | solvated μ_CO |
| CH3O_vacuum      | 2 | 1 | 0 0 0 0 1 (C H H H O) | off | radical |
| CH3O_vaspsol     | 2 | 1 | 0 0 0 0 1 (C H H H O) | on  | radical |
| CH3OH_vacuum     | 1 | – | – | off | closed-shell |
| CH3OH_vaspsol    | 1 | – | – | on  | closed-shell |
| H2_vacuum        | 1 | – | – | off | for CHE ½E(H2) |
| H2_vaspsol       | 1 | – | – | on  | for CHE ½E(H2) |
| CH3OH_vacuum_20A | 1 | – | – | off | box-size sanity (Δ vs 15 Å) |
| CH3O_vacuum_20A  | 2 | 1 | 0 0 0 0 1 | off | box-size sanity |

## Post-review changes (2026-07-17)

- Removed `LAMBDA_D_K` (belongs to electrolyte model, not neutral solvent).
- Added `NUPDOWN=1` to CH3O_* (enforces open-shell doublet).
- Explicit Γ-only KPOINTS file (safer than KSPACING for isolated molecules).
- Renamed `VASP_SOL_BIN` → `VASP_BIN` (VASP ≥5.4.1 standard builds usually
  include LSOL support; pilot-check on first run rather than assuming a
  separate binary is required).
- Added 20 Å variants for polar molecules to quantify finite-cell error.

## After computation

Extract from OUTCAR:
```bash
grep "energy(sigma->0)" OUTCAR | tail -1
grep -E "LSOL|EB_K|VASPsol|number of dipole" OUTCAR   # verify solvation active
grep "NUPDOWN" OUTCAR                                  # verify enforced moment
```

Record the DFT reference energies into `paper_data/03_mace_references.csv`
under a new `E_DFT_vacuum` / `E_DFT_vaspsol` column pair (currently only
MACE references are in that CSV).
"""
    (OUT/'README.md').write_text(readme)
    print(f'\n{len(dirs)} gas reference dirs under {OUT.relative_to(ROOT)}')

if __name__=='__main__':
    main()
