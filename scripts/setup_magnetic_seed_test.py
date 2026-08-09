"""Magnetic seed test — resolve initialization-vs-physical uncertainty.

Recipe from reviewer 2026-07-17 round 4:

  S1 clean + S1/CO idx=3 (both even-electron):
    ispin1              ISPIN=1, no MAGMOM
    lowmag              ISPIN=2, MAGMOM Pd=0.2, C/O=0
    highmag             ISPIN=2, MAGMOM Pd=1.0, C/O=0

  S3 clean + S3/CH3O idx=315 (S3 clean is even-electron; CH3O radical is odd):
    S3_clean:
      nonmag            ISPIN=1
      fm                ISPIN=2, MAGMOM Pd=0.5, O=0
      afm               ISPIN=2, MAGMOM alternating Pd=+0.5/-0.5 by fractional x
    S3/CH3O:
      doublet           ISPIN=2, NUPDOWN=1, MAGMOM substrate=0, adsorbate O=1
      unconstrained     ISPIN=2, MAGMOM Pd=0.3, adsorbate O=1

All: IBRION=-1, NSW=0 (STATIC). ISTART=0, ICHARG=2 (no WAVECAR reuse).
Geometry: T1.16 top-1 CONTCAR (ads) or G2 CONTCAR (clean).
LREAL and everything else identical to T1.17 for direct comparability.

Aim: judge whether each system's converged magnetic state is stable across
seeds (physical) or bounces around (metastable / artifact).
"""
import json
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
OUT = ROOT/'calculations/magnetic_seed_test'
OUT.mkdir(exist_ok=True)
G2 = ROOT/'calculations/G2_slab'
SDIRS = {'S1':'S1_Pd100','S3':'S3_PdO100'}
POT_LIB = Path('/home/hyunjin/POTENTIAL/potpaw_PBE')
POT_FOLDER = {'C':'C','H':'H','O':'O','Pd':'Pd_pv'}

INCAR_TEMPLATE = """SYSTEM = magnetic seed test (static single-point)
ENCUT = 520
PREC = Accurate
LASPH = .TRUE.
ADDGRID = .TRUE.
{spin_block}
IVDW = 12
EDIFF = 1e-06
NELM = 500
NELMIN = 5
ALGO = Normal
NCORE = 1
LREAL = Auto
LWAVE = .FALSE.
LCHARG = .FALSE.
LORBIT = 11
ISYM = 0
IBRION = -1
NSW = 0
ISMEAR = {ismear}
SIGMA = {sigma}
{ldipol_block}KSPACING = 0.25
ISTART = 0
ICHARG = 2
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
VASP_BIN=${VASP_BIN:-/home/hyunjin/vasp.6.4.3/bin/vasp_std}
NPROCS=${SLURM_NTASKS:-1}
mpirun --bind-to none -np ${NPROCS} ${VASP_BIN}
"""

def _spin_block(ispin, magmom=None, nupdown=None):
    lines = [f'ISPIN = {ispin}']
    if magmom: lines.append(f'MAGMOM = {magmom}')
    if nupdown is not None: lines.append(f'NUPDOWN = {nupdown}')
    return '\n'.join(lines)

def _magmom_for(species_order, atom_counts_per_species, per_species_mom):
    """Return VASP MAGMOM string. species_order = POSCAR line-6 order.
    per_species_mom: dict species → moment or callable(atom_index_within_species) → moment.
    """
    parts = []
    for sp, n in zip(species_order, atom_counts_per_species):
        m = per_species_mom.get(sp, 0.0)
        if callable(m):
            for i in range(n): parts.append(f'{m(i):.3f}')
        else:
            parts.append(f'{n}*{m:.3f}')
    return ' '.join(parts)

def build_seed(source_atoms, name, ismear, sigma, has_dipole, spin_block):
    dest = OUT/name
    dest.mkdir(exist_ok=True)
    write(str(dest/'POSCAR'), source_atoms, format='vasp', direct=True, sort=True, vasp5=True)
    ldipol = 'LDIPOL = .TRUE.\nIDIPOL = 3\n' if has_dipole else ''
    (dest/'INCAR').write_text(INCAR_TEMPLATE.format(
        spin_block=spin_block, ismear=ismear, sigma=sigma, ldipol_block=ldipol))
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    (dest/'POTCAR').write_text(''.join(
        (POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species))
    (dest/'submit_vasp.sh').write_text(SUBMIT)
    return dest

def afm_magmom_S3clean(atoms):
    """Split Pd atoms into two groups by fractional-x median → +/- 0.5.
    Returns MAGMOM in POSCAR-sort=True order (species-grouped: O, Pd for S3)."""
    ss = atoms.get_chemical_symbols()
    n_O = ss.count('O')
    n_Pd = ss.count('Pd')
    # after sort=True: O group first, then Pd group
    pd_indices_after_sort = list(range(n_O, n_O + n_Pd))
    # Get fractional-x of each Pd atom in the same sorted order.
    # We need atoms IN sort=True order to compute x-based sign.
    from ase.io import write as _w, read as _r
    tmp = OUT/'.tmp_sorted_S3clean.vasp'
    _w(str(tmp), atoms, format='vasp', sort=True, direct=True, vasp5=True)
    s_atoms = _r(tmp)
    fx = s_atoms.get_scaled_positions()[:,0]
    x_med = np.median(fx[pd_indices_after_sort])
    tmp.unlink()
    magmom_parts = []
    # O group: nonmagnetic
    magmom_parts.append(f'{n_O}*0.000')
    # Pd group: sign by fractional x
    for i in pd_indices_after_sort:
        m = 0.5 if fx[i] < x_med else -0.5
        magmom_parts.append(f'{m:+.3f}')
    return ' '.join(magmom_parts)

# ============================================================================
# Build seed dirs
# ============================================================================
manifest = []

# --- S1 clean seeds ---
S1_clean = read(G2/SDIRS['S1']/'CONTCAR')
S1_species = ['Pd']; S1_counts = [80]

for name, ispin, mom_pd in [('S1_clean__ispin1', 1, None),
                            ('S1_clean__lowmag',  2, 0.2),
                            ('S1_clean__highmag', 2, 1.0)]:
    if ispin == 1:
        sb = _spin_block(1)
    else:
        magmom = f'{S1_counts[0]}*{mom_pd:.3f}'
        sb = _spin_block(2, magmom=magmom)
    d = build_seed(S1_clean, name, ismear=1, sigma=0.10, has_dipole=True, spin_block=sb)
    manifest.append({'name':name,'source':'G2 S1 clean','ispin':ispin,
                     'magmom':(magmom if ispin==2 else 'default'),'note':'even electron'})
    print(f'  {name}')

# --- S1 CO idx=3 seeds — take T1.16 CONTCAR ---
S1CO_contcar = ROOT/'calculations/T1_16_DFT_L2/S1/CO/01_CO_idx00003/CONTCAR'
S1CO_atoms = read(S1CO_contcar)
# after sort=True species order = C, O, Pd; counts 1, 1, 80
for name, ispin, mom_pd, mom_C, mom_O in [
        ('S1_CO_idx3__ispin1',  1, None, None, None),
        ('S1_CO_idx3__lowmag',  2, 0.2,  0.0,  0.0),
        ('S1_CO_idx3__highmag', 2, 1.0,  0.0,  0.0)]:
    if ispin == 1:
        sb = _spin_block(1)
    else:
        magmom = f'1*{mom_C:.3f} 1*{mom_O:.3f} 80*{mom_pd:.3f}'
        sb = _spin_block(2, magmom=magmom)
    d = build_seed(S1CO_atoms, name, ismear=1, sigma=0.10, has_dipole=True, spin_block=sb)
    manifest.append({'name':name,'source':f'{S1CO_contcar.name}','ispin':ispin,
                     'magmom':(magmom if ispin==2 else 'default'),'note':'even electron'})
    print(f'  {name}')

# --- S3 clean seeds ---
S3_clean = read(G2/SDIRS['S3']/'CONTCAR')
n_O_S3 = sum(1 for s in S3_clean.get_chemical_symbols() if s=='O')
n_Pd_S3 = sum(1 for s in S3_clean.get_chemical_symbols() if s=='Pd')

# nonmagnetic (ISPIN=1)
d = build_seed(S3_clean, 'S3_clean__nonmag', ismear=0, sigma=0.05, has_dipole=True,
               spin_block=_spin_block(1))
manifest.append({'name':'S3_clean__nonmag','source':'G2 S3 clean','ispin':1,'magmom':'default'})
print('  S3_clean__nonmag')

# FM seed
magmom_fm = f'{n_O_S3}*0.000 {n_Pd_S3}*0.500'
d = build_seed(S3_clean, 'S3_clean__fm', ismear=0, sigma=0.05, has_dipole=True,
               spin_block=_spin_block(2, magmom=magmom_fm))
manifest.append({'name':'S3_clean__fm','source':'G2 S3 clean','ispin':2,'magmom':magmom_fm})
print('  S3_clean__fm')

# AFM seed (Pd atoms split by fractional-x)
magmom_afm = afm_magmom_S3clean(S3_clean)
d = build_seed(S3_clean, 'S3_clean__afm', ismear=0, sigma=0.05, has_dipole=True,
               spin_block=_spin_block(2, magmom=magmom_afm))
manifest.append({'name':'S3_clean__afm','source':'G2 S3 clean','ispin':2,'magmom':magmom_afm})
print('  S3_clean__afm')

# --- S3 CH3O idx=315 seeds — odd electron ---
S3CH3O_contcar = ROOT/'calculations/T1_16_DFT_L2/S3/CH3O/00_CH3O_idx00315/CONTCAR'
S3CH3O_atoms = read(S3CH3O_contcar)
# after sort=True: C H O Pd = 1, 3, N_O_all, 64
# adsorbate O = C-bonded O (one). Substrate O = rest.
# Species order in POSCAR: C, H, O, Pd (counts 1, 3, N_O_all, 64).
# In MAGMOM, we need to know which O within the O-block is the ads O.
# Simplest: unpaired μ on the first O entry in the O-block that is C-bonded.
# ASE sort=True preserves within-species order → substrate O (that was in G2) comes
# first, then adsorbate O appended last. So ads O is the LAST O in the block.
n_ads_O = 1
n_sub_O = sum(1 for s in S3CH3O_atoms.get_chemical_symbols() if s=='O') - n_ads_O
n_Pd_ads = 64

def build_ch3o_magmom(pd_mom, sub_O_mom, ads_O_mom, C_mom=0.0, H_mom=0.0):
    # species order after sort=True: C(1) H(3) O(N_O_all) Pd(64)
    # ads O is the LAST in the O-group.
    return (f'1*{C_mom:.3f} 3*{H_mom:.3f} '
            f'{n_sub_O}*{sub_O_mom:.3f} 1*{ads_O_mom:.3f} '
            f'{n_Pd_ads}*{pd_mom:.3f}')

# doublet seed
magmom_dbl = build_ch3o_magmom(pd_mom=0.0, sub_O_mom=0.0, ads_O_mom=1.0)
d = build_seed(S3CH3O_atoms, 'S3_CH3O_idx315__doublet',
               ismear=0, sigma=0.05, has_dipole=True,
               spin_block=_spin_block(2, magmom=magmom_dbl, nupdown=1))
manifest.append({'name':'S3_CH3O_idx315__doublet','source':f'{S3CH3O_contcar.name}',
                 'ispin':2,'nupdown':1,'magmom':magmom_dbl,'note':'odd electron doublet'})
print('  S3_CH3O_idx315__doublet')

# unconstrained seed
magmom_unc = build_ch3o_magmom(pd_mom=0.3, sub_O_mom=0.0, ads_O_mom=1.0)
d = build_seed(S3CH3O_atoms, 'S3_CH3O_idx315__unconstrained',
               ismear=0, sigma=0.05, has_dipole=True,
               spin_block=_spin_block(2, magmom=magmom_unc))
manifest.append({'name':'S3_CH3O_idx315__unconstrained','source':f'{S3CH3O_contcar.name}',
                 'ispin':2,'magmom':magmom_unc,'note':'odd electron, no NUPDOWN'})
print('  S3_CH3O_idx315__unconstrained')

# manifest
import csv
keys = ['name','source','ispin','nupdown','magmom','note']
with open(OUT/'manifest.csv','w',newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
    for r in manifest:
        r.setdefault('nupdown',''); r.setdefault('note','')
        w.writerow(r)

readme = f"""# Magnetic seed test — resolve MAGMOM initialization uncertainty

Reviewer 2026-07-17 round 4:
Current T1.16 candidates converged to significant magnetization (~21 μ_B for
S1/CO, ~11 μ_B for S3/CH3O) using VASP default MAGMOM=1.0/atom initialization.
Layer-wise decomposition shows moments distributed throughout the slab, not
surface-localized — suggests **metastable FM state from default seed**, not
physical surface Stoner magnetism.

This test computes STATIC single-point energies with different MAGMOM seeds on
the SAME geometry. Compare final `energy(sigma→0)` and final total magnetization.

## Seeds (11 dirs total)

| dir | ISPIN | NUPDOWN | MAGMOM | description |
|---|---|---|---|---|
| S1_clean__ispin1  | 1 | – | (none) | closed-shell, no spin polarization |
| S1_clean__lowmag  | 2 | – | 80*0.2 | small FM seed |
| S1_clean__highmag | 2 | – | 80*1.0 | default VASP seed reproduction |
| S1_CO_idx3__ispin1  | 1 | – | – | ditto with CO adsorbed |
| S1_CO_idx3__lowmag  | 2 | – | 1*0 1*0 80*0.2 | |
| S1_CO_idx3__highmag | 2 | – | 1*0 1*0 80*1.0 | |
| S3_clean__nonmag | 1 | – | – | S3 clean nonmagnetic |
| S3_clean__fm     | 2 | – | 64*0 64*0.5 | ferromagnetic seed on Pd |
| S3_clean__afm    | 2 | – | Pd atoms split ± by fractional-x | antiferromagnetic seed |
| S3_CH3O_idx315__doublet       | 2 | 1 | substrate 0, ads O = 1 | doublet ground state |
| S3_CH3O_idx315__unconstrained | 2 | – | Pd 0.3, ads O 1 | unconstrained radical |

## Common

- IBRION = -1, NSW = 0  (STATIC — no ionic movement)
- ISTART = 0, ICHARG = 2  (no WAVECAR/CHGCAR reuse — seed independence)
- KSPACING = 0.25 (matches T1.17)
- LDIPOL/IDIPOL = 3 (matches T1.17 slab)
- ENCUT/PREC/LASPH/IVDW/EDIFF/EDIFFG all match T1.17
- Same POSCAR geometry as the parent T1.16 CONTCAR (adsorbate) or G2 CONTCAR (clean)

## After submission

For each dir, extract from OUTCAR:
- `energy(sigma->0)` (grep last)
- Total magnetization (`grep "number of electron" OUTCAR | tail -1`)
- Per-atom moment (mag(x) table)
- Layer-wise Pd moment sum
- Electronic convergence status

## Decision rules

- All seeds converge to **same E and moment** → magnetic solution is stable → adopt.
- Seeds converge to **different moments but one has clearly lower E** → adopt the lowest-E state.
- Multiple states within **5–10 meV** → magnetic uncertainty; run NUPDOWN fixed-spin scan.
- Multiple magnetization basins with significant E span → propagate as an uncertainty in downstream ΔG.

Once T1.17 MAGMOM policy is decided from this test, regenerate T1.17 dirs with
that policy before mass submission.
"""
(OUT/'README.md').write_text(readme)
print(f'\n{len(manifest)} seed dirs created under {OUT.relative_to(ROOT)}')
