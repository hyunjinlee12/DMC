"""Mask sensitivity test — 40 fixed (T1.16 partial) vs 32 fixed (G2 complete).

Reviewer 2026-07-17 round 7:
"S1 대표 구조에서 `40 fixed`와 `32 fixed` 비교"

Compares full vacuum re-relaxation of S1/CO idx=3 (representative top-1) and
S1_clean under the two mask choices. All INCAR settings are identical to
T1.16 (vacuum, PBE-D3, ISPIN=2, default MAGMOM, NSW=200, IBRION=2, KSPACING=0.25);
only the FixAtoms constraint differs.

Result → ΔE_ads(mask=40) vs ΔE_ads(mask=32).
Acceptance: agreement within ~5-10 meV means partial-vs-complete mask choice
is not material for this system; divergence means we must use complete-layer.

4 dirs total, ~1 day of walltime.

No jobs submitted. Static-only INCAR could also be used, but only relaxation
reveals the mask's true effect (static at same geometry gives ~0 diff by
construction, since the fixed atoms don't move regardless).
"""
import json
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
OUT = ROOT/'calculations/mask_sensitivity_test'
OUT.mkdir(exist_ok=True)
G2 = ROOT/'calculations/G2_slab'
POT_LIB = Path('/home/hyunjin/POTENTIAL/potpaw_PBE')
POT_FOLDER = {'C':'C','H':'H','O':'O','Pd':'Pd_pv'}

INCAR = """SYSTEM = mask sensitivity test (vacuum, PBE-D3)
ENCUT = 520
PREC = Accurate
LASPH = .TRUE.
ADDGRID = .TRUE.
ISPIN = 2
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
IBRION = 2
NSW = 200
ISIF = 2
ISMEAR = 1
SIGMA = 0.1
EDIFFG = -0.03
LDIPOL = .TRUE.
IDIPOL = 3
KSPACING = 0.25
"""

SUBMIT = """#!/bin/bash
#SBATCH --partition=h200q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err
CONTAINER=${CONTAINER:-/scratch/taehun1/hyunjin/vasp_vaspsol.sif}
VASP_BIN=${VASP_BIN:-vasp_std}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
srun --mpi=pmix singularity exec --nv "${CONTAINER}" "${VASP_BIN}"
"""

def apply_mask_by_g2_indices(atoms_target, atoms_g2, g2_fixed_indices, tol=1.5):
    """Map G2 fixed indices to atoms_target by species+position; return new atoms
    with those indices fixed and their positions overwritten with G2 canonical."""
    a = atoms_target.copy()
    a.set_constraint()
    picked = []
    for gi in g2_fixed_indices:
        g_sym = atoms_g2.get_chemical_symbols()[gi]
        g_pos = atoms_g2.positions[gi]
        best_ai, best_d = None, float('inf')
        for ai in range(len(a)):
            if ai in picked: continue
            if a.get_chemical_symbols()[ai] != g_sym: continue
            d = float(np.linalg.norm(a.positions[ai] - g_pos))
            if d < best_d: best_d, best_ai = d, ai
        if best_ai is None or best_d > tol:
            raise ValueError(f'cannot map G2 atom {gi}')
        a.positions[best_ai] = g_pos
        picked.append(best_ai)
    a.set_constraint(FixAtoms(indices=picked))
    return a, picked

def apply_bottom_half_median_mask(atoms):
    """Reproduce T1.16 partial-layer rule: bottom half of substrate by median z."""
    from ase.constraints import FixAtoms
    a = atoms.copy()
    syms = a.get_chemical_symbols()
    # substrate = non-adsorbate (C, H, C-bonded O removed)
    c_i = [i for i,s in enumerate(syms) if s=='C']
    h_i = [i for i,s in enumerate(syms) if s=='H']
    o_i = [i for i,s in enumerate(syms) if s=='O']
    ads_o = []
    for oi in o_i:
        for ci in c_i:
            if a.get_distance(ci, oi, mic=True) < 1.5:
                ads_o.append(oi); break
    ads_set = set(c_i) | set(h_i) | set(ads_o)
    sub_indices = [i for i in range(len(a)) if i not in ads_set]
    z_sub = a.positions[sub_indices, 2]
    z_med = float(np.median(z_sub))
    fixed = [i for i in sub_indices if a.positions[i,2] < z_med]
    a.set_constraint()
    a.set_constraint(FixAtoms(indices=fixed))
    return a, fixed

def build_dir(name, atoms, meta):
    dest = OUT/name
    dest.mkdir(exist_ok=True)
    write(str(dest/'POSCAR'), atoms, format='vasp', direct=True, sort=True, vasp5=True)
    (dest/'INCAR').write_text(INCAR)
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    (dest/'POTCAR').write_text(''.join(
        (POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species))
    (dest/'submit_vasp.sh').write_text(SUBMIT)
    (dest/'metadata.json').write_text(json.dumps(meta, indent=2))
    return dest

# ============================================================================
# Build 4 dirs
# ============================================================================
S1_g2 = read(G2/'S1_Pd100/CONTCAR')
S1_g2_fixed = next((sorted(int(i) for i in c.index) for c in S1_g2.constraints
                    if isinstance(c, FixAtoms)), [])

S1_CO_l1 = read(ROOT/'calculations/T1_16_DFT_L2/S1/CO/01_CO_idx00003/CONTCAR')

# 1. S1_CO_idx3_40fix — reproduce T1.16 partial-layer mask
ads_40, fixed_40 = apply_bottom_half_median_mask(S1_CO_l1)
build_dir('S1_CO_idx3_40fix_partial', ads_40, {
    'name':'S1_CO_idx3_40fix_partial',
    'source':'T1_16_DFT_L2/S1/CO/01_CO_idx00003/CONTCAR',
    'mask':'T1.16 partial-layer (bottom half of substrate by median z)',
    'n_atoms':len(ads_40),'n_fixed':len(fixed_40),
    'note':'40 fixed; layer 3 split 8/8. Reproduces current T1.16 setting.'})
print(f'S1_CO_idx3_40fix_partial: {len(fixed_40)} fixed')

# 2. S1_CO_idx3_32fix — G2 complete-layer mask
ads_32, fixed_32 = apply_mask_by_g2_indices(S1_CO_l1, S1_g2, S1_g2_fixed)
build_dir('S1_CO_idx3_32fix_complete', ads_32, {
    'name':'S1_CO_idx3_32fix_complete',
    'source':'T1_16_DFT_L2/S1/CO/01_CO_idx00003/CONTCAR',
    'mask':'G2 complete-layer (L1+L2)',
    'n_atoms':len(ads_32),'n_fixed':len(fixed_32),
    'note':'32 fixed; layers 1+2 fully fixed, layer 3 fully free.'})
print(f'S1_CO_idx3_32fix_complete: {len(fixed_32)} fixed')

# 3. S1_clean_40fix — synthesize partial mask on G2 clean slab
def apply_partial_on_g2(g2_atoms):
    a = g2_atoms.copy()
    z = a.positions[:, 2]
    z_med = float(np.median(z))
    fixed = [i for i in range(len(a)) if z[i] < z_med]
    a.set_constraint()
    a.set_constraint(FixAtoms(indices=fixed))
    return a, fixed

clean_40, fixed_c40 = apply_partial_on_g2(S1_g2)
build_dir('S1_clean_40fix_partial', clean_40, {
    'name':'S1_clean_40fix_partial',
    'source':'G2_slab/S1_Pd100/CONTCAR',
    'mask':'partial-layer (bottom-half by median z applied to G2 clean slab)',
    'n_atoms':len(clean_40),'n_fixed':len(fixed_c40),
    'note':'Reproduces the T1.16-style partial mask on the clean slab. Only for '
           'ΔE_ads(40fix) computation.'})
print(f'S1_clean_40fix_partial: {len(fixed_c40)} fixed')

# 4. S1_clean_32fix — G2 CONTCAR verbatim
build_dir('S1_clean_32fix_complete', S1_g2, {
    'name':'S1_clean_32fix_complete',
    'source':'G2_slab/S1_Pd100/CONTCAR',
    'mask':'G2 original complete-layer (32 = L1+L2)',
    'n_atoms':len(S1_g2),'n_fixed':len(S1_g2_fixed),
    'note':'G2 CONTCAR verbatim, positions + constraint copied.'})
print(f'S1_clean_32fix_complete: {len(S1_g2_fixed)} fixed')

# README
(OUT/'README.md').write_text("""# Mask sensitivity test — S1 only

Reviewer 2026-07-17 round 7: compare partial-layer (T1.16) vs complete-layer
(G2) mask on S1 to quantify the effect on ΔE_ads.

## 4 dirs
| dir | mask | n_fixed | source |
|---|---|---|---|
| S1_CO_idx3_40fix_partial   | T1.16 partial (median z) | 40 | T1.16 top-1 CONTCAR |
| S1_CO_idx3_32fix_complete  | G2 complete-layer (L1+L2) | 32 | T1.16 top-1 CONTCAR |
| S1_clean_40fix_partial     | partial (median z) on G2  | 40 | G2 CONTCAR |
| S1_clean_32fix_complete    | G2 complete-layer         | 32 | G2 CONTCAR verbatim |

All 4: vacuum PBE-D3, IBRION=2 NSW=200, KSPACING=0.25, ISPIN=2 (default MAGMOM),
INCAR otherwise identical to T1.16.

## Analysis
After all 4 finish:
```
ΔE_ads(40) = E(S1_CO_40fix)    - E(S1_clean_40fix)    - μ_CO(gas)
ΔE_ads(32) = E(S1_CO_32fix)    - E(S1_clean_32fix)    - μ_CO(gas)
Δ           = ΔE_ads(40) − ΔE_ads(32)
```

- |Δ| < 10 meV → mask choice is immaterial for downstream (use either)
- |Δ| ≥ 10 meV → complete-layer mask (32) is the correct choice
""")
print(f'\n4 dirs in {OUT.relative_to(ROOT)}/')
