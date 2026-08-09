"""T1.17 setup — VASPsol Level-2 solvation with COMPLETE-LAYER mask.

Post-review 2026-07-17 round 7 (final decision):

The T1.16 production mask (bottom-half-of-substrate by median z) is
PARTIAL-LAYER: it splits e.g. S1 layer 3 into 8 fixed + 8 free within a
single atomic plane. This breaks in-plane symmetry — atoms with the same
z-value get different constraints because their z differs by <0.01 Å.

The G2 original mask is COMPLETE-LAYER by construction. Verified via
z-clustering: every atomic plane is either fully fixed or fully free
across all 5 surfaces.

  surface   G2 (complete-layer)   T1.16 (partial-layer)   T1.17 uses
  S1        32  (L1+L2)           40  (L1+L2 + half of L3)   32
  S2        40  (L1+L2)           56  (L1+L2 + partial L3)   40
  S3        32  (L1+L2)           64  (L1..L4)               32
  S3b       32  (L1+L2)           52  (L1..L3 + partial)     32
  S4        42  (L1..L4)          72  (L1..L4 + partial)     42

T1.17 uses G2's complete-layer mask.
Ads-slab dirs: G2 fixed indices mapped to T1.16 CONTCAR by species +
position (Euclidean, tol 1.5 Å); those atom positions are overwritten
with G2's canonical coordinates and marked fixed. Free substrate atoms
(including the 8 layer-3 atoms that T1.16 held fixed) keep T1.16-relaxed
positions and re-relax freely under VASPsol.
Clean slab: G2 CONTCAR verbatim (32/40/32/32/42 fixed).
Result: ads and clean dirs share identical fixed-atom count + identical
fixed-atom coordinates + identical cell + identical POTCAR substrate order.

T1.16 remains valid as SCREENING (relative ranking within a surface, all
candidates use same mask). Not used as final absolute-energy reference.

MAGMOM: left at VASP default — magnetic_seed_test/ resolves policy before
final T1.17 mass submission.
"""
import json, csv
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
OUT = ROOT/'calculations/T1_17_VASPsol'
OUT.mkdir(exist_ok=True)
G2 = ROOT/'calculations/G2_slab'
POT_LIB = Path('/home/hyunjin/POTENTIAL/potpaw_PBE')
POT_FOLDER = {'C':'C','H':'H','O':'O','Pd':'Pd_pv'}
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

INCAR_METAL = """SYSTEM = pddmc T1.17 VASPsol (methanol, complete-layer mask)
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
ISTART = 0
ICHARG = 2
LSOL = .TRUE.
EB_K = 32.6
TAU = 0
"""
INCAR_OXIDE = INCAR_METAL.replace("ISMEAR = 1\nSIGMA = 0.1",
                                   "ISMEAR = 0\nSIGMA = 0.05")

# H200 submit script — replaced by scripts/apply_h200_submit_and_fixes.py
# For fresh setup runs, we use the same H200 template inline.
SUBMIT_SCRIPT = """#!/bin/bash
#SBATCH --partition=h200q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err

CONTAINER=${CONTAINER:-/scratch/taehun1/hyunjin/vasp_vaspsol.sif}
VASP_BIN=${VASP_BIN:-vasp_std}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | Container: ${CONTAINER}"
echo "Start: $(date)"
srun --mpi=pmix singularity exec --nv "${CONTAINER}" "${VASP_BIN}"
echo "End: $(date)"
"""

# ============================================================================
# G2 canonical mask (complete-layer, verified via z-clustering)
# ============================================================================
def load_G2_mask(sid):
    g2 = read(G2/SDIRS[sid]/'CONTCAR')
    fx = next((sorted(int(i) for i in c.index) for c in g2.constraints
               if isinstance(c, FixAtoms)), [])
    if not fx:
        raise ValueError(f'{sid}: G2 CONTCAR has no FixAtoms constraint')
    return g2, fx

def apply_G2_mask_to_ads(ads_atoms, g2_atoms, g2_fixed_indices, tol=1.5):
    """Map each G2-fixed atom to an ads-slab atom by species+position match,
    overwrite the ads position with the canonical G2 coordinate, mark fixed.

    Free substrate atoms (including layer-3 atoms that T1.16 held fixed) and
    adsorbate atoms are untouched — they re-relax under VASPsol.
    """
    ads_new = ads_atoms.copy()
    ads_new.set_constraint()   # wipe existing T1.16 FixAtoms
    ads_new_fixed_idx = []
    for gi in g2_fixed_indices:
        g_sym = g2_atoms.get_chemical_symbols()[gi]
        g_pos = g2_atoms.positions[gi]
        best_ai, best_d = None, float('inf')
        for ai in range(len(ads_new)):
            if ai in ads_new_fixed_idx: continue
            if ads_new.get_chemical_symbols()[ai] != g_sym: continue
            d = float(np.linalg.norm(ads_new.positions[ai] - g_pos))
            if d < best_d: best_d, best_ai = d, ai
        if best_ai is None or best_d > tol:
            raise ValueError(f'Cannot map G2 atom {gi} ({g_sym}) to ads '
                             f'within {tol} Å (closest: {best_d} Å)')
        ads_new.positions[best_ai] = g_pos
        ads_new_fixed_idx.append(best_ai)
    ads_new.set_constraint(FixAtoms(indices=ads_new_fixed_idx))
    return ads_new, ads_new_fixed_idx

def write_common(dest, atoms, sid, extra_meta):
    write(str(dest/'POSCAR'), atoms, format='vasp', direct=True, sort=True, vasp5=True)
    (dest/'INCAR').write_text(INCAR_METAL if sid=='S1' else INCAR_OXIDE)
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    (dest/'POTCAR').write_text(''.join(
        (POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species))
    (dest/'submit_vasp.sh').write_text(SUBMIT_SCRIPT)
    (dest/'metadata.json').write_text(json.dumps({
        'level': 'T1.17 VASPsol (methanol, LSOL=T, EB_K=32.6, TAU=0)',
        'restart_kind': 'geometric only (ISTART=0, ICHARG=2)',
        'mask_convention': ('G2 complete-layer mask (32/40/32/32/42 per surface). '
                            'Applied by species+position mapping from G2 CONTCAR. '
                            'This is a COMPLETE-LAYER constraint — every atomic '
                            'plane is either fully fixed or fully free. Differs '
                            'from T1.16 which used a partial-layer bottom-half rule.'),
        'poscar_species': species,
        **extra_meta,
    }, indent=2))

def build_adsorbate_dir(sid, ads, idx, l1_contcar):
    dest = OUT/sid/ads/f'{ads}_idx{idx:05d}'
    if dest.exists(): return None, 'exists'
    dest.mkdir(parents=True)
    g2, g2_fixed = load_G2_mask(sid)
    ads_atoms = read(l1_contcar)
    ads_new, fixed_idx = apply_G2_mask_to_ads(ads_atoms, g2, g2_fixed)
    write_common(dest, ads_new, sid, extra_meta={
        'sid':sid,'ads':ads,'idx':int(idx),
        'source_L1_CONTCAR': str(l1_contcar.relative_to(ROOT)),
        'canonical_from_G2': str((G2/SDIRS[sid]/'CONTCAR').relative_to(ROOT)),
        'n_atoms_total': len(ads_new),
        'n_atoms_fixed': len(fixed_idx),
        'T1_16_partial_layer_atoms_now_freed': None,   # count elsewhere if needed
        'purpose': ('Level-2 solvated adsorbate slab with G2 complete-layer mask. '
                    'Layer-3 atoms that T1.16 fixed via partial-layer rule are '
                    'now FREE and re-relax under VASPsol.'),
    })
    return dest, 'created'

def build_clean_slab_dir(sid):
    """Clean slab = G2 CONTCAR verbatim (complete-layer mask + positions)."""
    dest = OUT/f'{sid}_clean'
    if dest.exists(): return None, 'exists'
    dest.mkdir(parents=True)
    g2 = read(G2/SDIRS[sid]/'CONTCAR')
    n_fixed = sum(len(c.index) for c in g2.constraints if isinstance(c, FixAtoms))
    write_common(dest, g2, sid, extra_meta={
        'sid':sid,'kind':'clean_slab',
        'source_CONTCAR': str((G2/SDIRS[sid]/'CONTCAR').relative_to(ROOT)),
        'n_atoms_total': len(g2),
        'n_atoms_fixed': n_fixed,
        'purpose': 'Clean-slab reference with same G2 complete-layer mask as ads.',
    })
    return dest, 'created'

def pick_top1_per_group():
    picks = {}
    for r in csv.DictReader(open(ROOT/'paper_data/07_dft_results.csv')):
        if r['DFT_status'] != 'DONE': continue
        key = (r['sid'], r['ads'])
        E = float(r['E_DFT_sigma0_eV'])
        if key not in picks or E < picks[key]['E']:
            picks[key] = {'idx': int(r['idx']), 'E': E,
                          'contcar_path': ROOT/r['contcar_path']}
    return picks

def auto_consistency_check():
    issues = []
    ads_dirs = list(OUT.glob('*/*/[Cc]*_idx*'))
    for ad in ads_dirs:
        sid = ad.parts[-3]
        clean = OUT/f'{sid}_clean'
        if not clean.exists():
            issues.append(f'{ad.relative_to(OUT)}: no {sid}_clean counterpart'); continue
        a_ads = read(ad/'POSCAR'); a_cl = read(clean/'POSCAR')
        if not np.allclose(a_ads.cell.array, a_cl.cell.array, atol=1e-6):
            issues.append(f'{sid}: cell differs')
        f_ads = sorted([int(i) for c in a_ads.constraints if isinstance(c,FixAtoms) for i in c.index])
        f_cl  = sorted([int(i) for c in a_cl.constraints if isinstance(c,FixAtoms) for i in c.index])
        if len(f_ads) != len(f_cl):
            issues.append(f'{sid}: mask count differs (ads {len(f_ads)} vs clean {len(f_cl)})')
        f_ads_pos = sorted(map(tuple, np.round(a_ads.positions[f_ads], 6).tolist()))
        f_cl_pos  = sorted(map(tuple, np.round(a_cl.positions[f_cl], 6).tolist()))
        if f_ads_pos != f_cl_pos:
            issues.append(f'{sid}: fixed-atom coordinates differ')
    return issues

def main():
    manifest = []; picks = pick_top1_per_group()
    created = skipped = 0
    for (sid, ads), p in sorted(picks.items()):
        try:
            dest, status = build_adsorbate_dir(sid, ads, p['idx'], p['contcar_path'])
        except Exception as e:
            print(f'  ERROR {sid}/{ads}: {e}')
            manifest.append({'kind':'ads','sid':sid,'ads':ads,'idx':p['idx'],
                             'E_L1':p['E'],'status':f'ERROR: {e}','dir':''})
            continue
        if dest is None: skipped += 1
        else:
            created += 1
            print(f'  ads   {sid}/{ads} idx={p["idx"]} → {dest.relative_to(ROOT)}')
        manifest.append({'kind':'ads','sid':sid,'ads':ads,'idx':p['idx'],
                         'E_L1':p['E'],'status':'created' if dest else 'exists',
                         'dir':str((OUT/sid/ads/f'{ads}_idx{p["idx"]:05d}').relative_to(ROOT))})
    for sid in ['S1','S2','S3','S3b','S4']:
        dest, status = build_clean_slab_dir(sid)
        if dest is None: skipped += 1
        else:
            created += 1
            print(f'  clean {sid} → {dest.relative_to(ROOT)}')
        manifest.append({'kind':'clean','sid':sid,'ads':'','idx':'',
                         'E_L1':'','status':'created' if dest else 'exists',
                         'dir':str((OUT/f'{sid}_clean').relative_to(ROOT))})
    # pending ads groups
    all_groups = [(s,a) for s in ['S1','S2','S3','S3b','S4']
                       for a in ['CO','CH3O','coads']
                       if not (s=='S4' and a=='coads')]
    for g in all_groups:
        if g not in picks:
            manifest.append({'kind':'ads','sid':g[0],'ads':g[1],'idx':'',
                             'E_L1':'','status':'PENDING L1','dir':''})
    with open(OUT/'manifest.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['kind','sid','ads','idx','E_L1','status','dir'])
        w.writeheader()
        for r in manifest: w.writerow(r)
    print(f'\nT1_17_VASPsol: created={created}, skipped={skipped}, '
          f'PENDING={sum(1 for r in manifest if "PENDING" in r["status"])}')
    issues = auto_consistency_check()
    if issues:
        print('\nCONSISTENCY ERRORS:')
        for e in issues: print(f'  {e}')
    else:
        print('Auto-consistency: all ads/clean pairs share IDENTICAL fixed count + '
              'coordinates + cell.')

if __name__=='__main__':
    main()
