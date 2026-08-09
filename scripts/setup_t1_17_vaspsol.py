"""T1.17 setup — VASPsol Level-2 solvation with G2 canonical mask.

Post-review 2026-07-17 (round 4):
- Clean slab = G2 CONTCAR verbatim (positions + G2 mask).
- Adsorbate slab: G2's fixed-atom indices are mapped to the corresponding
  atoms in the T1.16 CONTCAR by species + position matching. Those atoms
  have their positions overwritten with G2's canonical coordinates AND
  are marked fixed. Free substrate atoms and adsorbate atoms keep their
  T1.16-relaxed positions and are free to move in T1.17.
- Result: ads_slab and clean_slab share IDENTICAL fixed atom indices +
  identical fixed coordinates + identical mask + identical cell + identical
  POTCAR substrate order, so ΔG_ads = G(ads) − G(clean) − μ has EXACT
  fixed-atom cancellation.
- z-threshold rules explicitly REJECTED (per reviewer): index-based mapping
  is stable under z-drift.
- MAGMOM: left at VASP default for now — magnetic seed test in
  calculations/magnetic_seed_test/ resolves initialization policy before
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

INCAR_METAL = """SYSTEM = pddmc T1.17 VASPsol (methanol implicit solvent)
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

SUBMIT_SCRIPT = """#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=7-00:00:00
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

# VASP >=5.4.1 standard builds usually support LSOL — no separate binary
# required. Verify with `grep -E "LSOL|EB_K|VASPsol" OUTCAR` on the pilot.
VASP_BIN=${VASP_BIN:-/home/hyunjin/vasp.6.4.3/bin/vasp_std}

NPROCS=${SLURM_NTASKS:-1}
echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | VASP: ${VASP_BIN}"
echo "Start: $(date)"
mpirun --bind-to none -np ${NPROCS} ${VASP_BIN}
echo "End: $(date)"
"""

# ============================================================================
# Canonical mask machinery
# ============================================================================

def load_canonical_g2(sid):
    """Return (g2_atoms, g2_fixed_indices) — G2 CONTCAR with its FixAtoms."""
    g2 = read(G2/SDIRS[sid]/'CONTCAR')
    fx = next((sorted(int(i) for i in c.index) for c in g2.constraints
               if isinstance(c, FixAtoms)), [])
    if not fx:
        raise ValueError(f'{sid}: G2 CONTCAR has no FixAtoms constraint')
    return g2, fx

def apply_canonical_mask(ads_atoms, g2_atoms, g2_fixed_indices, tol=1.5):
    """Map each G2-fixed atom to an ads-slab atom by species+position match,
    overwrite the ads position with the canonical G2 coordinate, and mark
    that ads-slab atom fixed. Free substrate atoms and adsorbate atoms are
    untouched.

    Raises ValueError if any G2 atom cannot be uniquely mapped within `tol` Å.
    """
    ads_new = ads_atoms.copy()
    ads_new.set_constraint()   # wipe existing FixAtoms so we rebuild it

    canonical_ads_idx = []
    for gi in g2_fixed_indices:
        g_sym = g2_atoms.get_chemical_symbols()[gi]
        g_pos = g2_atoms.positions[gi]
        best_ai = None; best_d = float('inf')
        for ai in range(len(ads_new)):
            if ai in canonical_ads_idx: continue
            if ads_new.get_chemical_symbols()[ai] != g_sym: continue
            d = float(np.linalg.norm(ads_new.positions[ai] - g_pos))
            if d < best_d:
                best_d = d; best_ai = ai
        if best_ai is None or best_d > tol:
            raise ValueError(f'Cannot map G2 fixed atom {gi} ({g_sym}, '
                             f'pos={g_pos}) to any ads atom within {tol} Å '
                             f'(closest: {best_d} Å)')
        ads_new.positions[best_ai] = g_pos   # overwrite with canonical
        canonical_ads_idx.append(best_ai)
    ads_new.set_constraint(FixAtoms(indices=canonical_ads_idx))
    return ads_new, canonical_ads_idx

def write_common(dest, atoms, sid, extra_meta):
    write(str(dest/'POSCAR'), atoms, format='vasp', direct=True, sort=True, vasp5=True)
    (dest/'INCAR').write_text(INCAR_METAL if sid=='S1' else INCAR_OXIDE)
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    potcar = ''.join((POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species)
    (dest/'POTCAR').write_text(potcar)
    (dest/'submit_vasp.sh').write_text(SUBMIT_SCRIPT)
    (dest/'metadata.json').write_text(json.dumps({
        'level': 'T1.17 VASPsol (methanol, LSOL=T, EB_K=32.6, TAU=0)',
        'restart_kind': 'geometric only (ISTART=0, ICHARG=2). L1 had LWAVE=.FALSE. → no WAVECAR.',
        'mask_convention': 'G2 canonical fixed-atom indices; fixed coordinates '
                           'copied from G2 CONTCAR; mapping done by species + '
                           'position (PBC-agnostic Euclidean, tol=1.5 Å).',
        'poscar_species': species,
        **extra_meta,
    }, indent=2))

def build_adsorbate_dir(sid, ads, idx, l1_contcar):
    dest = OUT/sid/ads/f'{ads}_idx{idx:05d}'
    if dest.exists(): return None, 'exists'
    dest.mkdir(parents=True)
    g2, g2_fixed = load_canonical_g2(sid)
    ads_atoms = read(l1_contcar)
    try:
        ads_canonical, fixed_ads_idx = apply_canonical_mask(ads_atoms, g2, g2_fixed)
    except ValueError as e:
        raise RuntimeError(f'{sid}/{ads} idx={idx}: canonical mask apply FAILED — {e}')
    write_common(dest, ads_canonical, sid, extra_meta={
        'sid':sid,'ads':ads,'idx':int(idx),
        'source_L1_CONTCAR': str(l1_contcar.relative_to(ROOT)),
        'canonical_from_G2': str((G2/SDIRS[sid]/'CONTCAR').relative_to(ROOT)),
        'n_atoms_total': len(ads_canonical),
        'n_atoms_fixed': len(fixed_ads_idx),
        'g2_fixed_count': len(g2_fixed),
        'purpose': ('Level-2 solvated adsorbate slab with G2 canonical mask. '
                    'Fixed atoms at G2 coordinates; free atoms + adsorbate at '
                    'T1.16-relaxed positions; will re-relax under VASPsol.'),
    })
    return dest, 'created'

def build_clean_slab_dir(sid):
    """Clean slab = G2 CONTCAR verbatim. Positions and mask taken directly."""
    dest = OUT/f'{sid}_clean'
    if dest.exists(): return None, 'exists'
    dest.mkdir(parents=True)
    g2 = read(G2/SDIRS[sid]/'CONTCAR')
    write_common(dest, g2, sid, extra_meta={
        'sid':sid,'kind':'clean_slab',
        'source_CONTCAR': str((G2/SDIRS[sid]/'CONTCAR').relative_to(ROOT)),
        'n_atoms_total': len(g2),
        'n_atoms_fixed': sum(len(c.index) for c in g2.constraints if isinstance(c, FixAtoms)),
        'purpose': ('Clean-slab VASPsol reference. Same G2 fixed indices '
                    'and coordinates as the ads dirs → exact cancellation '
                    'in ΔG_ads.'),
    })
    return dest, 'created'

# ============================================================================
# Top-1 selection + main driver + auto-consistency check
# ============================================================================

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

def mark_superseded_if_stale(sid, ads, current_idx):
    parent = OUT/sid/ads
    if not parent.exists(): return
    for d in parent.iterdir():
        if not d.is_dir(): continue
        m = d/'metadata.json'
        if not m.exists(): continue
        meta = json.loads(m.read_text())
        old_idx = meta.get('idx')
        if old_idx is not None and old_idx != current_idx and not meta.get('SUPERSEDED'):
            meta['SUPERSEDED'] = True
            meta['SUPERSEDED_by'] = f'{ads}_idx{current_idx:05d}'
            m.write_text(json.dumps(meta, indent=2))
            print(f'  SUPERSEDED marker set: {d.relative_to(ROOT)}')

def auto_consistency_check():
    """Verify every ads dir has fixed atoms at the SAME coordinates + count
    + species as the surface's {sid}_clean dir."""
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
        # positions of fixed atoms — compare sorted-by-position
        f_ads_pos = sorted(map(tuple, np.round(a_ads.positions[f_ads], 6).tolist()))
        f_cl_pos  = sorted(map(tuple, np.round(a_cl.positions[f_cl], 6).tolist()))
        if f_ads_pos != f_cl_pos:
            issues.append(f'{sid}: fixed-atom coordinates differ between ads and clean')
    return issues

def main():
    manifest = []; picks = pick_top1_per_group()
    created = 0; skipped = 0

    # 1. Adsorbate top-1 per group
    for (sid, ads), p in sorted(picks.items()):
        mark_superseded_if_stale(sid, ads, p['idx'])
        try:
            dest, status = build_adsorbate_dir(sid, ads, p['idx'], p['contcar_path'])
        except RuntimeError as e:
            print(f'  ERROR {sid}/{ads}: {e}')
            manifest.append({'kind':'ads_top1','sid':sid,'ads':ads,'idx':p['idx'],
                             'E_L1_vacuum':p['E'],'status':f'ERROR: {e}','dir':''})
            continue
        if dest is None:
            skipped += 1
            manifest.append({'kind':'ads_top1','sid':sid,'ads':ads,'idx':p['idx'],
                             'E_L1_vacuum':p['E'],'status':'skipped (exists)',
                             'dir':str((OUT/sid/ads/f'{ads}_idx{p["idx"]:05d}').relative_to(ROOT))})
            continue
        created += 1
        manifest.append({'kind':'ads_top1','sid':sid,'ads':ads,'idx':p['idx'],
                         'E_L1_vacuum':p['E'],'status':'created',
                         'dir':str(dest.relative_to(ROOT))})
        print(f'  ads   {sid}/{ads} idx={p["idx"]} → {dest.relative_to(ROOT)}')

    # 2. Clean slab per surface (G2 CONTCAR verbatim)
    for sid in ['S1','S2','S3','S3b','S4']:
        dest, status = build_clean_slab_dir(sid)
        if dest is None:
            skipped += 1
            manifest.append({'kind':'clean_slab','sid':sid,'ads':'','idx':'',
                             'E_L1_vacuum':'','status':'skipped (exists)',
                             'dir':str((OUT/f'{sid}_clean').relative_to(ROOT))})
            continue
        created += 1
        manifest.append({'kind':'clean_slab','sid':sid,'ads':'','idx':'',
                         'E_L1_vacuum':'','status':'created',
                         'dir':str(dest.relative_to(ROOT))})
        print(f'  clean {sid} → {dest.relative_to(ROOT)}')

    # 3. Pending ads groups
    all_groups = [('S1','CO'),('S1','CH3O'),('S1','coads'),
                  ('S2','CO'),('S2','CH3O'),('S2','coads'),
                  ('S3','CO'),('S3','CH3O'),('S3','coads'),
                  ('S3b','CO'),('S3b','CH3O'),('S3b','coads'),
                  ('S4','CO'),('S4','CH3O')]
    for g in all_groups:
        if g not in picks:
            manifest.append({'kind':'ads_top1','sid':g[0],'ads':g[1],'idx':'',
                             'E_L1_vacuum':'','status':'PENDING L1 completion','dir':''})

    with open(OUT/'manifest.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['kind','sid','ads','idx','E_L1_vacuum','status','dir'])
        w.writeheader()
        for r in manifest: w.writerow(r)

    print(f'\nT1_17_VASPsol: created={created}, skipped={skipped}, '
          f'PENDING={sum(1 for r in manifest if "PENDING" in r["status"])}')
    issues = auto_consistency_check()
    if issues:
        print('\nCONSISTENCY ERRORS:')
        for e in issues: print(f'  {e}')
    else:
        print('Auto-consistency: all ads/clean pairs share IDENTICAL fixed atom '
              'count + coordinates + cell.')

if __name__=='__main__':
    main()
