"""T1.17 setup — VASPsol Level-2 solvation for each (sid, ads) group top-1.

Reads paper_data/07_dft_results.csv, picks the lowest-E_DFT candidate per
(sid, ads) group among DONE entries, and builds a T1_17_VASPsol dir with:
  POSCAR                  ← from L1 CONTCAR (already relaxed in vacuum PBE-D3)
  INCAR                   ← same as L1 + LSOL block
  POTCAR                  ← rebuilt from library to match POSCAR order
  submit_vasp_sol.sh      ← same env as L1 but expects VASP_SOL_BIN (VASPsol build)
  metadata.json           ← provenance

Idempotent: re-running only adds groups that newly turned DONE. Does not
touch existing L2 dirs or T1_16_DFT_L2/. No jobs submitted.
"""
import json, shutil, csv
from pathlib import Path
from ase.io import read, write
from ase.constraints import FixAtoms
import numpy as np

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
OUT = ROOT/'calculations/T1_17_VASPsol'
OUT.mkdir(exist_ok=True)
POT_LIB = Path('/home/hyunjin/POTENTIAL/potpaw_PBE')
POT_FOLDER = {'C':'C','H':'H','O':'O','Pd':'Pd_pv'}

# ---------- INCAR: same as L1 vacuum + VASPsol block ----------
INCAR_METAL = """SYSTEM = pddmc T1.17 VASPsol
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

# --- VASPsol (methanol) ---
LSOL = .TRUE.
EB_K = 32.6
TAU = 0
LAMBDA_D_K = 3.0
"""
INCAR_OXIDE = INCAR_METAL.replace("ISMEAR = 1\nSIGMA = 0.1",
                                   "ISMEAR = 0\nSIGMA = 0.05")

# ---------- submit script ----------
# same env as L1 but VASP_BIN → VASP_SOL_BIN (VASPsol build must exist on target server)
SUBMIT_SCRIPT = """#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err

# --- Clean conda/python from environment ---
unset PYTHONPATH PYTHONHOME CONDA_PREFIX CONDA_DEFAULT_ENV
unset CONDA_SHLVL CONDA_PROMPT_MODIFIER
export LD_LIBRARY_PATH=""

# --- NVHPC compiler + MPI ---
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

# ==============================================================
# IMPORTANT: VASPsol build (vasp_std_sol) required, not vanilla vasp_std
# ==============================================================
VASP_SOL_BIN=${VASP_SOL_BIN:-/home/hyunjin/vasp.6.4.3_sol/bin/vasp_std}

NPROCS=${SLURM_NTASKS:-1}
echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | MPI ranks: ${NPROCS} | GPUs: ${SLURM_GPUS_ON_NODE:-1}"
echo "VASP binary: ${VASP_SOL_BIN}"
echo "Start: $(date)"

mpirun --bind-to none -np ${NPROCS} ${VASP_SOL_BIN}

echo "End: $(date)"
"""

def build_l2_dir(sid, ads, idx, contcar_src, E_L1_vacuum):
    """Create one L2 VASPsol dir for a group's top-1 candidate."""
    dest = OUT/sid/ads/f'{ads}_idx{idx:05d}'
    if dest.exists():
        return None, 'already exists'
    dest.mkdir(parents=True)
    a = read(contcar_src)   # already has FixAtoms constraint from L1 relaxation
    write(str(dest/'POSCAR'), a, format='vasp', direct=True, sort=True, vasp5=True)
    (dest/'INCAR').write_text(INCAR_METAL if sid=='S1' else INCAR_OXIDE)
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    potcar = ''.join((POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species)
    (dest/'POTCAR').write_text(potcar)
    (dest/'submit_vasp_sol.sh').write_text(SUBMIT_SCRIPT)
    (dest/'metadata.json').write_text(json.dumps({
        'sid':sid,'ads':ads,'idx':int(idx),
        'source_L1_dir': str(contcar_src.parent.relative_to(ROOT)),
        'source_L1_CONTCAR': str(contcar_src.relative_to(ROOT)),
        'E_L1_vacuum_eV': E_L1_vacuum,
        'poscar_species': species,
        'level': 'T1.17 VASPsol (methanol, EB_K=32.6, TAU=0)',
        'purpose': 'Level-2 solvation single-point relaxation of L1 vacuum minimum '
                   'for adsorption energy in solvated environment (T1.18 input).',
        'note': 'Restart from L1 CONTCAR — do NOT re-relax from raw POSCAR.',
    }, indent=2))
    return dest, 'created'

def pick_top1_per_group():
    """From paper_data/07_dft_results.csv, pick lowest-E_DFT candidate per (sid, ads) DONE."""
    picks = {}
    for r in csv.DictReader(open(ROOT/'paper_data/07_dft_results.csv')):
        if r['DFT_status'] != 'DONE': continue
        key = (r['sid'], r['ads'])
        E = float(r['E_DFT_sigma0_eV'])
        if key not in picks or E < picks[key]['E']:
            picks[key] = {
                'idx': int(r['idx']),
                'E': E,
                'contcar_path': ROOT/r['contcar_path'],
            }
    return picks

def main():
    picks = pick_top1_per_group()
    manifest = []
    created = 0; skipped = 0
    for (sid, ads), p in sorted(picks.items()):
        contcar_src = p['contcar_path']
        if not contcar_src.exists():
            print(f'  MISSING CONTCAR: {contcar_src}')
            continue
        dest, status = build_l2_dir(sid, ads, p['idx'], contcar_src, p['E'])
        if dest is None:
            skipped += 1
            manifest.append({'sid':sid,'ads':ads,'idx':p['idx'],'E_L1':p['E'],
                             'dir':str((OUT/sid/ads/f'{ads}_idx{p["idx"]:05d}').relative_to(ROOT)),
                             'status':'skipped (exists)'})
            continue
        created += 1
        manifest.append({'sid':sid,'ads':ads,'idx':p['idx'],'E_L1':p['E'],
                         'dir':str(dest.relative_to(ROOT)), 'status':'created'})
        print(f'  {sid}/{ads} idx={p["idx"]} → {dest.relative_to(ROOT)}')

    # Placeholder rows for groups without L1 DONE
    groups_all = [('S1','CO'),('S1','CH3O'),('S1','coads'),
                  ('S2','CO'),('S2','CH3O'),('S2','coads'),
                  ('S3','CO'),('S3','CH3O'),('S3','coads'),
                  ('S3b','CO'),('S3b','CH3O'),('S3b','coads'),
                  ('S4','CO'),('S4','CH3O')]     # S4 coads excluded
    for g in groups_all:
        if g not in picks:
            manifest.append({'sid':g[0],'ads':g[1],'idx':'','E_L1':'',
                             'dir':'','status':'PENDING L1 completion'})

    with open(OUT/'manifest.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['sid','ads','idx','E_L1','dir','status']); w.writeheader()
        for r in manifest: w.writerow(r)

    readme = f"""# T1_17_VASPsol — Level-2 solvation setup

Regenerate: `python scripts/setup_t1_17_vaspsol.py` (idempotent).

Each candidate dir contains:
- POSCAR (copied from L1 CONTCAR — already relaxed in vacuum PBE-D3)
- INCAR (L1 settings + `LSOL=.TRUE., EB_K=32.6, TAU=0, LAMBDA_D_K=3.0`)
- POTCAR (rebuilt from library, matches POSCAR species order)
- submit_vasp_sol.sh (requires **VASPsol-enabled VASP build**; set
  `VASP_SOL_BIN` env var or edit the script before submitting)
- metadata.json (provenance + L1 vacuum energy)

## Current bundle status (auto-populated on rerun)

- L1-DONE groups get L2 dirs immediately.
- L1-PENDING groups are noted in `manifest.csv` with status=PENDING and no dir.
- Re-run this script after new L1 completions to add L2 dirs.

## Convention

- **1 candidate per (sid, ads) group** — the L1 global-minimum only.
- Restart from L1 CONTCAR to avoid convergence instability that occurs
  if LSOL is turned on from the raw initial guess.

## Naming

- v4 T1_16_DFT_L2/ dir kept as-is (misnamed — L2 here is shortlist version, not solvation level).
- T1_17_VASPsol/ is the actual solvation level.

Manifest: {len(manifest)} rows total, {created} L2 dirs created, {skipped} skipped (existed).
"""
    (OUT/'README.md').write_text(readme)
    print(f'\nT1_17_VASPsol: {created} created, {skipped} skipped, '
          f'{len([m for m in manifest if m["status"]=="PENDING L1 completion"])} PENDING')

if __name__=='__main__':
    main()
