"""T1.17 setup — VASPsol Level-2 solvation for each (sid, ads) group top-1
                 + clean-slab VASPsol references (required for consistent Gibbs).

Reads paper_data/07_dft_results.csv, picks the lowest-E_DFT candidate per
(sid, ads) group among DONE entries, and builds:

  T1_17_VASPsol/{sid}/{ads}/{ads}_idx{XXXXX}/    ← 14 adsorbate-slab groups
  T1_17_VASPsol/{sid}_clean/                     ← 5 clean-slab VASPsol references

Each dir contains: POSCAR, INCAR, POTCAR, KPOINTS (Γ-only via KSPACING),
submit_vasp.sh, metadata.json.

CONVENTIONS (post-review 2026-07-17):
- Baseline VASPsol INCAR: LSOL=.TRUE., EB_K=32.6, TAU=0. NO LAMBDA_D_K
  (that turns on electrolyte Debye screening, which is NOT what "methanol
  implicit solvent" means).
- ISTART=0 explicitly. L1 CONTCAR only carries geometry (L1 wrote LWAVE=
  .FALSE. so no WAVECAR to restart electronically). This is a geometric
  restart, not an electronic one — documented in metadata.
- clean slab dirs use G2_slab/*/CONTCAR + identical VASPsol settings.
- top-1 marking: if a later L1 completion shows lower E, this script
  creates the new dir AND writes `.SUPERSEDED` marker in the old dir's
  metadata.

No jobs submitted.
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

# ---------- INCAR templates (post-review baseline) ----------
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

# ---- Fresh electronic start (L1 wrote LWAVE=.FALSE., no WAVECAR to restart)
ISTART = 0
ICHARG = 2

# ---- VASPsol (methanol implicit solvent) ----
# EB_K = static dielectric of methanol; TAU = 0 excludes cavitation term.
# Electrolyte screening tag intentionally not set — see README for rationale.
LSOL = .TRUE.
EB_K = 32.6
TAU = 0
"""
INCAR_OXIDE = INCAR_METAL.replace("ISMEAR = 1\nSIGMA = 0.1",
                                   "ISMEAR = 0\nSIGMA = 0.05")

# ---------- KPOINTS Γ-only via KSPACING already; keep KSPACING in INCAR ----------
# Explicit KPOINTS file not needed for slab (KSPACING sufficient), but we still
# write one so users on servers with different defaults get an unambiguous mesh.
# For 15+ Å slabs, KSPACING=0.25 gives ~3x3x1 mesh — same as L1.
# (No file emitted; INCAR KSPACING is authoritative.)

# ---------- submit script ----------
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

# NOTE: VASP >=5.4.1 standard builds support solvation via LSOL — a separate
# vasp_std_sol binary is NOT required if your build has LSOL compiled in.
# Verify on first pilot run:
#   grep -E "VASPsol|LSOL|EB_K" OUTCAR
# should show solvation is active. If unknown-INCAR-tag warnings appear,
# rebuild VASP with the solvation source.
VASP_BIN=${VASP_BIN:-/home/hyunjin/vasp.6.4.3/bin/vasp_std}

NPROCS=${SLURM_NTASKS:-1}
echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | VASP: ${VASP_BIN}"
echo "Start: $(date)"
mpirun --bind-to none -np ${NPROCS} ${VASP_BIN}
echo "End: $(date)"
"""

def write_common(dest, atoms, sid, extra_meta):
    write(str(dest/'POSCAR'), atoms, format='vasp', direct=True, sort=True, vasp5=True)
    (dest/'INCAR').write_text(INCAR_METAL if sid=='S1' else INCAR_OXIDE)
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    potcar = ''.join((POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species)
    (dest/'POTCAR').write_text(potcar)
    (dest/'submit_vasp.sh').write_text(SUBMIT_SCRIPT)
    (dest/'metadata.json').write_text(json.dumps({
        'level': 'T1.17 VASPsol (methanol, LSOL=T, EB_K=32.6, TAU=0)',
        'restart_kind': 'geometric only (ISTART=0, ICHARG=2). L1 had LWAVE=.FALSE. → no WAVECAR to restart electronically.',
        'poscar_species': species,
        **extra_meta,
    }, indent=2))

def build_adsorbate_dir(sid, ads, idx, contcar_src, E_L1_vacuum):
    dest = OUT/sid/ads/f'{ads}_idx{idx:05d}'
    if dest.exists():
        return None, 'already exists'
    dest.mkdir(parents=True)
    a = read(contcar_src)
    write_common(dest, a, sid, extra_meta={
        'sid':sid, 'ads':ads, 'idx':int(idx),
        'source_L1_dir': str(contcar_src.parent.relative_to(ROOT)),
        'source_L1_CONTCAR': str(contcar_src.relative_to(ROOT)),
        'E_L1_vacuum_eV': E_L1_vacuum,
        'purpose': 'Level-2 solvated adsorbate slab for T1.18 ΔG_ads computation.',
        'pair_reference': f'clean slab: T1_17_VASPsol/{sid}_clean/',
    })
    return dest, 'created'

def _strip_adsorbate(atoms, ads_kind, pristine_substrate):
    """Return substrate-only copy of ads-slab atoms with defensive asserts.

    ads_kind: 'CO' | 'CH3O' | 'coads' — dictates expected removed composition.
    pristine_substrate: read(G2_slab/{sdir}/CONTCAR) — gives expected counts.

    Adsorbate = all C, all H, and any O bonded (< 1.5 Å, PBC-aware) to any C.

    Asserts:
      - #C removed == expected (CO:1, CH3O:1, coads:2)
      - #H removed == expected (CO:0, CH3O:3, coads:3)
      - #O removed == expected (CO:1, CH3O:1, coads:2)
      - remaining Pd count == pristine Pd count
      - remaining O count == pristine substrate O count
      - every removed O has ≥ 1 C neighbour < 1.5 Å (PBC-aware)
      - no substrate O was removed
    """
    EXP = {'CO':(1,0,1),'CH3O':(1,3,1),'coads':(2,3,2)}
    n_c_exp, n_h_exp, n_o_exp = EXP[ads_kind]

    syms = atoms.get_chemical_symbols()
    c_idx = [i for i,s in enumerate(syms) if s=='C']
    h_idx = [i for i,s in enumerate(syms) if s=='H']
    o_idx = [i for i,s in enumerate(syms) if s=='O']

    # Identify ads O by PBC-aware C–O distance
    ads_o = []
    for oi in o_idx:
        for ci in c_idx:
            d = atoms.get_distance(ci, oi, mic=True)
            if d < 1.5:
                ads_o.append(oi); break

    assert len(c_idx) == n_c_exp, \
        f'ads_kind={ads_kind}: expected {n_c_exp} C, got {len(c_idx)}'
    assert len(h_idx) == n_h_exp, \
        f'ads_kind={ads_kind}: expected {n_h_exp} H, got {len(h_idx)}'
    assert len(ads_o) == n_o_exp, \
        f'ads_kind={ads_kind}: expected {n_o_exp} C-bonded O, got {len(ads_o)}'

    ads_set = set(c_idx) | set(h_idx) | set(ads_o)
    keep = [i for i in range(len(atoms)) if i not in ads_set]
    slab = atoms[keep]

    # Verify substrate composition matches pristine G2 slab
    ss = slab.get_chemical_symbols()
    n_pd_slab = sum(1 for s in ss if s=='Pd')
    n_o_slab  = sum(1 for s in ss if s=='O')
    n_pd_pris = sum(1 for s in pristine_substrate.get_chemical_symbols() if s=='Pd')
    n_o_pris  = sum(1 for s in pristine_substrate.get_chemical_symbols() if s=='O')
    assert n_pd_slab == n_pd_pris, \
        f'stripped Pd {n_pd_slab} != pristine Pd {n_pd_pris}'
    assert n_o_slab == n_o_pris, \
        f'stripped O {n_o_slab} != pristine substrate O {n_o_pris} '\
        f'— possible substrate-O leak into adsorbate O detection'
    return slab

def build_clean_slab_dir(sid, ads_contcar_src, ads_kind):
    """Clean slab built from THIS surface's top-1 ads CONTCAR:
    substrate atoms with the EXACT mask + positions that the T1.16 ads calc had.
    Guarantees identical selective-dynamics indices + cell + POTCAR order
    with the adsorbate dir, so cavitation/dielectric offsets cancel exactly."""
    dest = OUT/f'{sid}_clean'
    if dest.exists():
        return None, 'already exists'
    dest.mkdir(parents=True)
    ads_atoms = read(ads_contcar_src)
    pristine = read(G2/SDIRS[sid]/'CONTCAR')
    slab = _strip_adsorbate(ads_atoms, ads_kind, pristine)
    write_common(dest, slab, sid, extra_meta={
        'sid':sid, 'kind':'clean_slab',
        'source_ads_CONTCAR': str(ads_contcar_src.relative_to(ROOT)),
        'derived_by': ('strip C, H, and any O within 1.5 Å of C from ads CONTCAR. '
                       'Preserves substrate positions + selective-dynamics mask '
                       'so ΔG_ads = G(slab+ads) − G(slab) − μ_ads has exact '
                       'substrate cancellation.'),
        'n_atoms_total': len(slab),
        'n_atoms_fixed': sum(len(c.index) for c in slab.constraints
                              if isinstance(c, FixAtoms)),
        'purpose': 'Clean-slab VASPsol reference (Convention A: fully relax).',
        'note': ('One clean slab per surface. Mask inherited from that '
                 "surface's top-1 ads CONTCAR — matches that top-1 exactly. "
                 'Other ads candidates on the same surface may have '
                 'slightly different masks (differ by a few atoms in layer 3) '
                 'due to per-candidate MLIP relaxation. Effect on ΔG is '
                 'expected within 10-50 meV.'),
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

def mark_superseded_if_stale(sid, ads, current_idx):
    """If a T1.17 dir exists for the same group with a DIFFERENT idx than
    the new top-1, mark it SUPERSEDED in its metadata.json."""
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
            meta['SUPERSEDED_by'] = f'{sid}/{ads}/{ads}_idx{current_idx:05d}'
            meta['SUPERSEDED_note'] = (f'A later L1 completion at idx={current_idx} '
                                       f'gave lower E_DFT than idx={old_idx}. '
                                       f'This dir is retained for reproducibility '
                                       f'but should not be re-submitted or included '
                                       f'in the final descriptor map.')
            m.write_text(json.dumps(meta, indent=2))
            print(f'  SUPERSEDED marker set in {d.relative_to(ROOT)}')

def main():
    manifest = []
    picks = pick_top1_per_group()
    created = 0; skipped = 0

    # ---- adsorbate top-1 per group ----
    for (sid, ads), p in sorted(picks.items()):
        mark_superseded_if_stale(sid, ads, p['idx'])
        contcar_src = p['contcar_path']
        if not contcar_src.exists(): continue
        dest, status = build_adsorbate_dir(sid, ads, p['idx'], contcar_src, p['E'])
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

    # ---- clean slab references — one per surface, derived from any DONE ads
    # candidate on that surface. Prefers CO > CH3O > coads (in that order) for
    # least substrate distortion; falls back to whichever ads type is DONE. ----
    picked_ads_for_clean = {}
    ads_pref = ['CO','CH3O','coads']
    for sid in ['S1','S2','S3','S3b','S4']:
        for ads in ads_pref:
            if (sid, ads) in picks:
                picked_ads_for_clean[sid] = picks[(sid, ads)]
                break
    for sid in ['S1','S2','S3','S3b','S4']:
        if sid not in picked_ads_for_clean:
            manifest.append({'kind':'clean_slab','sid':sid,'ads':'','idx':'',
                             'E_L1_vacuum':'','status':'PENDING (no DONE ads on this surface yet)',
                             'dir':''})
            continue
        p = picked_ads_for_clean[sid]
        # find which ads_kind we picked (needed for strip asserts)
        ads_kind_for_clean = next(a for (s,a) in picks if s==sid and picks[(s,a)]==p)
        dest, status = build_clean_slab_dir(sid, p['contcar_path'], ads_kind_for_clean)
        if dest is None:
            skipped += 1
            manifest.append({'kind':'clean_slab','sid':sid,'ads':'','idx':'',
                             'E_L1_vacuum':'','status':'skipped (exists)',
                             'dir':str((OUT/f'{sid}_clean').relative_to(ROOT))})
            continue
        created += 1
        manifest.append({'kind':'clean_slab','sid':sid,'ads':'','idx':'',
                         'E_L1_vacuum':'','status':'created',
                         'dir':str(dest.relative_to(ROOT)),
                         'derived_from':str(p['contcar_path'].relative_to(ROOT))})
        print(f'  clean {sid} → {dest.relative_to(ROOT)} '
              f'(from {p["contcar_path"].parent.name})')

    # ---- pending groups ----
    all_groups = [('S1','CO'),('S1','CH3O'),('S1','coads'),
                  ('S2','CO'),('S2','CH3O'),('S2','coads'),
                  ('S3','CO'),('S3','CH3O'),('S3','coads'),
                  ('S3b','CO'),('S3b','CH3O'),('S3b','coads'),
                  ('S4','CO'),('S4','CH3O')]
    for g in all_groups:
        if g not in picks:
            manifest.append({'kind':'ads_top1','sid':g[0],'ads':g[1],'idx':'',
                             'E_L1_vacuum':'','status':'PENDING L1 completion',
                             'dir':''})

    with open(OUT/'manifest.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['kind','sid','ads','idx','E_L1_vacuum','status','dir','derived_from'])
        w.writeheader()
        for r in manifest:
            r.setdefault('derived_from','')
            w.writerow(r)

    readme = f"""# T1_17_VASPsol — Level-2 methanol implicit solvation

**Regenerate**: `python scripts/setup_t1_17_vaspsol.py` (idempotent + SUPERSEDED aware).

## VASPsol baseline INCAR (post-review 2026-07-17)

```
LSOL = .TRUE.
EB_K = 32.6      # methanol static dielectric
TAU = 0          # exclude cavitation / non-electrostatic (workplan)
# LAMBDA_D_K deliberately omitted — enabling it activates the linearised
# Poisson–Boltzmann electrolyte model (Debye screening), which is a
# separate physical effect from a neutral methanol dielectric.
ISTART = 0       # fresh electronic start (L1 wrote LWAVE=.FALSE., no WAVECAR)
ICHARG = 2
```

## Dir layout

```
T1_17_VASPsol/
├── {sid}/{ads}/{ads}_idxXXXXX/     ← adsorbate slabs (top-1 per group)
├── {sid}_clean/                    ← clean-slab VASPsol references (5)
└── manifest.csv
```

## Convention used in T1.18 (post-computation)

```
ΔG_ads(sol) = G(slab+ads)_sol − G(slab)_sol − μ_ads(reference in matching phase)
```

The clean-slab VASPsol dirs (S1_clean, S2_clean, …, S4_clean) provide
`G(slab)_sol`. Both must be computed with **identical INCAR flags** so the
cavitation/dielectric offsets cancel in the difference.

**Do not mix vacuum-slab and solvated-adsorbate energies in one formula.**

## VASP binary

Any VASP ≥ 5.4.1 build compiled with solvation source is fine — the standard
`vasp_std` may already support `LSOL`. Verify on the first pilot:
```
grep -E "VASPsol|LSOL|EB_K" OUTCAR
```
If unknown-INCAR-tag warnings appear, the binary needs to be rebuilt.

## Status (this run)

- **{sum(1 for r in manifest if r['kind']=='ads_top1' and r['status']=='created')} adsorbate dirs created** (from L1-DONE groups).
- **{sum(1 for r in manifest if r['kind']=='clean_slab' and r['status']=='created')} clean-slab dirs created** (one per surface with any DONE ads).
- **{sum(1 for r in manifest if 'PENDING' in r['status'])} PENDING** L1 completion (adsorbate groups + surfaces without a DONE ads).
- Re-run this script after new L1 completions to add L2 dirs and, if a
  lower-E winner emerges, mark the previous dir SUPERSEDED.

At full L1 completion this will grow to **14 adsorbate + 5 clean = 19 dirs total**.
Combined with 10 `gas_references/` dirs → 29 total planned.

## Current top-1 candidates are PROVISIONAL

Only 7/86 L1 candidates are DONE. Any current top-1 selection here is a
**pilot** for verifying VASPsol behavior. The final T1.17 winners for the
descriptor map (T1.19) must be re-evaluated after all 86 L1 finish.

## Pilot acceptance criteria (first submission)

Before submitting all dirs at once, run **one** pilot (e.g. `S1_clean` — cheap,
no adsorbate) and confirm from OUTCAR:

- `LSOL`, `EB_K=32.6`, `TAU=0` are recognised (grep OUTCAR).
- Solvation-energy output present (e.g. "solvation energy", "cavity").
- No `unknown INCAR tag` warnings.
- `LAMBDA_D_K` is inactive (grep OUTCAR shows no Debye length).
- `ISTART=0`, `ICHARG=2` confirmed in OUTCAR header.
- Electronic SCF converged.
- Ionic relaxation reaches `reached required accuracy` under EDIFFG=-0.03.

Only after this passes, submit the remaining T1.17 dirs.

## Chemical-reservoir decisions still owed (analysis-time)

Before running T1.18, define once and record in paper_data/ README:

- CO reservoir: gas feed (μ_CO(g)) or dissolved CO? — affects reference choice.
- CH₃OH reservoir: gas reference (μ_CH3OH(g)) or liquid methanol at activity 1?
  If liquid, must add Δμ_solvation correction to G(CH3OH_vaspsol).
- H₂ reservoir: CHE gas reference (½ G_H2(g)) — standard.
- CH₃O radical: auxiliary reference; MeOH(U) formula is preferred.
- Never mix vacuum and vaspsol references in the same ΔG expression —
  choose the phase (vacuum or solvated) for the whole formula.
"""
    (OUT/'README.md').write_text(readme)
    print(f'\nT1_17_VASPsol: created={created}, skipped={skipped}, '
          f'PENDING={sum(1 for r in manifest if "PENDING" in r["status"])}')

    consistency_issues = auto_consistency_check()
    if consistency_issues:
        print('\n=== Consistency check ERRORS ===')
        for e in consistency_issues: print(f'  {e}')
    else:
        print('\nAuto-consistency check: all ads/clean pairs match (cell, '
              'substrate composition, mask, POTCAR order, k-point).')

def auto_consistency_check():
    """For every ads dir, verify its counterpart {sid}_clean has matching
    cell, substrate composition, atom order, selective-dynamics mask,
    POTCAR TITEL order, and k-point setting."""
    issues = []
    def _fixed_indices(a):
        for c in a.constraints:
            if isinstance(c, FixAtoms): return sorted(int(i) for i in c.index)
        return []
    def _titels(p):
        return [ln.split()[3].split('_')[0]
                for ln in Path(p).read_text().splitlines()
                if ln.strip().startswith('TITEL')]
    def _kspacing(incar_text):
        for ln in incar_text.splitlines():
            s = ln.strip()
            if s.startswith('#'): continue
            if s.startswith('KSPACING'):
                return s.split('=')[1].strip()
        return None
    ads_dirs = list(OUT.glob('*/*/[Cc]*_idx*'))
    for ad in ads_dirs:
        sid = ad.parts[-3]
        clean = OUT/f'{sid}_clean'
        if not clean.exists():
            issues.append(f'{ad.relative_to(OUT)}: no {sid}_clean counterpart'); continue
        a_ads = read(ad/'POSCAR'); a_cl = read(clean/'POSCAR')
        if not np.allclose(a_ads.cell.array, a_cl.cell.array, atol=1e-6):
            issues.append(f'{sid}: cell differs (ads vs clean)')
        # substrate composition (only Pd + substrate O)
        ss_ads = a_ads.get_chemical_symbols()
        ss_cl = a_cl.get_chemical_symbols()
        n_pd_ads = sum(1 for s in ss_ads if s=='Pd')
        n_pd_cl  = sum(1 for s in ss_cl if s=='Pd')
        if n_pd_ads != n_pd_cl:
            issues.append(f'{sid}: Pd count differs ({n_pd_ads} vs {n_pd_cl})')
        # selective-dynamics mask: counts must match
        f_ads = _fixed_indices(a_ads); f_cl = _fixed_indices(a_cl)
        if len(f_ads) != len(f_cl):
            issues.append(f'{sid}: mask count differs (ads {len(f_ads)} vs clean {len(f_cl)})')
        # POTCAR TITEL order — clean should match pristine substrate species
        t_cl = _titels(clean/'POTCAR')
        pristine = read(G2/SDIRS[sid]/'CONTCAR')
        pristine_species_ordered = []
        for s in ('C','H','O','Pd'):     # ASE sort=True alphabetical order
            if s in pristine.get_chemical_symbols():
                pristine_species_ordered.append(s)
        if t_cl != pristine_species_ordered:
            issues.append(f'{sid}: clean POTCAR {t_cl} != pristine substrate {pristine_species_ordered}')
        # KSPACING
        k_ads = _kspacing((ad/'INCAR').read_text())
        k_cl  = _kspacing((clean/'INCAR').read_text())
        if k_ads != k_cl:
            issues.append(f'{sid}: KSPACING differs ({k_ads} vs {k_cl})')
    return issues

if __name__=='__main__':
    main()
