# Project Status — Pd/PdO/PdO₂ DMC Formation DFT Study

Last updated: 2026-06-04

## Current position: Pre-G1 (환경 세팅 완료, bulk 최적화 대기)

## Environment checklist

| Item | Status | Note |
|------|--------|------|
| conda env `pddmc` (Python 3.11) | OK | all 10 packages import verified |
| ase | 3.28.0 | |
| pymatgen | 2026.5.4 | |
| mp-api (MPRester) | 0.46.1 | |
| rdkit | 2026.03.2 | |
| autoadsorbate | 0.2.5 | Surface, Fragment OK |
| mace-torch | 0.3.16 | mace_mp calculator OK |
| numpy | 2.4.6 | |
| scipy | 1.17.1 | |
| pandas | 3.0.3 | |
| matplotlib | 3.10.9 | |
| MP_API_KEY | SET | value not displayed |
| sbatch / squeue | /usr/bin/sbatch, /usr/bin/squeue | SLURM available |
| VASP 6.4.3 (vanilla) | `/home/hyunjin/vasp.6.4.3/bin/vasp_std` | std/gam/ncl all present |
| VASP 6.4.3 + VTST | `/home/hyunjin/vasp.6.4.3_vtst/bin/vasp_std` | NEB/dimer/climbing confirmed |
| VASPsol | **AVAILABLE** in both builds | solvation.o compiled; `Solvation Ediel_sol` confirmed |
| VTST scripts | `/home/hyunjin/VTST/vtstscripts-1040/` | nebmake.pl, dimer scripts, etc. |
| git repo | initialized | initial commit ffde43c |
| GitHub repo | git@github.com:hyunjinlee12/DMC.git | SSH, pushed to main |

## Gate status

| Gate | Status | Criteria |
|------|--------|----------|
| G1 | PENDING | 3 bulk structures converged + lattice vs experiment |
| G2 | BLOCKED by G1 | 4 clean slabs validated |
| G3 | BLOCKED by G2 | descriptor map + Phase-2 surface selection |
| G4 | BLOCKED by G3 | DMC Gibbs profile + conclusion |

## Next: G1 plan (T1.1–T1.4)

See below for the detailed execution plan. Awaiting researcher approval to start DFT.
