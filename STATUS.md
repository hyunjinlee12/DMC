# Project Status — Pd/PdO/PdO₂ DMC Formation DFT Study

Last updated: 2026-06-16

## Current position: G2 PASSED ✓ — G3 sampling phase (T1.13 done, T1.14 advisor checkpoint)

## G2 slab status (T1.5–T1.9) — ALL CONVERGED

| Slab | Model | Atoms | E (eV) | Status |
|------|-------|-------|--------|--------|
| S1 | Pd(100) | 80 | -434.408 | CONVERGED ✓ |
| S2 | 1 ML PdO(101)/Pd(100) (√5×√5)R27° | 112 | -618.565 | CONVERGED ✓ |
| S3 | PdO(100) O-term | 128 | -724.152 | CONVERGED ✓ |
| S3b | PdO(100) PdO-term | 104 | -570.770 | CONVERGED ✓ |
| S4 | PdO₂(110) | 144 | -788.493 | CONVERGED ✓ |

T1.9 validation passed for all 5 surfaces (rumpling, Pd–O, coordination, dipole, termination stability) → **checkpoint B cleared**.

## G3 sampling status (T1.10–T1.13) — HEURISTIC CANDIDATES READY

Heuristic AutoAdsorbate enumeration done for 5 surfaces.

| Slab | sites | CO* | CH₃O* | Set A | Set TS | Set B | side | rejected |
|------|------:|----:|------:|------:|-------:|------:|-----:|---------:|
| S1   | 112   | 112 | 336   | 6,338 | 1,097  | 9,600 | 95   | 3,099    |
| S2   | 166   | 166 | 498   | 11,598| 1,625  | 29,147| 183  | 4,935    |
| S3   | 135   | 135 | 405   | 8,281 | 1,441  | 16,110| 177  | 3,998    |
| S3b  | 87    | 87  | 261   | 3,764 | 469    | 7,157 | 78   | 1,301    |
| S4   | 129   | 129 | 387   | 7,975 | 1,182  | 14,270| 120  | 4,412    |
| **Σ**| **629** | **629** | **1,887** | **37,956** | **5,814** | **76,284** | **653** | **17,745** |

Set A count exceeds guide target 150–500 by 10–25× — advisor decision pending on sampling option (A: keep + MLIP dedup, B: stride, C: distance-bin stratified). Awaiting advisor reply before T1.14.

**G3 gate NOT YET passed** — need T1.14 MLIP rank + T1.15 DFT shortlist + T1.16–17 DFT + T1.18 ads E table + T1.19 descriptor map + T1.20 Case A–D classification.

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
| MP_API_KEY | SET | value not displayed |
| sbatch / squeue | /usr/bin/sbatch, /usr/bin/squeue | SLURM available |
| VASP 6.4.3 (vanilla) | `/home/hyunjin/vasp.6.4.3/bin/vasp_std` | std/gam/ncl all present |
| VASP 6.4.3 + VTST | `/home/hyunjin/vasp.6.4.3_vtst/bin/vasp_std` | NEB/dimer/climbing confirmed |
| VASPsol | **AVAILABLE** in both builds | solvation.o compiled |
| VTST scripts | `/home/hyunjin/VTST/vtstscripts-1040/` | nebmake.pl, dimer scripts, etc. |
| git repo | initialized | branch `g2-passed-checkpoint` preserves G2 state |
| GitHub repo | git@github.com:hyunjinlee12/DMC.git | SSH, branches pushed |

## Gate status

| Gate | Status | Criteria |
|------|--------|----------|
| G1 | PASSED ✓ | 3 bulk structures converged + lattice validated |
| G2 | PASSED ✓ | 5 clean slabs converged + T1.9 validation passed (checkpoint B) |
| G3 | IN PROGRESS | T1.10–T1.13 done (candidates); T1.14–T1.20 pending advisor decision on Set A reduction |
| G4 | BLOCKED by G3 | DMC Gibbs profile + conclusion |

## Next

1. Advisor reply on Set A sampling option (A / B / C).
2. T1.14: MACE relax + post-relax dedup (E + structural fingerprint) → ranking.
3. T1.15: DFT shortlist 3–5 per surface (S2: 5–8).
4. T1.16–17: Level 1 (vacuum) → Level 2 (VASPsol) DFT.
5. T1.18–20: descriptor map + Case A–D classification → G3 gate.
