# gas_references — isolated-molecule DFT references for T1.18

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
