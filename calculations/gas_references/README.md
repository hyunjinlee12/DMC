# gas_references — isolated-molecule DFT (T1.18 references)

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
