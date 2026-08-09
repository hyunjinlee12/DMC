# T1_17_VASPsol — Level-2 solvation setup

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

Manifest: 14 rows total, 2 L2 dirs created, 0 skipped (existed).
