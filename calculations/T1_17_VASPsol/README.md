# T1_17_VASPsol — Level-2 methanol implicit solvation

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
├── S4/CH3O/CH3O_idxXXXXX/     ← adsorbate slabs (top-1 per group)
├── S4_clean/                    ← clean-slab VASPsol references (5)
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

## Status

- **2 adsorbate dirs created** (from L1-DONE groups).
- **5 clean-slab dirs created**.
- **12 PENDING** L1 completion.
- Re-run this script after new L1 completions to add L2 dirs and, if a
  lower-E winner emerges, mark the previous dir SUPERSEDED.

## Current top-1 candidates are PROVISIONAL

Only 7/86 L1 candidates are DONE. Any current top-1 selection here is a
**pilot** for verifying VASPsol behavior (LSOL recognized, solvation output
present, no unknown-tag warnings). The final T1.17 winners for the descriptor
map (T1.19) must be re-evaluated after all 86 L1 finish.
