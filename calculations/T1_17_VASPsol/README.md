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
├── S4/coads/coads_idxXXXXX/     ← adsorbate slabs (top-1 per group)
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

## Status (this run)

- **2 adsorbate dirs created** (from L1-DONE groups).
- **2 clean-slab dirs created** (one per surface with any DONE ads).
- **15 PENDING** L1 completion (adsorbate groups + surfaces without a DONE ads).
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
