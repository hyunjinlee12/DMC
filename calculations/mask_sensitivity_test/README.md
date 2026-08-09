# Mask sensitivity test — S1 only

Reviewer 2026-07-17 round 7: compare partial-layer (T1.16) vs complete-layer
(G2) mask on S1 to quantify the effect on ΔE_ads.

## 4 dirs
| dir | mask | n_fixed | source |
|---|---|---|---|
| S1_CO_idx3_40fix_partial   | T1.16 partial (median z) | 40 | T1.16 top-1 CONTCAR |
| S1_CO_idx3_32fix_complete  | G2 complete-layer (L1+L2) | 32 | T1.16 top-1 CONTCAR |
| S1_clean_40fix_partial     | partial (median z) on G2  | 40 | G2 CONTCAR |
| S1_clean_32fix_complete    | G2 complete-layer         | 32 | G2 CONTCAR verbatim |

All 4: vacuum PBE-D3, IBRION=2 NSW=200, KSPACING=0.25, ISPIN=2 (default MAGMOM),
INCAR otherwise identical to T1.16.

## Analysis
After all 4 finish:
```
ΔE_ads(40) = E(S1_CO_40fix)    - E(S1_clean_40fix)    - μ_CO(gas)
ΔE_ads(32) = E(S1_CO_32fix)    - E(S1_clean_32fix)    - μ_CO(gas)
Δ           = ΔE_ads(40) − ΔE_ads(32)
```

- |Δ| < 10 meV → mask choice is immaterial for downstream (use either)
- |Δ| ≥ 10 meV → complete-layer mask (32) is the correct choice
