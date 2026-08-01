# v5 addendum (2026-07-17 refinements)

## 1. Errata

- **v4 S1/CO composition** was mis-summarized as "5 br_PdPd" in the earlier text.
  The v4 picks were **atop_Pd 1 (idx 64) + br_PdPd 4 (idx 3, 8, 28, 30)**
  per `picks_analysis.csv`. The v5 selection log itself was correct; only the
  narrative summary was wrong.

## 2. Pool count clarification

| sid | ads | raw (AutoAdsorbate) | unique_json (MLIP dedup) | converged | valid (final v5 pool) |
|---|---|---|---|---|---|
| S1 | CO | 112 | 12 | 12 | 12 |
| S1 | CH3O | 336 | 85 | 85 | 85 |
| S1 | coads | 6338 | 1985 | 1985 | 1985 |
| S2 | CO | 166 | 26 | 26 | 26 |
| S2 | CH3O | 498 | 321 | 321 | 311 |
| S2 | coads | 11598 | 5524 | 5524 | 5334 |
| S3 | CO | 135 | 10 | 10 | 10 |
| S3 | CH3O | 405 | 172 | 172 | 168 |
| S3 | coads | 8281 | 3614 | 3614 | 3484 |
| S3b | CO | 87 | 18 | 18 | 18 |
| S3b | CH3O | 261 | 189 | 189 | 189 |
| S3b | coads | 3764 | 1581 | 1581 | 1579 |
| S4 | CO | 129 | 52 | 50 | 26 |
| S4 | CH3O | 387 | 291 | 290 | 165 |

**Drop reasons**:

- `unique_json → converged`: MLIP relaxation did not reach EDIFFG stop → excluded.
- `converged → valid`: intramolecular geometry check (`valid_CO/CH3O/coads`) rejects
  broken adsorbates — C–O bond outside allowed range, CH₃ lost an H, or a 2nd O
  is within 2.0 Å of the C (CO₂-like collapse). This is the dominant loss on **S4/CO
  (52 → 26)** — half of the "converged" candidates collapsed to CO₂-like on the O-rich
  PdO₂(110) surface. `S2/CH3O 321→311`, `S3/CH3O 172→168` = a handful of dissociated methoxides.

## 3. Refined additions (proposed_additions_only.csv)

Priority sorted: MUST → RECOMMENDED → REVIEW-NEEDED → OPTIONAL.
Each row includes E_MLIP, ΔE(global), site+region, neighbor coordination,
distance bin (coads), xy MIC distance to nearest v4 pick, and neighbor distance table.

### MUST (4)
- **S3/coads idx=3481** — E=-713.2831 eV, ΔE=0.0414 eV, d_reactive=2.988 Å (reactive-close), combo atop_Pd+atop_Pd. xy dist to nearest v4 = 0.21 Å.
    reason: distance-bin gap (reactive) not covered by v4
- **S3/coads idx=5161** — E=-713.3037 eV, ΔE=0.0209 eV, d_reactive=3.249 Å (reactive-loose), combo atop_Pd+atop_Pd. xy dist to nearest v4 = 0.386 Å.
    reason: distance-bin gap (reactive) not covered by v4
- **S3b/coads idx=2051** — E=-574.1964 eV, ΔE=0.5542 eV, d_reactive=2.895 Å (reactive-close), combo br_PdPd+atop_Pd. xy dist to nearest v4 = 3.473 Å.
    reason: distance-bin gap (reactive) not covered by v4
- **S3b/coads idx=2754** — E=-574.6289 eV, ΔE=0.1217 eV, d_reactive=3.663 Å (reactive-loose), combo br_PdPd+br_PdPd. xy dist to nearest v4 = 4.138 Å.
    reason: distance-bin gap (reactive) not covered by v4

### RECOMMENDED (4)
- **S1/coads idx=487** — E=-432.0193 eV, ΔE=0.0318 eV, d_reactive=2.989 Å (reactive-close), combo br_PdPd+br_PdPd. xy dist to nearest v4 = 0.596 Å.
    reason: distinct-site br_PdPd+br_PdPd + reactive-close
- **S1/coads idx=4021** — E=-432.0405 eV, ΔE=0.0107 eV, d_reactive=5.402 Å (thermodynamic), combo br_PdPd+br_PdPd. xy dist to nearest v4 = 1.974 Å.
    reason: thermodynamic representative
- **S2/coads idx=3633** — E=-600.9994 eV, ΔE=0.2289 eV, d_reactive=5.186 Å (thermodynamic), combo br_PdPd+atop_Pd. xy dist to nearest v4 = 2.943 Å.
    reason: thermodynamic representative
- **S2/coads idx=8079** — E=-600.7868 eV, ΔE=0.4414 eV, d_reactive=2.858 Å (reactive-close), combo atop_Pd+atop_Pd. xy dist to nearest v4 = 1.361 Å.
    reason: distinct-site atop_Pd+atop_Pd + reactive-close

### REVIEW-NEEDED (auto-added coads picks NOT on researcher list)
- S1/coads idx=1696 — separated combo br_PdPd+br_PdPd, ΔE=0.0028 eV. *Not a reactive rep — likely redundant with global min or another separated/thermo pick. RECOMMENDED action: drop unless the site combo is genuinely novel.*
- S1/coads idx=3255 — thermodynamic combo 4f_4Pd0O+br_PdPd, ΔE=0.0483 eV. *Not a reactive rep — likely redundant with global min or another separated/thermo pick. RECOMMENDED action: drop unless the site combo is genuinely novel.*
- S2/coads idx=9507 — separated combo br_PdPd+atop_Pd, ΔE=0.2008 eV. *Not a reactive rep — likely redundant with global min or another separated/thermo pick. RECOMMENDED action: drop unless the site combo is genuinely novel.*
- S2/coads idx=11165 — separated combo atop_Pd+br_PdPd, ΔE=0.2187 eV. *Not a reactive rep — likely redundant with global min or another separated/thermo pick. RECOMMENDED action: drop unless the site combo is genuinely novel.*
- S3/coads idx=4593 — separated combo physi+atop_Pd, ΔE=0.0691 eV. *Not a reactive rep — likely redundant with global min or another separated/thermo pick. RECOMMENDED action: drop unless the site combo is genuinely novel.*
- S3b/coads idx=2484 — separated combo br_PdPd+br_PdPd, ΔE=0.0063 eV. *Not a reactive rep — likely redundant with global min or another separated/thermo pick. RECOMMENDED action: drop unless the site combo is genuinely novel.*
- S3b/coads idx=2841 — thermodynamic combo atop_Pd+br_PdPd, ΔE=0.1335 eV. *Not a reactive rep — likely redundant with global min or another separated/thermo pick. RECOMMENDED action: drop unless the site combo is genuinely novel.*

### OPTIONAL (single-ads v5-new — 다양성 이유이지만 ΔE 검토 필요)
- S1/CH3O idx=283 — site=atop_Pd/metal, ΔE=0.0024 eV, d_anchor_surf=1.987 Å. neighbors: Pd:1.99;Pd:3.31;Pd:3.35...
- S1/CO idx=43 — site=4f_4Pd0O/metal, ΔE=0.0464 eV, d_anchor_surf=2.154 Å. neighbors: Pd:2.15;Pd:2.17;Pd:2.18;Pd:2.18;Pd:3.09...
- S2/CH3O idx=217 — site=br_PdO/interface, ΔE=0.2877 eV, d_anchor_surf=2.032 Å. neighbors: Pd:2.03;O:2.59;O:2.98...
- S2/CH3O idx=496 — site=atop_O/oxide, ΔE=0.3223 eV, d_anchor_surf=2.25 Å. neighbors: O:2.25;Pd:3.33...
- S2/CO idx=14 — site=br_PdPd/metal, ΔE=0.2951 eV, d_anchor_surf=2.02 Å. neighbors: Pd:2.02;Pd:2.02;O:2.71;O:2.73;O:2.87;O:2.95;Pd:3.24;Pd:3.47...
- S3/CH3O idx=385 — site=br_OO/oxide, ΔE=0.3678 eV, d_anchor_surf=2.523 Å. neighbors: O:2.52;O:2.53;Pd:3.48...
- S3/CH3O idx=395 — site=physi/physi, ΔE=0.0292 eV, d_anchor_surf=2.632 Å. neighbors: O:2.63;O:2.91...
- S3/CO idx=76 — site=physi/physi, ΔE=0.105 eV, d_anchor_surf=3.196 Å. neighbors: O:3.20;O:3.25...
- S3b/CO idx=6 — site=br_PdPd/metal, ΔE=0.134 eV, d_anchor_surf=1.918 Å. neighbors: Pd:1.92;Pd:2.01;O:3.32;O:3.44...
- S3b/CO idx=84 — site=physi/physi, ΔE=0.0051 eV, d_anchor_surf=3.59 Å. neighbors: ...
- S4/CH3O idx=329 — site=br_PdO/interface, ΔE=0.0419 eV, d_anchor_surf=1.998 Å. neighbors: Pd:2.00;O:2.55;O:2.72;O:2.94;O:3.01;O:3.36;O:3.45...
- S4/CH3O idx=383 — site=physi/physi, ΔE=0.0063 eV, d_anchor_surf=2.809 Å. neighbors: O:2.81;O:2.95;O:3.25...

## 4. Site verification (anchor-neighbor distance table)

Each MUST/RECOMMENDED candidate below shows all slab atoms within 3.5 Å of
the anchor, so the site label is grounded in actual coordination:

**S3/coads idx=3481** (site_combo=atop_Pd+atop_Pd)
  neighbors within 3.5 Å: C_CO:[Pd:2.47;O:3.07;O:3.08;O:3.12;O:3.12] O_CH3:[Pd:2.00;O:2.75;O:2.76;O:2.86;O:2.87]

**S3/coads idx=5161** (site_combo=atop_Pd+atop_Pd)
  neighbors within 3.5 Å: C_CO:[Pd:2.49;O:2.98;O:3.08;O:3.15;O:3.25] O_CH3:[Pd:2.00;O:2.80;O:2.81;O:2.82;O:2.84]

**S3b/coads idx=2051** (site_combo=br_PdPd+atop_Pd)
  neighbors within 3.5 Å: C_CO:[Pd:1.89;Pd:2.11;O:3.18;O:3.26;Pd:3.31] O_CH3:[Pd:2.01;Pd:2.97]

**S3b/coads idx=2754** (site_combo=br_PdPd+br_PdPd)
  neighbors within 3.5 Å: C_CO:[Pd:1.93;Pd:1.98] O_CH3:[Pd:2.01;Pd:2.18;O:3.25]

**S1/coads idx=487** (site_combo=br_PdPd+br_PdPd)
  neighbors within 3.5 Å: C_CO:[Pd:1.98;Pd:1.98;Pd:3.18;Pd:3.19] O_CH3:[Pd:2.08;Pd:2.10;Pd:3.37;Pd:3.42]

**S1/coads idx=4021** (site_combo=br_PdPd+br_PdPd)
  neighbors within 3.5 Å: C_CO:[Pd:1.97;Pd:1.97;Pd:3.40;Pd:3.41;Pd:3.44;Pd:3.45] O_CH3:[Pd:2.10;Pd:2.11;Pd:3.41;Pd:3.42]

**S2/coads idx=3633** (site_combo=br_PdPd+atop_Pd)
  neighbors within 3.5 Å: C_CO:[Pd:2.01;Pd:2.03;O:2.73;O:2.74;O:2.89;O:2.94;Pd:3.17] O_CH3:[Pd:2.02;O:2.70;O:2.81;Pd:3.19;Pd:3.30]

**S2/coads idx=8079** (site_combo=atop_Pd+atop_Pd)
  neighbors within 3.5 Å: C_CO:[Pd:1.90;O:2.68;O:2.83;Pd:3.49] O_CH3:[Pd:2.03;O:2.73;O:2.79;Pd:3.21;Pd:3.23]

## 5. S4 coadsorption

Still intentionally excluded per researcher decision.

## 6. What was NOT done (guardrails)

- No file under `T1_16_DFT_L2/` or the pending65 bundle was touched.
- No new calc directories were created.
- No jobs submitted or cancelled.
- Existing S1/CO 5 completed candidates (idx 64, 3, 8, 28, 30) untouched;
  v4-only idx=30 is NOT replaced.
