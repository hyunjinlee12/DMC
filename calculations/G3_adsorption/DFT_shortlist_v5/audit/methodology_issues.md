# Methodology issues found during audit

Read-only survey. Ordered by severity. Nothing was fixed automatically —
each item is a candidate for discussion / a follow-up PR / a caveat in
the final report.

## Severity: high (affects interpretation)

### H1. Two different `site_type` fields are called the same name

- `MLIP_phase1/unique_*.json:site_type` is produced by
  `scripts/run_mace_phase1.py::classify_site()` using an **adaptive**
  neighbor buffer `overall_min + 0.5 Å`. Categories: `atop / bridge /
  hollow / unknown`.
- `picks_analysis.csv:site` and everything downstream (`select_top5_v5.py`,
  `refine_v5_additions.py`, `audit_selection.py`) uses a **fixed 2.6 Å**
  neighbor cutoff with more granular labels (`atop_Pd / atop_O / br_PdPd /
  br_PdO / br_OO / h3(…) / 4f(…) / physi`).
- The two labels can disagree — e.g. an S1 CH₃O candidate is `unknown`
  in the raw pool (adaptive rule fired 0 neighbors within `overall_min +
  0.5 Å` because the min distance was >2.6 Å) but relabelled as `br_PdPd`
  by v5.
- **Action**: rename the raw MLIP field to `site_type_mlip` (or
  `raw_site_hint`) to disambiguate. Document that v5's site label is the
  authoritative one for selection.

### H2. Reference for CH₃O binding is the *radical*, not the MeOH(U) form

- The reference file
  `calculations/G3_adsorption/mace_d3_references.json` stores both
  `CH3O_radical` and `CH3O_ref = E(CH3OH) − ½ E(H2)`.
- Every violin plot and E_bind CSV in the audit uses `CH3O_radical`.
- Workplan §T1.18 requires the descriptor map to be reported against
  **both** references (radical and MeOH(U)); the MeOH(U) version needs
  the applied potential `eU` term, which is not yet computed.
- **Action**: extend the audit / final report to emit `E_bind_MLIP_rad`
  **and** `E_bind_MLIP_MeOH(U)` (=`E_bind_MLIP_rad − E(CH3O_ref) +
  E(CH3O_radical) + eU`) side-by-side. Not required for selection but
  required for T1.19.

### H3. `d_CH3O_bond` column mapping bug in the original picks_analysis

Fixed in v5's `refine_v5_additions.py` output but the source
`calculations/G3_adsorption/DFT_shortlist_v3/picks_analysis.csv` still
carries the CH₃O C–O bond length in the `d_CO_bond` column instead of a
separate `d_CH3O_bond`. The v5 audit CSVs are already correct; the v3
file is preserved as historical evidence.

## Severity: medium (may bias future extensions)

### M1. No adsorbate orientation descriptor is used in selection

- CO tilt vs surface normal, CH₃O tilt, and O–C=O bend for coads are not
  recorded anywhere. Two candidates with the same anchor xy and same
  E_MLIP could differ in orientation, but the selector treats them as
  identical (see §6 of the methodology audit).
- Impact so far: minimal, because MLIP LBFGS relaxation from the same
  initial site usually produces the same orientation. But a Phase 2
  extension will need this: reactive TS1 endpoints require O–CH₃O axis
  pointing towards C_CO, not away.
- **Action**: add `θ_CO_axis`, `θ_CH3O_tilt`, `θ_CO_CH3O_alignment` to
  the pool descriptor when doing T2.1 sampling.

### M2. Distance-bin boundary at 3.0 Å has no external reference

- `docs/DMC_Pd_workplan.md` §9 defines the reactive range as `2.1–4.0 Å`
  and the thermodynamic reference as `≥5.0 Å`.
- The 3.0 Å split into `reactive-close` and `reactive-loose` is a v5
  internal convention; it happens to fall near a natural bimodal
  minimum in the S3/S3b pools but is not literature-motivated.
- Impact: shifting the boundary to 2.8 or 3.2 Å moves 2/8 coads picks
  between bins (see `selection_sensitivity.csv:threshold-sensitive-bin`
  count = 2).
- **Action**: report the bin distribution using both `2.8/3.2` and
  `3.0` boundaries in the final descriptor map; either choose one and
  justify, or report as `d_reactive` continuous variable.

### M3. `energy_cap = 0.5 eV` for single-ads diversity is arbitrary

- The 0.5 eV threshold was chosen to be well above the ~0.05 eV MLIP
  noise floor but well below the 1 eV chemisorbed-vs-physi gap on
  PdO₂(110). It is not derived from a published benchmark.
- Consequence: S4/CO chemisorbed candidates at ΔE ≈ 2.5–2.7 eV are
  excluded. If DFT re-relaxes some physi candidate into a chemisorbed
  state, we'd have missed the correct sample.
- **Action**: after DFT results come back, compare `E_bind_DFT` for the
  physi picks to what would have been the chemisorbed candidates; if
  the gap remains large, the 0.5 eV cap was correct.

### M4. xy dedup does not compare orientation-swapped pairs

- For coads, xy MIC of both anchors must be <1.5 Å for a duplicate call.
  A pair `(C_CO at A, O_CH3 at B)` and `(C_CO at B, O_CH3 at A)` would
  NOT be flagged as duplicate even though it is the same physical
  coadsorption topologically (swap CO and CH₃O positions).
- Impact: could bloat the coads pool with swap-equivalents. Not
  currently a problem because none of the 16 v5-new picks appears twice
  under swap (checked visually), but the code has no protection.
- **Action**: symmetrize the coads dedup key over adsorbate swap.

## Severity: low (documentation clarity)

### L1. `manifest.csv` `nearest_v4_idx` is only meaningful for v5-new rows

The v4 rows carry blank in this column. Fine, but easy to miss.

### L2. Priority tags in `combined_summary.csv` — some ambiguity

For v4-existing rows, `priority` is set to `v4-baseline` (in the
manifest.csv output of `setup_v5add16.py`). For v5-new rows it is
`MUST / RECOMMENDED / MUST-diagnostic / OPTIONAL`. There is no
`REVIEW-NEEDED` row in v5add16 because those 7 candidates were dropped
by the researcher before bundle build. Documented in the audit for
future readers.

### L3. `submit_vasp_gpu.sh` paths still point at the local NVHPC cluster

The bundle ships the same submit script as `pending65`. Users on other
servers (H200 etc.) must edit `VASP_BIN`, `NVHPC=$HOME/nvhpc`,
`--partition=debug`, `--gres=gpu:1`. This is by design (pending65 uses
the same convention) but should be repeated in the README of any
future add-on bundle.

## Assumptions not empirically verified

1. That MACE-MH+D3(BJ) `oc20_usemppbe` head is a good ranking surrogate
   for VASP PBE+D3(BJ) on Pd/PdO/PdO₂ chemistry. Small-scale benchmark
   (`reports/G3/…`) supports this on S1 for CO but coverage on S3/S4 is
   thinner.
2. That LBFGS `fmax=0.05` is a tight enough convergence for MLIP-side
   dedup. Anecdotally OK; not systematically compared to `fmax=0.01`.
3. That the CH₃O radical reference is the "correct" reference for
   binding-energy convention comparisons across the DMC literature.
   Cross-check needed with Shi 2024 Angew.

## Remaining known limitations of the current shortlist

1. **S1 pool depth is small** (CO*: 12 unique, CH₃O*: 85, coads: 1985).
   Diversity is bounded above by the pool, not by the selector.
2. **All CH₃O raw pool site labels are `unknown`** in the raw JSON — the
   MLIP adaptive classifier fails on tilted methoxide anchors. v5's
   re-classification is authoritative but reflects only the frozen
   MLIP-relaxed geometry; DFT re-relaxation can produce different tilt.
3. **S3b coads reactive-close only has one candidate** (idx 2051) and
   it is 0.55 eV above the global min. If DFT rejects it, T2.5 has no
   S3b starting endpoint.
4. **`product-like` bin is empty on every surface** — MLIP LBFGS does
   not spontaneously form the C–O bond. This is expected but means Phase
   2 will need dedicated dimer / NEB sampling to find CH₃OCO* precursors.
5. **No cross-check with an independent MLIP** (SevenNet was
   benchmarked once but not used for parallel ranking of the final
   shortlist). If the descriptor map shows surprising Case A–D
   assignments after DFT, an MLIP-vs-MLIP disagreement analysis would
   help diagnose whether the pool coverage was the bottleneck.

## What was NOT touched by this audit

- No calc dir was created, deleted, or modified.
- No POSCAR/INCAR/POTCAR/submit script was regenerated.
- No job was submitted or cancelled.
- No v4/v5 CSV or existing script was edited.
- Read-only files created under `calculations/G3_adsorption/DFT_shortlist_v5/audit/`:
  - `selection_methodology_audit.md`
  - `selection_formulae.md`
  - `selection_candidate_rationale.csv`
  - `selection_sensitivity.csv`
  - `selection_group_summary.csv`
  - `methodology_issues.md`
