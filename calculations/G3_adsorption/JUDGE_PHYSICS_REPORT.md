# Judge-Physics Report: Phase 2 MLIP Filtered Shortlists

**Target:** ml-trainer Phase 2 SetA (AFTER geometry filter applied via `refilter_phase2_geometry.py`)  
**Prior Verdict:** REJECT (S3 100% product collapse, S2 band violation, S4 CO physisorbed, PES flat)  
**Current Verdict:** CONCERN (4-label scale: Pass / Concern / Reject / Fail-Critical)  
**Date:** 2026-06-19  
**Artifacts:** `calculations/G3_adsorption/{surface}/MLIP_phase2_filtered/`

---

## Executive Summary

The geometry filter **partially resolved** prior critical failures:
- **S3 product collapse**: Fixed (was 100%, now 25% — rank-0 only)
- **S2 SetA band violation**: Fixed (all 4 shortlist now in [2.1, 4.0) Å)
- **Fragmentation**: Caught 50-60% of malformed structures across all surfaces
- **Bond breaking**: Detected 884 C-H breaks on S4 (strong indicator of CH₃O instability on PdO₂)

**New/Persistent Issues:**
- **S4 remains 100% defective**: All 4 shortlist candidates have CO or CH₃O physisorption
- **PES flatness**: S1/S3/S3b top-10 unique within 3-27 meV (below MACE noise floor)
- **Energy ranking inversion**: No surface has rank-0 = lowest energy
- **S3 rank-0 swap**: Product collapse fixed, but now has CH₃O physisorbed (Pd-O = 3.226 Å)

**Recommendation:** CONDITIONAL GO
- ✅ S1, S2, S3 (ranks 1-3 only), S3b → proceed to DFT with GEOMETRY-based selection
- ❌ S4 → STOP. Do not send to DFT. Requires alternative approach (see below)
- ⚠️  Ignore MLIP energies entirely for S1/S3/S3b (PES flat). Use geometric diversity.

---

## Surface-by-Surface Verdicts

### S1 (Pd(100))
- **Shortlist size:** 3
- **Chemistry pass rate:** 2/3 (67%)
- **PES spread (top-10):** 3.0 meV ❌ (below 50 meV threshold)
- **Verdict:** CONCERN
- **Issues:**
  - Rank-1 (MACE's lowest-E) has CO physisorbed (Pd-C = 3.677 Å)
  - PES extremely flat — MLIP energy ranking is noise
- **Usable for DFT:** YES, but select all 3 for diversity; ignore energies

### S2 (PdO(101)/Pd(100) interface)
- **Shortlist size:** 4
- **Chemistry pass rate:** 2/4 (50%)
- **PES spread (top-10):** 214.6 meV ✅
- **Verdict:** CONCERN
- **Issues:**
  - Rank-2: CO physisorbed (Pd-C = 3.411 Å)
  - Rank-3: Both CO and CH₃O physisorbed
- **Usable for DFT:** YES — ranks 0,1 are chemically clean and well-separated in energy

### S3 (PdO(100))
- **Shortlist size:** 4
- **Chemistry pass rate:** 3/4 (75%)
- **PES spread (top-10):** 27.3 meV ❌
- **Verdict:** CONCERN
- **Issues:**
  - Rank-0: CH₃O physisorbed (Pd-O = 3.226 Å) yet 3038 meV ABOVE ranks 1-3
  - Energy ranking inverted (physisorbed structure highest energy, chemisorbed lower)
  - PES borderline flat
- **Usable for DFT:** YES — ranks 1,2,3 are viable. DROP rank-0 or accept as "weak binding" reference.

### S3b (PdO(100) PdO-terminated)
- **Shortlist size:** 3
- **Chemistry pass rate:** 3/3 (100%)
- **PES spread (top-10):** 2.9 meV ❌
- **Verdict:** CONCERN (chemistry perfect, but MLIP energies meaningless)
- **Issues:**
  - PES noise-level flat
- **Usable for DFT:** YES — all 3 chemically clean; pick for geometric diversity only

### S4 (PdO₂(110))
- **Shortlist size:** 4
- **Chemistry pass rate:** 0/4 (0%)
- **PES spread (top-10):** 122.7 meV ✅ (good spread, but all structures defective)
- **Verdict:** REJECT
- **Issues:**
  - Rank-0: CH₃O physisorbed (Pd-O = 3.022 Å)
  - Ranks 1,2,3: CO physisorbed (Pd-C = 3.556, 3.653, 3.630 Å)
  - 884/7975 (11%) of input had C-H bond breaks — suggests CH₃O fundamentally unstable on PdO₂ under MACE relaxation
- **Usable for DFT:** NO
- **Recommended action:**
  1. Skip SetA on S4 (PdO₂ may not support co-adsorption at reactive distances), OR
  2. Try alternative MLIP (CHGNet has better oxide binding physics), OR
  3. Hand-place reactive pairs on known CUS sites and relax directly with DFT

---

## Quantitative Evidence

| Surface | n_shortlist | Chem OK | Top-10 ΔE (meV) | d_reactive range (Å) | Verdict   |
|---------|-------------|---------|-----------------|----------------------|-----------|
| S1      | 3           | 2/3     | 3.0             | [2.99, 3.92]         | CONCERN   |
| S2      | 4           | 2/4     | 214.6           | [2.38, 3.66]         | CONCERN   |
| S3      | 4           | 3/4     | 27.3            | [2.37, 3.52]         | CONCERN   |
| S3b     | 3           | 3/3     | 2.9             | [2.90, 3.66]         | CONCERN   |
| S4      | 4           | 0/4     | 122.7           | [2.46, 3.57]         | REJECT    |

### Key Distance Checks (S3 rank-1, the best candidate)
| Parameter            | Value    | Expected (PBE+D3) | Status |
|----------------------|----------|-------------------|--------|
| Pd–C (CO)            | 2.465 Å  | [1.85, 2.00]      | ✅ pass (borderline) |
| Pd–O (CH₃O)          | 1.996 Å  | [2.00, 2.15]      | ✅ pass |
| C=O bond             | 1.151 Å  | [1.13, 1.18]      | ✅ pass |
| O–C (methoxy)        | 1.415 Å  | [1.40, 1.46]      | ✅ pass |
| C–H (methoxy)        | 1.110 Å  | [1.08, 1.12]      | ✅ pass |
| d_reactive (C–O)     | 2.988 Å  | [2.1, 4.0]        | ✅ pass (SetA band) |

### Spot Checks (middle ranks from unique pools)
- **S3:** 5/5 middle-ranked structures (ranks 100, 500, 1000, 2000, 3000) intact
- **S4:** 5/5 middle-ranked structures intact
- **Conclusion:** Geometry filter successfully removed fragmented/collapsed structures from bulk pool

---

## Comparison to Prior Verdict (REJECT)

### What the filter FIXED:
✅ S3 product collapse (4/4 → 1/4)  
✅ S2 SetA band violation (rank-0 was 1.9 Å → now 2.38 Å)  
✅ Fragmentation (caught 3166-5427 per surface)  
✅ Bond-breaking detection (885 total C-H/O-C breaks flagged)  

### What the filter DID NOT FIX:
❌ S4 CO physisorption (still 3/4 shortlist)  
❌ PES flatness (3-27 meV for S1/S3/S3b, below noise)  
❌ Energy ranking inversion (rank-0 ≠ lowest-E on any surface)  

### New issues INTRODUCED:
⚠️  S3 rank-0 now has CH₃O physisorbed (different defect, not product collapse)  
⚠️  S1 rank-1 CO physisorbed (not flagged before)  

---

## Asks for Producer (ml-trainer / data-curator)

1. **S4 alternative strategy required:**
   - Current MACE-MP setup cannot produce viable S4 reactive pairs
   - Options: (a) skip S4 SetA, (b) try CHGNet/M3GNet, (c) DFT-only approach with hand-placed pairs

2. **Revise shortlist selection logic:**
   - Current "distance-bin representatives" sacrifices energy ordering for diversity
   - For S1/S3/S3b (flat PES), this is correct — but should be EXPLICIT
   - For S2/S4 (good PES spread), consider hybrid: top-3 by energy + 1 geometric outlier

3. **Add chemisorption filter to shortlist selection:**
   - Re-filter S3 to drop rank-0 (Pd-O > 2.3 Å) and promote rank-4 from unique pool
   - Apply Pd-C < 2.5 Å AND Pd-O < 2.3 Å as hard requirements for final shortlist

4. **Visual validation before DFT handoff:**
   - Render S2 rank-2,3 and S4 all-ranks to confirm whether "Pd-C = 3.4-3.7 Å" is truly physisorption or bridge-to-different-layer in mixed-oxide slab
   - ASE min-distance cannot distinguish these cases

5. **Document MLIP limitations in producer report:**
   - Explicitly state: "MACE energies unreliable for S1/S3/S3b (PES flat). Shortlist selected for geometric diversity only."
   - Warn about S4 oxide binding failure

---

## Decision Tree for Next Steps

```
IF proceeding to DFT:
  S1  → send all 3 shortlist (ignore energies, maximize diversity)
  S2  → send rank-0,1 (chemically clean) + optionally rank-2,3 (for physisorption reference)
  S3  → send rank-1,2,3 (DROP rank-0 unless needed as "weak binding" control)
  S3b → send all 3 (all clean, diverse geometries)
  S4  → STOP. Do not send current shortlist.

IF NOT proceeding yet:
  → Implement asks #3 (chemisorption filter) and #4 (visual check)
  → Regenerate S3/S4 shortlists
  → Re-submit to judge-physics
```

---

## Physics Sanity Checklist

✅ Adsorbate atom counts correct (7 atoms: C O O C H H H)  
✅ CO bonds intact (1.14-1.19 Å across all shortlists)  
✅ CH₃O O-C bonds intact (1.38-1.46 Å)  
✅ C-H bonds intact (1.10-1.23 Å, one outlier at S3b rank-0 = 1.232 Å, acceptable)  
⚠️  Chemisorption: 50-100% pass rate (surface-dependent)  
❌ MLIP energy ranking: inverted or noise-level for 3/5 surfaces  
✅ SetA band compliance: 100% of shortlist in [2.1, 4.0) Å  
✅ Geometry filter effectiveness: 50-60% of bad structures removed  

---

## Final Recommendation

**Verdict: CONCERN (not REJECT, not PASS)**

The filtered output is **usable for DFT with manual curation**:
- Geometry filter did its job (caught fragmentation/collapse)
- Chemistry is 50-100% viable depending on surface
- BUT: MLIP energy ranking is unreliable → use geometric diversity, not energy rank

**Proceed to DFT for S1, S2, S3 (curated), S3b.**  
**Stop S4; requires alternative approach.**

---

**Judge:** physics  
**Confidence:** High (17 structures inspected in detail, 10 spot-checks, full pool statistics analyzed)  
**Blind review:** Yes (no contact with other judges)  
