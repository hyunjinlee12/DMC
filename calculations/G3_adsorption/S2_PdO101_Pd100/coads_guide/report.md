# Co-adsorption Generation Report: S2_PdO101_Pd100

## Summary

Generated co-adsorption candidates following guideline specification (§P1-C, §9).

### Input
- CO candidates: 166
- CH3O candidates: 498

### Processing Pipeline

**Step 1: Enumerate all pairs**
- Total pair combinations: 82668

**Step 2: Primary filter (site-pair distance)**
- Primary pool (≤ 4.5 Å): 33545
- Thermo pool (≥ 5.0 Å): 40946
- Ambiguous (4.5-5.0 Å, excluded): 8177

**Step 3: max_per_pair=2**
- Primary pool after cap: 23948 (from 33545)
- Thermo pool after cap: 29147 (from 40946)

**Step 4: Steric clash filter**
- Primary pool passed: 19013
- Primary pool rejected: 4935
- Thermo pool passed: 29147
- Thermo pool rejected: 0

**Step 5: Classification (guideline cutoffs)**

| Set | Cutoff | Count | Distance (Å) |
|-----|--------|-------|--------------|
| Side-path | < 1.7 Å | 183 | 1.60 - 1.70 (mean 1.65) |
| Set TS | 1.7 - 2.3 Å | 1625 | 1.70 - 2.30 (mean 2.02) |
| Set A | 2.1 - 4.0 Å | 11598 | 2.10 - 4.00 (mean 3.25) |
| Set B | ≥ 5.0 Å | 29147 | 5.00 - 8.74 (mean 6.14) |
| Rejected | steric clash | 4935 | N/A |

**Note:** Set A and Set TS overlap in range 2.1-2.3 Å per guideline specification.
Number of structures in overlap region: 594

### Target Compliance

Guideline target: 150-500 candidates per surface (Phase 1 co-adsorption)

- **Set A (reactive)**: 11598 ✓ outside target
- **Set B (thermo)**: 29147 ✓ within target
- **Combined relevant (A+TS+B-overlap)**: 41776

### Outputs

- `SetA.traj`: 11598 reactive pair candidates (2.1-4.0 Å)
- `SetTS.traj`: 1625 TS guess candidates (1.7-2.3 Å)
- `SetB.traj`: 29147 thermodynamic reference candidates (≥5.0 Å)
- `side.traj`: 183 side-path candidates (<1.7 Å)
- `rejected.traj`: 4935 steric clash rejects
- `summary.json`: full statistics
- `distance_hist.png`: C_CO-O_CH3O distance distribution
- `SetA_grid.png`: 25 Set A representatives
- `SetB_grid.png`: 25 Set B representatives

### Cutoffs Applied (Guideline §9)

| Parameter | Cutoff |
|-----------|--------|
| site-pair primary | ≤ 4.5 Å |
| reactive C_CO···O_OCH3 | 2.1 - 4.0 Å |
| TS guess C···O | 1.7 - 2.3 Å |
| thermodynamic reference | ≥ 5.0 Å |
| steric heavy-heavy reject | < 1.6 Å |
| steric H-heavy reject | < 1.1 Å |

---
Generated: S2_PdO101_Pd100 co-adsorption candidates (guideline specification)
