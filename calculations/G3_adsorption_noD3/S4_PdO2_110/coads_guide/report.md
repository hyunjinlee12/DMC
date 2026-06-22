# Co-adsorption Generation Report: S4_PdO2_110

## Summary

Generated co-adsorption candidates following guideline specification (§P1-C, §9).

### Input
- CO candidates: 129
- CH3O candidates: 387

### Processing Pipeline

**Step 1: Enumerate all pairs**
- Total pair combinations: 49923

**Step 2: Primary filter (site-pair distance)**
- Primary pool (≤ 4.5 Å): 24654
- Thermo pool (≥ 5.0 Å): 20051
- Ambiguous (4.5-5.0 Å, excluded): 5218

**Step 3: max_per_pair=2**
- Primary pool after cap: 17439 (from 24654)
- Thermo pool after cap: 14270 (from 20051)

**Step 4: Steric clash filter**
- Primary pool passed: 13027
- Primary pool rejected: 4412
- Thermo pool passed: 14270
- Thermo pool rejected: 0

**Step 5: Classification (guideline cutoffs)**

| Set | Cutoff | Count | Distance (Å) |
|-----|--------|-------|--------------|
| Side-path | < 1.7 Å | 120 | 1.60 - 1.70 (mean 1.65) |
| Set TS | 1.7 - 2.3 Å | 1182 | 1.70 - 2.30 (mean 2.02) |
| Set A | 2.1 - 4.0 Å | 7975 | 2.10 - 4.00 (mean 3.24) |
| Set B | ≥ 5.0 Å | 14270 | 5.00 - 8.00 (mean 6.05) |
| Rejected | steric clash | 4412 | N/A |

**Note:** Set A and Set TS overlap in range 2.1-2.3 Å per guideline specification.
Number of structures in overlap region: 442

### Target Compliance

Guideline target: 150-500 candidates per surface (Phase 1 co-adsorption)

- **Set A (reactive)**: 7975 ✓ outside target
- **Set B (thermo)**: 14270 ✓ within target
- **Combined relevant (A+TS+B-overlap)**: 22985

### Outputs

- `SetA.traj`: 7975 reactive pair candidates (2.1-4.0 Å)
- `SetTS.traj`: 1182 TS guess candidates (1.7-2.3 Å)
- `SetB.traj`: 14270 thermodynamic reference candidates (≥5.0 Å)
- `side.traj`: 120 side-path candidates (<1.7 Å)
- `rejected.traj`: 4412 steric clash rejects
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
Generated: S4_PdO2_110 co-adsorption candidates (guideline specification)
