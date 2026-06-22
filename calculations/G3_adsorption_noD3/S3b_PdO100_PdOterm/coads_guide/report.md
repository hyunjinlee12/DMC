# Co-adsorption Generation Report: S3b_PdO100_PdOterm

## Summary

Generated co-adsorption candidates following guideline specification (§P1-C, §9).

### Input
- CO candidates: 87
- CH3O candidates: 261

### Processing Pipeline

**Step 1: Enumerate all pairs**
- Total pair combinations: 22707

**Step 2: Primary filter (site-pair distance)**
- Primary pool (≤ 4.5 Å): 10342
- Thermo pool (≥ 5.0 Å): 9664
- Ambiguous (4.5-5.0 Å, excluded): 2701

**Step 3: max_per_pair=2**
- Primary pool after cap: 7500 (from 10342)
- Thermo pool after cap: 7157 (from 9664)

**Step 4: Steric clash filter**
- Primary pool passed: 6199
- Primary pool rejected: 1301
- Thermo pool passed: 7157
- Thermo pool rejected: 0

**Step 5: Classification (guideline cutoffs)**

| Set | Cutoff | Count | Distance (Å) |
|-----|--------|-------|--------------|
| Side-path | < 1.7 Å | 78 | 1.60 - 1.70 (mean 1.65) |
| Set TS | 1.7 - 2.3 Å | 469 | 1.70 - 2.30 (mean 2.02) |
| Set A | 2.1 - 4.0 Å | 3764 | 2.10 - 4.00 (mean 3.27) |
| Set B | ≥ 5.0 Å | 7157 | 5.00 - 8.14 (mean 5.91) |
| Rejected | steric clash | 1301 | N/A |

**Note:** Set A and Set TS overlap in range 2.1-2.3 Å per guideline specification.
Number of structures in overlap region: 173

### Target Compliance

Guideline target: 150-500 candidates per surface (Phase 1 co-adsorption)

- **Set A (reactive)**: 3764 ✓ outside target
- **Set B (thermo)**: 7157 ✓ within target
- **Combined relevant (A+TS+B-overlap)**: 11217

### Outputs

- `SetA.traj`: 3764 reactive pair candidates (2.1-4.0 Å)
- `SetTS.traj`: 469 TS guess candidates (1.7-2.3 Å)
- `SetB.traj`: 7157 thermodynamic reference candidates (≥5.0 Å)
- `side.traj`: 78 side-path candidates (<1.7 Å)
- `rejected.traj`: 1301 steric clash rejects
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
Generated: S3b_PdO100_PdOterm co-adsorption candidates (guideline specification)
