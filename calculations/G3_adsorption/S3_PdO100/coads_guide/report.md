# Co-adsorption Generation Report: S3_PdO100

## Summary

Generated co-adsorption candidates following guideline specification (§P1-C, §9).

### Input
- CO candidates: 135
- CH3O candidates: 405

### Processing Pipeline

**Step 1: Enumerate all pairs**
- Total pair combinations: 54675

**Step 2: Primary filter (site-pair distance)**
- Primary pool (≤ 4.5 Å): 25412
- Thermo pool (≥ 5.0 Å): 22592
- Ambiguous (4.5-5.0 Å, excluded): 6671

**Step 3: max_per_pair=2**
- Primary pool after cap: 17887 (from 25412)
- Thermo pool after cap: 16110 (from 22592)

**Step 4: Steric clash filter**
- Primary pool passed: 13889
- Primary pool rejected: 3998
- Thermo pool passed: 16110
- Thermo pool rejected: 0

**Step 5: Classification (guideline cutoffs)**

| Set | Cutoff | Count | Distance (Å) |
|-----|--------|-------|--------------|
| Side-path | < 1.7 Å | 177 | 1.60 - 1.70 (mean 1.65) |
| Set TS | 1.7 - 2.3 Å | 1441 | 1.70 - 2.30 (mean 2.01) |
| Set A | 2.1 - 4.0 Å | 8281 | 2.10 - 4.00 (mean 3.23) |
| Set B | ≥ 5.0 Å | 16110 | 5.00 - 8.16 (mean 5.97) |
| Rejected | steric clash | 3998 | N/A |

**Note:** Set A and Set TS overlap in range 2.1-2.3 Å per guideline specification.
Number of structures in overlap region: 489

### Target Compliance

Guideline target: 150-500 candidates per surface (Phase 1 co-adsorption)

- **Set A (reactive)**: 8281 ✓ outside target
- **Set B (thermo)**: 16110 ✓ within target
- **Combined relevant (A+TS+B-overlap)**: 25343

### Outputs

- `SetA.traj`: 8281 reactive pair candidates (2.1-4.0 Å)
- `SetTS.traj`: 1441 TS guess candidates (1.7-2.3 Å)
- `SetB.traj`: 16110 thermodynamic reference candidates (≥5.0 Å)
- `side.traj`: 177 side-path candidates (<1.7 Å)
- `rejected.traj`: 3998 steric clash rejects
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
Generated: S3_PdO100 co-adsorption candidates (guideline specification)
