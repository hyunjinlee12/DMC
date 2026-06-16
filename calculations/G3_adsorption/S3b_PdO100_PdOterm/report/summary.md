# S3b_PdO100_PdOterm G3 Adsorption Sampling Report

Generated: 2026-06-15T15:46:40

## 1. Site Map

- Total sites detected: 87
- Slab atoms: 104
- Formula: O48Pd56

### By site type:
```
{
  "5-fold": 4,
  "3-fold-hollow": 23,
  "bridge": 44,
  "atop": 16
}
```

## 2. CO* Candidates

- Total: 87
- SMILES: `Cl[C]=O`
- Fragment atoms: {'C': 1, 'O': 1}
- Note: guideline Cl[C-]#[O+] failed RDKit valence; corrected to Cl[C]=O (brackets prevent implicit H)

## 3. CH3O* Candidates

- Total: 261
- SMILES: `ClOC`
- Fragment atoms: {'O': 1, 'C': 1, 'H': 3}
- Generation params: {'to_initialize': 20, 'sample_rotation': True, 'conformers_per_site_cap': 3}

## 4. Co-adsorption

### Classification (guideline-compliant cutoffs):

- **Set A (reactive, 2.1–4.0 Å)**: 100
- **Set B (thermo, ≥5.0 Å)**: 162
- **Set TS (1.7–2.1 Å)**: 15
- **side-path (<1.7 Å)**: 2
- **rejected (steric)**: 36

### Distance statistics:

- Min: 0.62 Å
- Max: 7.78 Å
- Mean: 4.44 Å

### Cutoffs applied:
```
{
  "side": 1.7,
  "TS_low": 1.7,
  "TS_high": 2.1,
  "A_low": 2.1,
  "A_high": 4.0,
  "amb_high": 5.0,
  "steric_heavy": 1.6,
  "steric_H": 1.1
}
```

## 5. Output Files

- Site map: `site_map/sites.json`, `site_map/overlay.png`
- CO*: `CO/candidates.traj`, `CO/summary.json`, `CO/grid.png`
- CH3O*: `CH3O/candidates.traj`, `CH3O/summary.json`, `CH3O/grid.png`
- Co-ads: `coads/SetA.traj`, `coads/SetB.traj`, `coads/SetTS.traj`, `coads/side.traj`, `coads/rejected.traj`
- Co-ads plots: `coads/distance_hist.png`, `coads/SetA_grid.png`, `coads/SetB_grid.png`
- Report: `report/summary.md`

## 6. Notes

- Pipeline: AutoAdsorbate mode='all', overlap_thr=1.25
- CO* SMILES corrected from guideline to avoid RDKit valence error
- Steric rejection: heavy-heavy < 1.6 Å OR H-heavy < 1.1 Å
- Ambiguous range (4.0–5.0 Å) excluded from all sets

**Checkpoint stop: NO MLIP, NO DFT performed.**
