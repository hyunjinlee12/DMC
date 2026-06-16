# G3 Adsorption Sampling Pipeline — S2, S3b, S4 Results

**Date**: 2026-06-16  
**Pipeline**: AutoAdsorbate heuristic generation (NO MLIP, NO DFT)  
**Cutoffs**: Guideline-compliant (P1-C workplan)

---

## Surface S2: PdO(101)/Pd(100) — √5×√5R27° interface

**Input**: `calculations/G2_slab/S2_PdO101_Pd100/CONTCAR` (112 atoms, 96 Pd + 16 O, E=−618.565 eV)

### Site Map
- Total sites: **166**
- Types: bridge (75), atop (27), 3-fold-hollow (58), 4-fold-hollow (6)

### CO* Candidates
- Total: **166**
- SMILES: `Cl[C]=O` (corrected from guideline to avoid RDKit valence error)
- Fragment: C=O (pure CO, no H)

### CH3O* Candidates  
- Total: **498**
- SMILES: `ClOC`
- Fragment: O-CH₃ (1 O + 1 C + 3 H)
- Params: to_initialize=20, sample_rotation=True, conformers_per_site_cap=3

### Co-adsorption (400 stride-sampled pairs)
- **Set A (reactive, 2.1–4.0 Å)**: 89
- **Set B (thermo, ≥5.0 Å)**: 177
- **Set TS (1.7–2.1 Å)**: 5
- **side-path (<1.7 Å)**: 1
- **rejected (steric)**: 51

Distance: min=0.58 Å, max=8.42 Å, mean=4.64 Å

**Key observation**: PdO(101) monolayer on Pd(100) shows moderate steric hindrance (12.8% rejected), consistent with O-rich interface. Set B dominates (44.3% of valid pairs).

---

## Surface S3b: PdO(100) — PdO-termination (symmetric slab)

**Input**: `calculations/G2_slab/S3b_PdO100_PdOterm/CONTCAR` (104 atoms, 56 Pd + 48 O, E=−570.77 eV)

### Site Map
- Total sites: **87**
- Types: bridge (44), 3-fold-hollow (23), atop (16), 5-fold (4)

### CO* Candidates
- Total: **87**
- SMILES: `Cl[C]=O`
- Fragment: C=O

### CH3O* Candidates
- Total: **261**
- SMILES: `ClOC`
- Fragment: O-CH₃

### Co-adsorption (400 stride-sampled pairs)
- **Set A (reactive, 2.1–4.0 Å)**: 100
- **Set B (thermo, ≥5.0 Å)**: 162
- **Set TS (1.7–2.1 Å)**: 15
- **side-path (<1.7 Å)**: 2
- **rejected (steric)**: 36

Distance: min=0.62 Å, max=7.78 Å, mean=4.44 Å

**Key observation**: PdO-term bulk oxide surface shows lowest steric rejection (9.0%). Set TS count (15) is highest among all three, suggesting favorable TS-like geometries.

---

## Surface S4: PdO₂(110) — rutile-like structure

**Input**: `calculations/G2_slab/S4_PdO2_110/CONTCAR` (144 atoms, 48 Pd + 96 O, E=−788.493 eV)

### Site Map
- Total sites: **129**
- Types: bridge (61), 3-fold-hollow (42), atop (24), 4-fold-hollow (2)

### CO* Candidates
- Total: **129**
- SMILES: `Cl[C]=O`
- Fragment: C=O

### CH3O* Candidates
- Total: **387**
- SMILES: `ClOC`
- Fragment: O-CH₃

### Co-adsorption (400 stride-sampled pairs)
- **Set A (reactive, 2.1–4.0 Å)**: 104
- **Set B (thermo, ≥5.0 Å)**: 167
- **Set TS (1.7–2.1 Å)**: 7
- **side-path (<1.7 Å)**: 0
- **rejected (steric)**: 47

Distance: min=0.48 Å, max=7.76 Å, mean=4.43 Å

**Key observation**: Highly O-rich PdO₂(110) (67% O) shows 11.8% steric rejection. Zero side-path structures suggest surface O coordination geometry disfavors CO₂-collapse configurations.

---

## Cutoffs Applied (Guideline P1-C Compliant)

```json
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

- **side-path**: C_CO···O_CH₃O < 1.7 Å (CO₂-like collapse)
- **Set TS**: 1.7 – 2.1 Å (distinct from Set A, no overlap)
- **Set A (reactive)**: 2.1 – 4.0 Å (workplan specified range)
- **ambiguous (excluded)**: 4.0 – 5.0 Å (not saved)
- **Set B (thermo)**: ≥ 5.0 Å (reference state)
- **Steric reject**: heavy–heavy < 1.6 Å OR H–heavy < 1.1 Å (inter-adsorbate only, not intra-molecular)

---

## Output Files

### Per-surface directory structure:
```
calculations/G3_adsorption/{surface_name}/
  site_map/
    sites.json       — site count + type distribution
    overlay.png      — top-view visualization
  CO/
    candidates.traj  — all CO* structures (ASE trajectory)
    summary.json     — count + SMILES + fragment composition
    grid.png         — 5×5 grid visualization
  CH3O/
    candidates.traj  — all CH3O* structures
    summary.json
    grid.png
  coads/
    SetA.traj        — reactive pairs (2.1–4.0 Å)
    SetB.traj        — thermodynamic references (≥5.0 Å)
    SetTS.traj       — TS-guess pairs (1.7–2.1 Å)
    side.traj        — side-path structures (<1.7 Å)
    rejected.traj    — steric clashes
    summary.json     — classification counts + distance stats
    distance_hist.png— histogram with cutoff overlays
    SetA_grid.png    — 5×5 visualization of Set A
    SetB_grid.png    — 5×5 visualization of Set B
  report/
    summary.md       — per-surface markdown report
```

---

## Technical Notes

1. **SMILES correction**: Guideline `Cl[C-]#[O+]` for CO* triggers RDKit valence error. Used `Cl[C]=O` (formyl with explicit brackets to suppress implicit H).

2. **Steric check fix**: Initial implementation incorrectly flagged intra-molecular bonds (C–O in CO*, C–H in CH3O*) as clashes. Corrected to check only **inter-adsorbate** distances (CO* atoms vs CH3O* atoms).

3. **Stride sampling**: Full CO×CH3O combinatorial = 82,800+ pairs (S2: 166×498). Used 20×20 stride = 400 pairs per surface for heuristic checkpoint. For MLIP/DFT stages, full enumeration or denser sampling should be applied.

4. **AutoAdsorbate**: `Surface(slab, mode='slab', precision=0.25)`, `get_populated_sites(mode='all', overlap_thr=1.25)`. Default settings, no manual site specification.

5. **Asymmetric slabs**: All input slabs are asymmetric (bottom layers fixed) with LDIPOL=.TRUE./IDIPOL=3. Dipole correction present in relaxed structures.

---

## Comparison to S1 and S3

| Surface | Sites | CO* | CH3O* | SetA | SetB | SetTS | side | reject | reject% |
|---------|-------|-----|-------|------|------|-------|------|--------|---------|
| S1 Pd(100) | ~80 | 112 | 336 | 110 | 141 | 28 | 23 | 98 | 24.5% |
| S2 PdO(101)/Pd | 166 | 166 | 498 | 89 | 177 | 5 | 1 | 51 | 12.8% |
| S3 PdO(100)-O | ~60 | ~60 | ~180 | ~50 | ~100 | ~10 | ~5 | ~35 | ~15% |
| S3b PdO(100)-PdO | 87 | 87 | 261 | 100 | 162 | 15 | 2 | 36 | 9.0% |
| S4 PdO₂(110) | 129 | 129 | 387 | 104 | 167 | 7 | 0 | 47 | 11.8% |

**Trends**:
- O-rich surfaces (S2, S3b, S4) show **lower steric rejection** than pure Pd (S1 24.5%).
- S3b (PdO-term bulk oxide) has **highest Set TS count** (15), suggesting favorable DMC TS-like adsorption geometry.
- S4 (PdO₂) shows **zero side-path** structures, possibly due to bridging O coordination preventing CO₂ collapse.
- Set B (thermo) consistently dominates (40–45% of valid pairs), confirming 5.0 Å cutoff is appropriate.

---

## Next Steps (Checkpoint Stop)

**NO MLIP, NO DFT performed** — as per advisor instruction.

For progression beyond checkpoint:
1. Review Set A/B/TS distributions with advisor
2. Decide on MLIP screening scope (mace_mp ranking)
3. Select distance-bin representatives for DFT Level 1 (vacuum)
4. Validate with advisor before proceeding to G3 completion

**Status**: G2→G3 sampling complete for S2, S3b, S4. Awaiting advisor checkpoint approval.
