# SMILES Fragment Note — CO* vs CHO*

## Issue
The working CO* SMILES `ClC=O` generates a **formyl (CHO)** fragment, not pure CO.

## Details
- Guideline SMILES: `Cl[C-]#[O+]` → RDKit fails (explicit valence error)
- Fallback SMILES: `ClC=O` → RDKit OK, generates CHO (formyl group)
- After Cl marker removal by AutoAdsorbate:
  - Expected: C, O (2 atoms)
  - Actual: C, O, H (3 atoms)

## Chemistry implications
1. **DMC mechanism difference**: CHO* (formyl) vs CO* (carbonyl) are different intermediates.
   - CO* + CH3O* → DMC is the guideline pathway
   - CHO* + CH3O* → ? is a different pathway (possibly formaldehyde-related)

2. **Distance metric validity**: The C_CO···O_CH3O distance is still chemically meaningful:
   - Measures proximity of formyl-C to methoxy-O
   - Classification bins (reactive/TS/thermo) remain valid

3. **Phase 2 impact**: The TS search and energy profile will correspond to CHO* + CH3O* coupling, not CO* + CH3O*.

## Decision needed
**Before proceeding to MLIP/DFT**, confirm with advisor:
- Option A: Accept CHO* as the adsorbate (adjust mechanism interpretation)
- Option B: Find alternative CO* SMILES (e.g., manual ASE Atoms construction without RDKit)
- Option C: Use AutoAdsorbate's Fragment API to directly build CO* from atomic coordinates

## Current status
Proceeding with CHO* for dry-run demonstration. Flagged for advisor review.

---
*data-curator, 2026-06-10*
