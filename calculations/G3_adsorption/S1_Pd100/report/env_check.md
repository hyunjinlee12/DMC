# Environment Check - S1 Pd(100)

- AutoAdsorbate: 0.2.5
- ASE: OK
- RDKit: OK

## SMILES
- CO* (guideline `Cl[C-]#[O+]` FAILED in RDKit)
- CO* (working): `ClC=O`
- CH3O*: `ClOC`

## Deviation from guideline
The guideline SMILES `Cl[C-]#[O+]` for CO* causes explicit valence errors in RDKit.
Used `ClC=O` instead (formyl group, Cl marker removed yields C=O adsorbate).
