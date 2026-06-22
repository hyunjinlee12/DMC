# PdO Slab Structure Audit — S3 vs S3b

## Bulk reference (PdO tetragonal P4₂/mmc)
- a = 3.054 Å (exp 3.043, +0.35%)
- c = 5.406 Å (exp 5.336, +1.31%)
- Pd-O bond mean: 2.039 Å (n=4)
- Each Pd 4-coord to O, each O 4-coord to Pd (square-planar PdO)

## Slabs (T1.9 validation re-audit)

### S3
- atoms: 128 (Pd:64, O:64, O/Pd = 1.00)
- layers: 8, thickness 10.62 Å, vacuum 20.07 Å
- TOP layer (z=12.07, rumpling=0.14 Å): O:16, Pd:8
- 2nd layer (z=10.75): Pd:8
- Pd-O bonds: n=240, mean **2.034 Å** (bulk reference 2.039 Å), range [1.966, 2.088]
- Top-layer Pd coord to O: [4, 4, 4, 4, 4, 4, 4, 4] (bulk = 4)
- Charge-weighted dipole-z proxy: -97.73

### S3b
- atoms: 104 (Pd:56, O:48, O/Pd = 0.86)
- layers: 7, thickness 9.28 Å, vacuum 19.97 Å
- TOP layer (z=10.78, rumpling=0.00 Å): Pd:8
- 2nd layer (z=9.14): O:16, Pd:8
- Pd-O bonds: n=192, mean **2.046 Å** (bulk reference 2.039 Å), range [2.027, 2.126]
- Top-layer Pd coord to O: [2, 2, 2, 2, 2, 2, 2, 2] (bulk = 4)
- Charge-weighted dipole-z proxy: +1.01

## Interpretation

- **S3 (claimed O-term)**: top layer is `8 Pd + 16 O` (mixed, O-rich). The 16 O is double the 8 Pd → every top-Pd is bridged by 2 extra O above. This makes the surface effectively O-terminated (no bare Pd on top — O atoms protrude above the Pd plane). Naming "O-term" is correct in spirit (O-rich top with O protrusion) but technically it's a stoichiometric PdO(100) slab with the O-on-top motif.
- **S3b (claimed Pd-term)**: top layer is `8 Pd` only (pure Pd). One extra layer of O is removed compared to S3, making it Pd-rich (O/Pd = 0.857). The Pd on top is undercoordinated (3-fold, vs bulk 4-fold).
- **Pd-O bond preservation**: both slabs have Pd-O bond lengths matching bulk (~2.04 Å mean), within PBE+D3 typical error. No bond reconstruction issues.
- **Dipole**: S3 has larger charge-weighted dipole proxy because of the asymmetric Pd-bottom / O-top composition (charge imbalance). S3b is more symmetric in the z-direction (uniform Pd/O distribution).
- **Which is more physical for DMC chemistry?** 
  - S3 (O-rich): represents oxygen-saturated PdO conditions (high O coverage limit).   Likely Case C/D candidate (CO weakly binds, methoxy possible).
  - S3b (Pd-rich): represents partially reduced PdO. Top-Pd more reactive (3-coord).   Likely Case A/B (DMC favorable). Closer to Pd/PdO interface chemistry.

## Validity for project

Both are physically reasonable terminations of PdO(100). Their inclusion together gives the project a **stoichiometry axis** (O-rich vs Pd-rich) which adds dimensionality to the descriptor map beyond the single-termination Shi 2024 comparison.
