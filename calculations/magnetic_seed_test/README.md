# Magnetic seed test — resolve MAGMOM initialization uncertainty

Reviewer 2026-07-17 round 4:
Current T1.16 candidates converged to significant magnetization (~21 μ_B for
S1/CO, ~11 μ_B for S3/CH3O) using VASP default MAGMOM=1.0/atom initialization.
Layer-wise decomposition shows moments distributed throughout the slab, not
surface-localized — suggests **metastable FM state from default seed**, not
physical surface Stoner magnetism.

This test computes STATIC single-point energies with different MAGMOM seeds on
the SAME geometry. Compare final `energy(sigma→0)` and final total magnetization.

## Seeds (11 dirs total)

| dir | ISPIN | NUPDOWN | MAGMOM | description |
|---|---|---|---|---|
| S1_clean__ispin1  | 1 | – | (none) | closed-shell, no spin polarization |
| S1_clean__lowmag  | 2 | – | 80*0.2 | small FM seed |
| S1_clean__highmag | 2 | – | 80*1.0 | default VASP seed reproduction |
| S1_CO_idx3__ispin1  | 1 | – | – | ditto with CO adsorbed |
| S1_CO_idx3__lowmag  | 2 | – | 1*0 1*0 80*0.2 | |
| S1_CO_idx3__highmag | 2 | – | 1*0 1*0 80*1.0 | |
| S3_clean__nonmag | 1 | – | – | S3 clean nonmagnetic |
| S3_clean__fm     | 2 | – | 64*0 64*0.5 | ferromagnetic seed on Pd |
| S3_clean__afm    | 2 | – | Pd atoms split ± by fractional-x | antiferromagnetic seed |
| S3_CH3O_idx315__doublet       | 2 | 1 | substrate 0, ads O = 1 | doublet ground state |
| S3_CH3O_idx315__unconstrained | 2 | – | Pd 0.3, ads O 1 | unconstrained radical |

## Common

- IBRION = -1, NSW = 0  (STATIC — no ionic movement)
- ISTART = 0, ICHARG = 2  (no WAVECAR/CHGCAR reuse — seed independence)
- KSPACING = 0.25 (matches T1.17)
- LDIPOL/IDIPOL = 3 (matches T1.17 slab)
- ENCUT/PREC/LASPH/IVDW/EDIFF/EDIFFG all match T1.17
- Same POSCAR geometry as the parent T1.16 CONTCAR (adsorbate) or G2 CONTCAR (clean)

## After submission

For each dir, extract from OUTCAR:
- `energy(sigma->0)` (grep last)
- Total magnetization (`grep "number of electron" OUTCAR | tail -1`)
- Per-atom moment (mag(x) table)
- Layer-wise Pd moment sum
- Electronic convergence status

## Decision rules

- All seeds converge to **same E and moment** → magnetic solution is stable → adopt.
- Seeds converge to **different moments but one has clearly lower E** → adopt the lowest-E state.
- Multiple states within **5–10 meV** → magnetic uncertainty; run NUPDOWN fixed-spin scan.
- Multiple magnetization basins with significant E span → propagate as an uncertainty in downstream ΔG.

Once T1.17 MAGMOM policy is decided from this test, regenerate T1.17 dirs with
that policy before mass submission.
