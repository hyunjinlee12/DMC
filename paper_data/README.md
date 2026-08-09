# Paper-ready raw data

Consolidated from repo state; regenerable via `scripts/build_paper_data.py`.
All energies in eV. All distances in Å. All angles in degrees.

## Files

- **01_bulk_data.csv** — G1 bulk relaxation results for Pd, PdO, PdO₂ (and one PdO₂ w/o D3 for comparison). Includes lattice a/b/c, angles, volume, total energy, per-atom energy, experimental lattice references, and deviation from experiment.
- **02_slab_data.csv** — G2 clean-slab results for the 5 modeled surfaces (S1 Pd(100), S2 PdO(101)/Pd(100), S3 PdO(100) O-rich, S3b PdO(100) Pd-terminated, S4 PdO₂(110)). Includes atom counts, top-layer composition, rumpling, min/avg Pd–O bond length in the top region, total energy, max force on free atoms.
- **03_mace_references.csv** — reference energies used to compute E_bind_MACE_D3 in downstream tables. Same MACE-MH+D3(BJ)+cueq calculator as the pool relaxations.
- **04_mlip_pool_singles.csv** — every MLIP-relaxed unique CO*/CH₃O* candidate across the 5 surfaces (~1200 rows). Columns: sid, ads, idx, E_MACE_D3, E_bind_MACE_D3, dE_rel_meV, converged, n_steps, d_min_ads_sub, site_type (raw MLIP label), fingerprint.
- **05_mlip_pool_coads.csv** — every MLIP-relaxed unique coadsorption candidate (Set A) for S1/S2/S3/S3b (~14600 rows). Adds d_reactive (C_CO ↔ O_CH3O, Å, MIC) and distance_bin classification.
- **06_dft_shortlist.csv** — 86 candidates selected for DFT L2 (70 from v4 selector + 16 from v5 additions). Includes site/region/distance-bin labels, priority (v4-baseline / MUST / RECOMMENDED / OPTIONAL / MUST-diagnostic), selection reason, DFT hypothesis.
- **07_dft_results.csv** — 86 rows in one-to-one correspondence with 06; DFT_status is one of DONE / RUNNING_OR_UNCONVERGED / PENDING. Once a given candidate finishes, `scripts/build_paper_data.py` re-run automatically fills E_DFT_sigma0, F_max_free, reached_accuracy.

## Provenance / methods (short form)

- Functional: PBE + D3(BJ), IVDW=12, ENCUT=520 eV, PREC=Accurate, LASPH, ADDGRID.
- Bulk: ISIF=3 (full relaxation), ISMEAR=1 for Pd (σ=0.10), ISMEAR=0 for PdO/PdO₂ (σ=0.05).
- Slab: ISIF=2 (ionic only, cell fixed from bulk), bottom-half atoms FixAtoms, vacuum 20 Å, LDIPOL=True, IDIPOL=3.
- MLIP: MACE-MH mh-1 head oc20_usemppbe, float64, cueq enabled, D3(BJ, xc=pbe) dispersion active. LBFGS fmax=0.05 eV/Å, 200–300 steps.
- DFT L2: same slab INCAR as L1 (vacuum); EDIFFG=-0.03 eV/Å; ISPIN=2.
- Convention: **E_bind is defined against isolated gas monomers** (CO gas + CH₃O_radical gas). More negative = stronger binding.

## Sign conventions & unit conventions

| quantity | unit | sign |
|---|---|---|
| E_bulk, E_slab, E_MLIP, E_DFT | eV | absolute (arbitrary offset, only diffs meaningful) |
| E_bind | eV | more negative = stronger binding |
| ΔE (relative to global min) | eV | ≥0 |
| lattice a/b/c | Å | positive |
| force | eV/Å | positive magnitude |
| d_reactive, d_PdO, d_anchor_surf | Å | positive |
| rumpling | Å | positive |

## Regenerate

```bash
python scripts/build_paper_data.py
```
The script is idempotent and re-parses OUTCARs each time so completed DFT jobs are automatically picked up in **07_dft_results.csv**.
