# Selection formulae — E_MLIP, E_bind_MLIP, ΔE, interaction energy

Read-only audit of what the code actually computes. All formulae reference
actual source lines.

## 1. E_MLIP (raw MACE total energy)

`E_MLIP` = MACE-MH+D3(BJ)+cueq PBE-flavor **total electronic energy** in **eV**
of the fully relaxed slab+adsorbate system, evaluated with:

- model = `mh-1`, head = `oc20_usemppbe`
- `default_dtype='float64'`, `enable_cueq=True`
- `dispersion=True, damping='bj', dispersion_xc='pbe'`
- LBFGS relaxation, `fmax=0.05 eV/Å`, `steps=200–300`
- bottom 50 % of slab atoms (`z < median(z)`) frozen via `FixAtoms`

Source: `scripts/run_mace_phase1.py:104-183` (singles) and
`scripts/run_mace_phase2.py:151` (coads). Stored under key `"E"` in each
`MLIP_phase1/unique_{CO,CH3O}.json` and
`MLIP_phase2_filtered/unique_SetA.json`.

**Sign convention**: lower `E_MLIP` = lower total energy. Absolute values
have no physical meaning (they include arbitrary reference offsets from
the PBE pseudopotential). **Cross-surface comparison of raw `E_MLIP` is
meaningless** because each surface has a different atom count and clean-slab
reference.

## 2. Clean-slab and gas reference energies

Computed once with the **same MACE+D3+cueq calculator** in
`scripts/compute_mace_d3_references.py`; stored in
`calculations/G3_adsorption/mace_d3_references.json`:

| key | value (eV) | how computed |
|---|---|---|
| slab.S1 | −389.5571 | MACE relax of `G2_slab/S1_Pd100/CONTCAR`, bottom half fixed |
| slab.S2 | −559.9475 | same for S2 |
| slab.S3 | −673.2548 | same for S3 |
| slab.S3b | −531.0268 | same for S3b |
| slab.S4 | −738.9077 | same for S4 |
| gas.CO | −14.3639 | isolated CO in vacuum cell, MACE relax |
| gas.CH3O_radical | −24.1549 | isolated CH₃O radical (spin-polarized inside MACE) |
| gas.CH3OH | −29.8192 | isolated methanol |
| gas.H2 | −6.9266 | isolated H₂ |
| gas.CH3O_ref | −26.3559 | = E(CH3OH) − ½ E(H2) — the workplan §T1.18 MeOH(U) reference |

All references were computed with **exactly the same MACE calculator, D3(BJ)
dispersion, cueq acceleration, and float64 precision** as the slab+ads
relaxations, so `E_bind_MLIP` cancels systematic MLIP offsets between the
gas molecule and the adsorbed state.

## 3. E_bind_MLIP (binding / adsorption energy)

### Singles (CO*, CH₃O*)

```
E_bind_MLIP(CO*, sid, idx)   = E_MLIP(sid, CO,  idx)   − E_slab(sid) − E_gas(CO)
E_bind_MLIP(CH3O*, sid, idx) = E_MLIP(sid, CH3O, idx) − E_slab(sid) − E_gas(CH3O_radical)
```

Units: **eV**. **Sign convention: more negative = stronger binding**
(convention consistent with catalysis literature and workplan §T1.18).

Source: `scripts/build_violin_v5add.py:59-63` and the auditor
`scripts/audit_selection.py:E_bind()`.

Note: the workplan defines two thermodynamic references for CH₃O* —
radical (`ΔG_CH3O*^rad`) and MeOH(U) (`ΔG_CH3O*^MeOH(U)`). **v4/v5
selection and violin plots use the radical reference** (`CH3O_radical`),
not `CH3O_ref = E(CH3OH) − ½ E(H2)`. Both are stored in the references
file; only `CH3O_radical` is currently consumed downstream.

### Coadsorption

```
E_bind_MLIP(CO*+CH3O*, sid, idx) = E_MLIP(sid, coads, idx)
                                   − E_slab(sid)
                                   − E_gas(CO)
                                   − E_gas(CH3O_radical)
```

This is **definition A** (binding energy against isolated gas monomers),
not an interaction energy. Interaction energy against two separately
adsorbed monomers,

```
E_interaction = E(slab+CO+CH3O) + E(slab) − E(slab+CO) − E(slab+CH3O),
```

**is not computed anywhere in the current code**. This is by design at
Phase 1: we want a single reference for both single- and coadsorption to
allow the descriptor map in T1.19. Interaction energy will be relevant in
Phase 2 (TS analysis) and can be reconstructed from paired DFT results
post-hoc.

## 4. ΔE (relative energy inside a group)

```
ΔE(sid, ads, idx) = E_MLIP(sid, ads, idx) − min_{idx' ∈ valid_pool(sid, ads)} E_MLIP(sid, ads, idx')
```

Applied **inside a single (sid, ads) group** — never across surfaces.
Values in `picks_analysis.csv` and `proposed_additions_only.csv` come
from this formula (`scripts/refine_v5_additions.py:build_addition_row`
line ~236 and `scripts/select_top5_v5.py:main` line ~278).

## 5. Ranking equivalence

Because inside a single (sid, ads) group the references (`E_slab(sid)`,
`E_gas(CO)`, `E_gas(CH3O)`) are constants, subtracting them shifts every
`E_MLIP` by the **same amount**. Therefore:

- **rank(E_MLIP)  ≡  rank(E_bind_MLIP)  ≡  rank(ΔE)** inside one group.
- The selector consumes `E_MLIP` (or equivalently `ΔE`) because it is
  cheaper and unambiguous.
- The violin plot consumes `E_bind_MLIP` because it needs a physically
  meaningful axis that can be compared **across surfaces** on the same
  y-scale.

Both values are exposed side-by-side in
`audit/selection_candidate_rationale.csv` for verification.

## 6. What selection uses vs. what visualisation uses

| stage                  | value used            | source line              |
|------------------------|-----------------------|--------------------------|
| `select_top5_v4.py`    | `E_MLIP` (`r['E']`)   | `picks.sort(key=lambda r: r['E'])` |
| `select_top5_v5.py`    | `E_MLIP` (`r['E']`)   | `pool.sort(key=lambda r: r['E'])` |
| `refine_v5_additions`  | `E_MLIP` + ΔE         | `dE = E − E_min` |
| `build_violin_*.py`    | `E_bind_MLIP`         | `d = [r['E']−E_SLAB[sid]−E_ads_sum for r in recs]` |

**Answer to the direct question**: "Why does the selector use `E_MLIP` and
the plot use `E_bind_MLIP`?" — because ranking inside one (sid, ads)
group is unaffected by any constant shift, so it is safe and cheap for
the selector to use `E_MLIP`; but the plot compares **different surfaces
side-by-side**, so it must subtract each surface's own reference to give
a physically meaningful axis. The two orderings are identical **within**
a group; they only differ **between** groups.
