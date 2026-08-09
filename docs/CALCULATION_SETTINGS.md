# Calculation settings — full audit trail

All settings currently in use across the pipeline, extracted from the actual
INCAR/KPOINTS/POTCAR files in the repo. This is the document to review
before submitting any batch to another server.

Common convention:
- **Functional**: PBE + D3(BJ) (`IVDW = 12`)
- **POTCAR**: PAW_PBE, **Pd_pv** (semicore 4p, 16 valence electrons), C, H, O
  (all standard). Library at `/home/hyunjin/POTENTIAL/potpaw_PBE/`.
- **Cutoff**: `ENCUT = 520 eV`, `PREC = Accurate`, `LASPH`, `ADDGRID`, `LORBIT = 11`
- **Spin**: `ISPIN = 2` everywhere (workplan mandate; per session Q&A —
  MAGMOM only set for CH₃O radical). NUPDOWN=1 on CH₃O.
- **D3 flavour**: Becke-Johnson damping (`IVDW = 12` is D3-BJ in VASP 6.x).

---

## 1. G1 Bulk optimisation (T1.1–T1.4)

Fully-relaxed cell + ions, ISIF=3.

| tag | Pd (metal) | PdO (oxide) | PdO₂ (oxide) | PdO₂ (no-D3 control) |
|---|---|---|---|---|
| ENCUT | 520 | 520 | 520 | 520 |
| PREC | Accurate | Accurate | Accurate | Accurate |
| LASPH | .TRUE. | .TRUE. | .TRUE. | .TRUE. |
| ADDGRID | .TRUE. | .TRUE. | .TRUE. | .TRUE. |
| ISPIN | 2 | 2 | 2 | 2 |
| IVDW | **12** (D3-BJ) | **12** | **12** | **absent** (control) |
| EDIFF | 1e-06 | 1e-06 | 1e-06 | 1e-06 |
| NELM | 200 | 200 | 200 | 200 |
| NELMIN | 5 | 5 | 5 | 5 |
| ALGO | Normal | Normal | Normal | Normal |
| NCORE | 1 | 1 | 1 | 1 |
| LREAL | .FALSE. | .FALSE. | .FALSE. | .FALSE. |
| LWAVE | .FALSE. | .FALSE. | .FALSE. | .FALSE. |
| LCHARG | .FALSE. | .FALSE. | .FALSE. | .FALSE. |
| ISYM | **2** (symmetry ON) | 2 | 2 | 2 |
| IBRION | 2 | 2 | 2 | 2 |
| NSW | 200 | 200 | 200 | 200 |
| **ISIF** | **3** (full cell relax) | 3 | 3 | 3 |
| ISMEAR | **1** (Methfessel-Paxton) | **0** (Gaussian) | 0 | 0 |
| SIGMA | 0.10 | 0.05 | 0.05 | 0.05 |
| EDIFFG | **-0.01** (tight, bulk) | -0.01 | -0.01 | -0.01 |
| k-mesh | 12×12×12 Γ | 8×8×6 Γ | 6×6×8 Γ | 6×6×8 Γ |

Notes:
- Metal uses Methfessel-Paxton smearing (ISMEAR=1) — appropriate for
  metallic Fermi surface.
- Oxides use Gaussian smearing (ISMEAR=0) with sharper SIGMA — appropriate for
  band-gap systems.
- Bulk gets the tightest force criterion (0.01 eV/Å); slab & adsorbate loosen to
  0.03 eV/Å.

---

## 2. G2 Clean slab (T1.5–T1.9) — ✅ CONVERGED

Ionic relaxation only (`ISIF = 2`), cell fixed from bulk equilibrium lattice.

| tag | S1 Pd(100) | S2 PdO(101)/Pd(100) | S3 PdO(100) | S3b PdO(100) Pd-term | S4 PdO₂(110) |
|---|---|---|---|---|---|
| ENCUT | 520 | 520 | 520 | 520 | 520 |
| ISPIN | 2 | 2 | 2 | 2 | 2 |
| IVDW | 12 | 12 | 12 | 12 | 12 |
| EDIFF | 1e-06 | 1e-06 | 1e-06 | 1e-06 | 1e-06 |
| NELM | 200 | 500 | 500 | 500 | 500 |
| NELMIN | 5 | 5 | 5 | 5 | 5 |
| LREAL | Auto | Auto | Auto | Auto | Auto |
| LWAVE / LCHARG | F / F | F / F | F / F | F / F | F / F |
| **ISYM** | **0** (symmetry OFF) | 0 | 0 | 0 | 0 |
| IBRION | 2 | 2 | 2 | 2 | 2 |
| NSW | 300 | 300 | 300 | 300 | 300 |
| ISIF | 2 | 2 | 2 | 2 | 2 |
| ISMEAR | 1 | 0 | 0 | 0 | 0 |
| SIGMA | 0.10 | 0.05 | 0.05 | 0.05 | 0.05 |
| EDIFFG | -0.03 | -0.03 | -0.03 | -0.03 | -0.03 |
| LDIPOL | .TRUE. | .TRUE. | .TRUE. | .TRUE. | .TRUE. |
| IDIPOL | 3 | 3 | 3 | 3 | 3 |
| KSPACING | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 |
| natoms | 80 | 112 | 128 | 104 | 144 |
| vacuum | 20 Å | 20 Å | 20 Å | 20 Å | 20 Å |
| fix mask (G2 rule) | bottom 2 layers (32) | (40) | (32) | (32) | (42) |
| E_slab (eV, E₀) | -434.380 | -618.565 | -724.103 | -570.772 | -788.493 |
| F_max free (eV/Å) | 0.019 | 0.026 | 0.023 | 0.028 | 0.025 |

Common:
- Asymmetric slab: bottom-half fixed via selective dynamics.
- Vacuum 20 Å along surface normal.
- Dipole correction along z (IDIPOL=3) since asymmetric slab has net dipole.
- KSPACING = 0.25 Å⁻¹ (project standard, matches ~3×3×1 mesh on typical
  supercells).
- ISYM=0: symmetry off for adsorbate-carrying / defective slabs.

---

## 3. MLIP screening (T1.10–T1.14, via MACE)

**Not DFT**, but included for completeness. `scripts/run_mace_phase*.py`.

| tag | value |
|---|---|
| model | mace_mp `mh-1` |
| head | `oc20_usemppbe` |
| default_dtype | float64 |
| enable_cueq | True |
| device | cuda |
| dispersion | True |
| damping | 'bj' |
| dispersion_xc | 'pbe' |
| optimiser | LBFGS |
| fmax (singles) | 0.05 eV/Å |
| max steps (singles) | 200 |
| fmax (coads) | 0.05 eV/Å |
| max steps (coads) | **300** (bumped from 200 per S4/S3 convergence concern) |
| fix mask | same as T1.16 (bottom-half by median z) |
| dedup key (singles) | `(round(E,2), intra-ads geometric fingerprint)` |
| dedup key (coads filtered) | `(round(E,2), round(d_reactive,1), loose fingerprint)` |

Reference energies (in `mace_d3_references.json`) computed with **exactly the
same MACE calculator** — all E_bind values in downstream CSVs cancel MLIP
systematic offset:

| kind | value (eV) |
|---|---|
| slab.S1 | -389.5571 |
| slab.S2 | -559.9475 |
| slab.S3 | -673.2548 |
| slab.S3b | -531.0268 |
| slab.S4 | -738.9077 |
| gas.CO | -14.3639 |
| gas.CH3O_radical | -24.1549 |
| gas.CH3OH | -29.8192 |
| gas.H2 | -6.9266 |
| gas.CH3O_ref (=E_CH3OH − ½ E_H2) | -26.3559 |

---

## 4. T1.16 DFT Level-1 (vacuum) — 86 candidates

Same core INCAR as G2 slab, adapted for adsorbate + slab.
`calculations/T1_16_DFT_L2/{sid}/{ads}/{ads}_idx{XXXXX}/`.

| tag | metal (S1) | oxide (S2, S3, S3b, S4) |
|---|---|---|
| ENCUT | 520 | 520 |
| PREC | Accurate | Accurate |
| LASPH / ADDGRID | .TRUE. | .TRUE. |
| ISPIN | 2 | 2 |
| IVDW | 12 | 12 |
| EDIFF | 1e-06 | 1e-06 |
| NELM | 500 | 500 |
| NELMIN | 5 | 5 |
| ALGO | Normal | Normal |
| NCORE | 1 | 1 |
| LREAL | Auto | Auto |
| LWAVE / LCHARG | F / F | F / F |
| LORBIT | 11 | 11 |
| ISYM | 0 | 0 |
| IBRION | 2 | 2 |
| NSW | 300 | 300 |
| ISIF | 2 | 2 |
| ISMEAR | 1 | 0 |
| SIGMA | 0.10 | 0.05 |
| EDIFFG | -0.03 | -0.03 |
| LDIPOL | .TRUE. | .TRUE. |
| IDIPOL | 3 | 3 |
| KSPACING | 0.25 | 0.25 |

**MAGMOM**: not set → VASP default 1.0 μ_B/atom initialisation. Empirically
converges to small values for these Pd/PdO/PdO₂ systems (S1 CO* 5 candidates
already verified).

**fix mask**: bottom-half of substrate atoms by median z (rule applied at
setup time in `scripts/select_top5_v4.py::fix_bottom_half`). Because this
mask includes the (high-z) adsorbate atoms in the median computation, the
resulting fixed count differs from G2 by up to ~10 atoms (documented in
methodology audit).

**POTCAR order**: after `sort=True`, POSCAR groups species alphabetically —
`[C, (H), O, Pd]`. POTCAR concatenated in this exact order.

**Status**: 7/86 DONE (S1/CO ×5, S3/CH3O ×2), 79 PENDING.

---

## 5. T1.17 DFT Level-2 (VASPsol solvation)

`calculations/T1_17_VASPsol/{sid}/{ads}/{ads}_idxXXXXX/` (adsorbate slabs)
and `T1_17_VASPsol/{sid}_clean/` (clean-slab references).

Same INCAR as T1.16 **plus** the VASPsol block and explicit restart flags:

| tag | value | comment |
|---|---|---|
| ENCUT | 520 | same as L1 |
| PREC | Accurate | same |
| ISPIN | 2 | same |
| IVDW | 12 | same |
| EDIFF | 1e-06 | same |
| NELM | 500 | same |
| ALGO | Normal | same |
| NCORE | 1 | same |
| LREAL | Auto | same |
| LWAVE / LCHARG | F / F | same |
| ISYM | 0 | same |
| IBRION | 2 | same — clean slab also relaxes (per reviewer, not static) |
| NSW | **200** | slightly reduced from 300 (starting from L1 CONTCAR) |
| ISIF | 2 | same |
| ISMEAR | 1 (S1) / 0 (oxides) | same |
| SIGMA | 0.10 / 0.05 | same |
| EDIFFG | -0.03 | same |
| LDIPOL | .TRUE. | same |
| IDIPOL | 3 | same |
| KSPACING | 0.25 | same |
| **ISTART** | **0** | **fresh electronic start** — L1 wrote LWAVE=.FALSE. → no WAVECAR to reuse. Restart is geometric only (CONTCAR → POSCAR). |
| **ICHARG** | **2** | initial CHG from superposition of atomic charge densities. |
| **LSOL** | **.TRUE.** | activates VASPsol |
| **EB_K** | **32.6** | methanol static dielectric ε |
| **TAU** | **0** | exclude cavitation/non-electrostatic term (workplan) |
| ~~LAMBDA_D_K~~ | **absent** | deliberately not set — would activate the linearised Poisson–Boltzmann electrolyte model (Debye screening), which is a separate physical effect. |

**Clean slab constraints**:
- Built by stripping adsorbate (C, all H, and any O within 1.5 Å of C, PBC-aware)
  from that surface's top-1 T1.16 CONTCAR.
- Selective-dynamics mask + POTCAR order + cell inherited exactly from T1.16
  ads calc → auto-consistency check verifies cell/Pd count/mask count/POTCAR
  TITEL/KSPACING match ads dir (currently 100 % pass for S1, S3).

**One clean slab per surface**, using the surface's top-1 ads (prefers CO >
CH₃O > coads for least substrate distortion). Small mask mismatch to other
ads types on the same surface is a documented ~10-50 meV limitation.

**Status**: 4 dirs generated (S1: CO+clean; S3: CH₃O+clean). 15 more will be
created as more T1.16 completes.

---

## 6. Gas-phase references (for T1.18 μ_gas)

`calculations/gas_references/`. Isolated molecule in a cubic vacuum box.

Common:
- ENCUT 520, PREC Accurate, LASPH, ADDGRID, IVDW=12 — **identical to slab
  settings** so systematic PBE-D3 error cancels in ΔG_ads formula.
- NCORE 1, LREAL Auto, ISYM 0.
- **ISMEAR = 0, SIGMA = 0.01** (very sharp — isolated molecule).
- **EDIFFG = -0.01** (tighter than slab, since molecule geometry sensitive).
- KPOINTS **explicit Γ-only file** (1×1×1 Γ, no KSPACING active).
- box side 15 Å cubic (primary); 20 Å variants for polar molecules (CH₃OH, CH₃O)
  as paired-static box-size sanity checks.

### Per-molecule config

| molecule | box | variants | ISPIN | NUPDOWN | MAGMOM | LSOL | IBRION / NSW |
|---|---|---|---|---|---|---|---|
| CO | 15 Å | vacuum + vaspsol | 1 | – | – | off / on | 2 / 200 |
| CH₃OH | 15 Å | vacuum + vaspsol | 1 | – | – | off / on | 2 / 200 |
| H₂ | 15 Å | vacuum + vaspsol | 1 | – | – | off / on | 2 / 200 |
| **CH₃O (radical)** | 15 Å | vacuum + vaspsol | **2** | **1** | `0 0 0 0 1` (C,H,H,H,O) | off / on | 2 / 200 |
| CH₃O (paired static) | 15 Å + 20 Å | `_15A_static` + `_20A_static` | 2 | 1 | `0 0 0 0 1` | off | **-1 / 0** |
| CH₃OH (paired static) | 15 Å + 20 Å | `_15A_static` + `_20A_static` | 1 | – | – | off | -1 / 0 |

**VASPsol block** (identical to T1.17): `LSOL=.TRUE., EB_K=32.6, TAU=0`. No
`LAMBDA_D_K`.

**Restart**: gas molecules relax freely, no ISTART override.

**Paired-static workflow**:
1. Run primary relax dirs (`CH3O_vacuum`, `CH3OH_vacuum`) first.
2. When each finishes: `python scripts/prepare_static_from_relax.py` — reads
   the 15 Å CONTCAR, extracts coordinates ONLY (not the cell), writes into
   both `_15A_static` and `_20A_static` dirs' POSCARs with the correct box.
3. Submit the 4 static dirs. Target: |E(20 Å) − E(15 Å)| < 5 meV finite-cell
   error.

---

## 7. What is NOT set anywhere (deliberate)

- `LAMBDA_D_K` — electrolyte Debye length. Would activate linearised
  Poisson–Boltzmann. Not needed for neutral methanol solvation.
- `LSOL` in vacuum dirs — solvation off.
- `NUPDOWN` — only in CH₃O radical dirs. Slab calcs let VASP find the
  moment.
- Dipole correction (`LDIPOL`, `IDIPOL`) — **on** for asymmetric slabs
  (G2, T1.16, T1.17). **off** for gas (no dipole issue in cubic box).
- `KSPACING` — on for slabs (0.25). Off for gas (explicit Γ KPOINTS instead).
- POTCAR "Pd" (10-electron) — not used; project uses **Pd_pv** (16-electron
  semicore) throughout, per project decision.

---

## 8. Convention summary for ΔG (T1.18 post-processing)

```
ΔG_CO*(vac)         = G(slab+CO)_vac      − G(slab)_vac      − μ_CO(gas, vacuum)
ΔG_CO*(sol)         = G(slab+CO)_sol      − G(slab)_sol      − μ_CO(gas, vaspsol)

ΔG_CH3O*(rad, vac)  = G(slab+CH3O)_vac    − G(slab)_vac      − G(CH3O_radical, vacuum)
ΔG_CH3O*(rad, sol)  = G(slab+CH3O)_sol    − G(slab)_sol      − G(CH3O_radical, vaspsol)

ΔG_CH3O*(MeOH(U))   = G(slab+CH3O) + ½ G(H2) − G(slab) − G(CH3OH) − eU
                      (CHE, all references in matching phase)

ΔG_coads            = G(slab+CO+CH3O) − G(slab) − μ_CO − G(CH3O)
                      (definition A: against isolated gas monomers)
```

**Do not mix vacuum-slab and solvated-adsorbate energies in one formula.**
Choose a phase (vacuum or vaspsol) per formula and use references from that
phase throughout.

Free-energy corrections to be added at analysis time (ASE Thermochemistry):
- Adsorbate: **harmonic approximation** (finite-difference frequency calc, IBRION=5)
- Gas molecule: **ideal-gas** (rotational + translational + vibrational, with
  molecular symmetry number)

---

## 9. Open decisions before T1.18 (analysis-time)

- CO reservoir: gas feed (`μ_CO(g)`) or dissolved CO?
- CH₃OH reservoir: gas or liquid methanol at activity 1? If liquid → add
  Δμ_solvation correction to G(CH3OH_vaspsol).
- H₂ reservoir: standard CHE `½ G_H2(g)`.
- Applied potential eU: value + sign convention.

These do NOT affect the DFT setup, only the post-processing formula.

---

## 10. Regeneration

Every DFT dir under `calculations/` is idempotently regenerable:

| stage | script |
|---|---|
| G1 bulk | manually set up (historical) |
| G2 slab | manually set up (historical) |
| MLIP pool | `scripts/run_mace_phase{1,2,3}.py`, `refilter_phase{2,3}_geometry.py` |
| T1.15 shortlist | `scripts/select_top5_v{4,5}.py`, `refine_v5_additions.py` |
| T1.16 v4 (70) | `scripts/setup_dft_v3.py` |
| T1.16 v5-add (16) | `scripts/setup_v5add16.py` |
| T1.17 VASPsol | `scripts/setup_t1_17_vaspsol.py` |
| gas references | `scripts/setup_gas_references.py` |
| paired-static POSCARs | `scripts/prepare_static_from_relax.py` (after primary relax done) |
| paper CSVs | `scripts/build_paper_data.py` |
| methodology audit | `scripts/audit_selection.py` |
