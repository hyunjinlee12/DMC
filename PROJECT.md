# Project: Pd / PdO / PdO₂ surface DMC formation — DFT study

<!-- The Director reads this as its brief, every session. This is a DFT MECHANISTIC
     study, NOT an ML property-prediction project: there are no train/val/test
     splits and no held-out MAE. The deliverables are a descriptor map (Phase 1)
     and a DMC Gibbs free-energy profile (Phase 2). -->

## Goal
Determine whether anodic oxidation of Pd nanocube {100} facets (Pd⁰ → Pd²⁺ → Pd⁴⁺)
**promotes or suppresses** electrochemical DMC formation. Do this by computing how the
CO* / CH₃O* adsorption balance and the DMC reaction barriers change across four
surfaces of increasing oxidation state.
Reaction: 2 CH₃OH + CO − 2e⁻ → (CH₃O)₂CO + 2H⁺.

## Target property & system
Surfaces along the oxidation axis S1 → S4:
- S1  Pd(100)                 — Pd⁰ metallic baseline
- S2  1 ML PdO(101)/Pd(100)   — Pd²⁺ surface oxide / Pd–PdO interface
- S3  bulk / O-rich PdO(100)  — Pd²⁺ thermodynamic / O-rich limit
- S4  PdO₂(110)               — Pd⁴⁺ high-valence reference

Quantities computed:
- ΔG_CO* and ΔG_CH₃O* (radical reference + CHE / MeOH(U) reference)
- Descriptor map: ΔG_CO* vs ΔG_CH₃O*^MeOH(U) across S1–S4  (Phase 1 deliverable)
- DMC Gibbs free-energy profile + TS1 / TS2 barriers          (Phase 2 deliverable)
- Side-path barriers: CO + O_lattice → CO₂ + V_O, methanol oxidation

## Approach
Per structure class, the 3-stage pipeline:
  AutoAdsorbate (generate candidates) → MACE foundation MLIP (relax + RANK, screening only)
  → VASP PBE+D3 (+VASPsol) (final judgment). TS via VTST CI-NEB / ASE-NEB (+ dimer).
Two phases:
- Phase 1 (T1.1–T1.20): bulk opt → clean slab build+validate → adsorption sampling →
  MLIP ranking → DFT adsorption energies → descriptor map → select Phase-2 surfaces.
- Phase 2 (T2.1–T2.8): reactive-pair + intermediate sampling → endpoint opt →
  TS1/TS2 → side-path → DMC Gibbs profile + cross-surface comparison.

## Data / resources
- Initial structures: Materials Project via mp-api (MP_API_KEY in env) / pymatgen —
  fcc Pd, tetragonal PdO, rutile/hydrophilite-like PdO₂.
- Cluster VASP (+VASPsol) via SLURM for all DFT. GPU node for MACE (mace-torch).
- Reference docs (read as authoritative project context, in docs/):
  DMC_Pd_workplan.md · DMC_Pd_package_guideline.md ·
  DMC_Pd_references_summary.md · DMC_Pd_student_guide.md

## Constraints
- All slabs asymmetric: adsorption on top side only; bottom 2 layers (or bottom
  30–40%) fixed; vacuum 20 Å; LDIPOL=.TRUE., IDIPOL=3; ISYM=0.
- DFT: PBE+D3 (IVDW=12), ENCUT=520, PREC=Accurate, LASPH=.TRUE., ADDGRID=.TRUE.,
  ISPIN=2. Smearing: Pd metal ISMEAR=1/SIGMA=0.10; oxides ISMEAR=0/SIGMA=0.05.
  VASPsol (LSOL=.TRUE., EB_K=32.6, TAU=0): apply ONLY after vacuum convergence.
- MLIP is RANKING ONLY — never use MACE absolute energies for a conclusion; always
  cross-validate distance-bin representatives with DFT.
- AutoAdsorbate CO*/CH₃O* must use mode='all' (heuristic = top site only).
  Co-adsorption uses the custom site-pair wrapper, not bare get_populated_sites().
- Decision gates G1–G4 are HARD STOPS. Do not start sampling before G2, or Phase 2
  before G3. Show the researcher the matching checkpoint (A–E) before crossing a gate.
- Reproducibility: record INCAR/KPOINTS/POTCAR + package versions (esp. the MACE
  model checkpoint) with every result.

## Definition of done
- G1: 3 converged + validated bulk structures (lattice vs experiment; ENCUT/k-mesh).
- G2: 4 validated clean slabs (rumpling, Pd–O bond length, coordination, dipole,
  termination stability).
- G3: adsorption-energy table + descriptor map + selected Phase-2 surfaces (Case A–D).
- G4: per-surface DMC Gibbs profile + side-path comparison + a conclusion on whether
  anodic oxidation promotes or suppresses DMC, benchmarked against Shi 2024 Angew
  (pure Pd TS1/TS2 = 1.08/0.85 eV, Pd₃Cu = 0.86/0.79 eV).
