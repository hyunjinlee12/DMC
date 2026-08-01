# Selection methodology audit — v4 (70) + v5-add (16)

Read-only audit. Every claim below cites the actual selection code (path
+ line range) so anyone can reproduce.

> **Ground rule this audit follows**: code first, narrative second. Where
> the report and the code disagree, the code is what actually ran; the
> divergence is called out in `methodology_issues.md`.

---

## 1. Full selection workflow

Step-by-step pipeline from raw AutoAdsorbate output to DFT inputs. Numbers
in parentheses are the counts I re-derived from
`audit/selection_group_summary.csv`.

| # | stage | script | input | pass / drop | output |
|---|---|---|---|---|---|
| 1 | Initial adsorbate placement (AutoAdsorbate `mode='all'`; coads = custom pair wrapper) | `data-curator` agent + `coads_guide/*.py` | slab CONTCAR + AutoAdsorbate `Fragment(SMILES)` | 87–498 raw candidates per (sid, ads); 3 764–11 598 raw per coads group | `{surface}/{ads}/candidates.traj`, `coads_guide/SetA.traj` |
| 2 | MLIP relaxation (singles) | `scripts/run_mace_phase1.py` `relax_candidate()` (L104-183) | `candidates.traj`, MACE-MH+D3+cueq, LBFGS `fmax=0.05 eV/Å`, `steps=200`, bottom-half fixed | not converged → dropped (typical drop 0–2 %); n_steps recorded | `MLIP_phase1/relaxed_{CO,CH3O}.traj` + `summary.json` |
| 3 | MLIP relaxation (coads) | `scripts/run_mace_phase2.py` `L151-190` | `coads_guide/SetA.traj`, same MACE settings, `fmax=0.05`, `steps=300` (bumped from 200 per S4/S3 committee retry) | fingerprint + E-round dedup, `groups[(round(E,2), fingerprint)]` | `MLIP_phase2/relaxed_SetA.traj` + `ranked_SetA.json` |
| 4a | Convergence + intramolecular validity (singles) | `scripts/select_top5_v4.py` `valid_CO/CH3O()` (L32-56) | `unique_{ads}.json` | rejects C–O outside [1.05, 1.30] Å (CO*) / [1.30, 1.55] Å (CH₃O*); rejects if 2nd O within 2.0 Å (CO₂-like); rejects if any C–H outside [0.90, 1.25] Å | in-memory valid pool |
| 4b | Convergence + geometry filter (coads) | `scripts/refilter_phase2_geometry.py` `geometry_valid()` (L51-71) | `MLIP_phase2/relaxed_SetA.traj` | drops fragmented (max intra-pair >2.5 Å from ADS_DIRECT_MAX), collapsed (min <1.0 Å from ADS_DIRECT_MIN), CO_bond_broken, methoxy_OC_broken, CH_broken, unconverged | `MLIP_phase2_filtered/unique_SetA.json` |
| 5 | Site + fingerprint classification (Phase-1) | `scripts/run_mace_phase1.py` `classify_site()` (L54-85) — uses `overall_min+0.5 Å` buffer (NOT fixed 2.6 Å) | relaxed atoms | atop/bridge/hollow/unknown label stored per candidate | `site_type` field in each `unique_{CO,CH3O}.json` entry |
| 6 | Structural dedup | `run_mace_phase1.deduplicate()` (L185-199) singles; `run_mace_phase2.py` L184-190 coads; `refilter_phase2_geometry.py` `loose_fingerprint()` (L90-99) re-dedup coads | relaxed pool | grouped by `(round(E,2), fingerprint)` → one representative per group | `unique_{CO,CH3O}.json`, `unique_SetA.json` (deduplicated) |
| 7a | **v4 selector** (E + xy dedup only, no site diversity) | `scripts/select_top5_v4.py` `pick_top5()` (L111-148) | `unique_*.json` + `relaxed_*.traj` | for each candidate in ascending E: keep if no prior pick has \|ΔE\|<0.03 eV **AND** xy MIC dist <1.5 Å from anchor(s); stop at 5 | `DFT_shortlist_v3/{sid}/{ads}/{rank}_{ads}_idx{idx}.vasp` + `summary.json` |
| 7b | **v5 selector** (diversity-aware) | `scripts/select_top5_v5.py` `select_single()` (L189-227) + `select_coads()` (L229-269) | same pool | 1) global E min → slot 0. 2) then walk sorted pool; accept if new `site_class` or new `region` and ΔE≤0.5 eV (single) / ≤0.8 eV (coads) and target `d_bin` not yet covered (coads). 3) fill remaining with lowest-E non-dup. | `DFT_shortlist_v5/summary.csv` + selection log in `report.md` |
| 8 | Priority triage (v5) | `scripts/refine_v5_additions.py` `refine()` (L214-249) | v5 picks minus v4 picks (=27 v5-new) | MUST for reactive-close/loose gap-fillers on S3/S3b; RECOMMENDED for S1/S2 reactive+thermo reps; OPTIONAL for single-ads news; REVIEW-NEEDED for coads news outside researcher list | `proposed_additions_only.csv` |
| 9 | DFT input build (v4=70; v5add16=16 approved subset) | `scripts/setup_dft_v3.py` (70) + `scripts/setup_v5add16.py` (16) | approved shortlist + reference POTCARs from `/home/hyunjin/POTENTIAL/potpaw_PBE/` | POSCAR (`write(..., sort=True, direct=True, vasp5=True)`), INCAR (S1 metal ISMEAR=1 / others oxide ISMEAR=0), POTCAR concatenated in POSCAR species order, `submit_vasp_gpu.sh` copied from a reference dir | `T1_16_DFT_L2/*` (70) + `T1_16_DFT_L2_v5add16/*` (16) |
| 10 | DFT PBE-D3 vacuum relaxation | VASP 6.4.3 + IVDW=12; EDIFFG=−0.03 eV/Å, KSPACING=0.25, ENCUT=520, LDIPOL, ISPIN=2 | POSCAR/INCAR/POTCAR | `reached required accuracy` on OUTCAR = pass | OUTCAR/CONTCAR/OSZICAR |

**Bottom line**: the current shortlist funnel is `raw → MLIP relax → validity filter → structural dedup → diversity-aware pick → DFT`. Every step is scripted and each output file is deterministic given the same input+random seed. `run_mace_phase1.py` uses a **different** site-classification rule (`overall_min+0.5 Å` buffer) than v5's re-classifier (fixed 2.6 Å cutoff); the field `unique_*.json:site_type` and the field `picks_analysis.csv:site` are **not the same quantity** — see §4 and `methodology_issues.md`.

---

## 2. E_MLIP vs E_bind_MLIP

See the companion file **`selection_formulae.md`** for the exact formulae,
sign conventions, references, and the answer to *"why selector uses
E_MLIP and plot uses E_bind_MLIP"*. Short version:

- `E_MLIP` = MACE-MH+D3+cueq total electronic energy of the relaxed
  slab+ads system (eV). No physical baseline; only comparable **within a
  single (sid, ads) group**.
- `E_bind_MLIP` = `E_MLIP − E_slab(sid) − Σ E_gas(ads_i)`, with all
  references computed with the **same MACE calculator**. Lower = stronger
  binding.
- Selector (`select_top5_v[45].py`) uses `E_MLIP` because inside a group
  `rank(E_MLIP) ≡ rank(E_bind_MLIP)` and the subtraction is redundant.
- Violin plots (`build_violin_*.py`) use `E_bind_MLIP` because they
  compare surfaces side-by-side.
- Coadsorption uses **binding definition A** (against isolated gas
  monomers CO + CH₃O_radical), not interaction energy. Interaction energy
  is not computed anywhere yet; it can be reconstructed post-DFT if
  needed.

---

## 3. Why energy is a selection criterion — and why not alone

1. **Ranking**: MLIP+D3 was calibrated per project guide as a *ranking*
   surrogate for PBE-D3 DFT (`docs/DMC_Pd_workplan.md` §T1.14).
   Global minimum candidates are the most likely to survive DFT
   re-relaxation without qualitative reorganization, so slot 0 = global
   min is **always** included (`select_top5_v5.py:189`).

2. **Foundation-MLIP uncertainty**: absolute MACE energies carry ~50 meV
   systematic error vs DFT for adsorption problems (project MLIP
   validation, `reports/G3/…`). Rank inversions of a few meV are
   therefore *not physical* — candidates within ~30 meV of each other
   should be treated as a **degenerate cluster**, not a strict ordering.
   `select_top5_v4.py:DELTA_E_DUP = 0.03 eV` was chosen with this in
   mind, and `select_top5_v5.py` re-uses the same cluster tolerance.

3. **Why not just take the top-5 by E**: a set of 5 near-degenerate
   candidates at the same site fails the workplan requirement (§T1.15
   "거리 bin별 대표 구조 포함") — you can't build a descriptor map from
   five copies of the same physical minimum. Site/xy diversity is
   necessary to *sample the local PES*, not just pin down one minimum.

4. **When we accept high-ΔE picks**: only when they represent a
   qualitatively different chemistry that could plausibly matter for
   the reaction — e.g. `S3b/coads idx=2051` (ΔE=0.55 eV, but only
   reactive-close candidate in the pool). Even 0.55 eV in MLIP is well
   within the range where DFT can invert ordering, and TS1/TS2 barriers
   for DMC formation are 0.7–1.0 eV — so a "high" reactive endpoint may
   turn out to be a viable TS approach vector after DFT.

5. **Different ΔE cap for singles vs coads**: single-ads uses 0.5 eV cap,
   coads uses 0.8 eV — because coads pools are 100× larger and the
   binding-energy range is wider (multiple site + separation
   combinations). A tighter cap would exclude interface pairs on S2.

### The `S3b/coads idx=2051` case
- ΔE = 0.554 eV — high, but the **only reactive-close candidate** in the
  entire filtered S3b coads pool of 1 579 valid structures.
- Kept because a Phase 2 TS1 NEB requires *some* reactive-close endpoint.
  If we exclude it, we effectively concede there is no DMC-relevant
  coadsorption on S3b — a conclusion that must be earned from DFT, not
  from an MLIP ranking.
- Flagged **`MUST-diagnostic`** in `metadata.json`, meaning: run the DFT
  to answer "does this configuration remain bound at PBE-D3 level, or
  does it relax back to a lower-E separated arrangement?". Both outcomes
  are informative.

---

## 4. Site classification

**Anchors**:
- `CO*` anchor = C (the carbon that bonds to the surface; project convention).
- `CH₃O*` anchor = O (the oxygen that binds to the surface; the methyl C
  points into vacuum).
- coads anchors = `C_CO` and `O_CH3O` (the two surface-facing atoms).

**Cutoff**: fixed **2.6 Å** in all v5 code paths (`select_top5_v5.py`,
`refine_v5_additions.py`, `audit_selection.py`). This is a `Pd–C_bond +
~0.5 Å` cushion — captures chemisorbed C above Pd (typical d≈1.9–2.1 Å)
and O above Pd (typical d≈2.0–2.2 Å) while rejecting van-der-Waals
neighbors.

**Class definitions** (v5 `site_label()`):

| class | n_Pd within 2.6 Å | n_O within 2.6 Å | extra condition | example (v5-new idx) |
|---|---|---|---|---|
| `physi` | 0 | 0 | — | S3/CH3O 395, S4/CH3O 383 |
| `atop_Pd` | 1 | 0 | — | S1/CH3O 283, S3b/CH3O picks |
| `atop_O` | 0 | 1 | — | S2/CH3O 496 |
| `br_PdPd` | 2 | 0 | — | S3b/CO 6 |
| `br_PdO` | 1 | 1 | — | S2/CH3O 217, S4/CH3O 329 |
| `br_OO` | 0 | 2 | — | S3/CH3O 385 (OPTIONAL, dropped by researcher) |
| `h3Pd` | 3 | 0 | — | — |
| `h3O` | 0 | 3 | — | — |
| `h3(nPdmO)` | n | m (n+m=3) | — | — |
| `4f_4Pd0O` | 4 | 0 | — | **S1/CO 43** (new 4-fold hollow) |
| any 4+f | n | m (n+m=4+) | — | — |
| `unknown` | — | — | not produced by v5 — MLIP `run_mace_phase1.classify_site` produces this when its `overall_min+0.5` buffer captures 0 substrate atoms | present in raw pool JSONs |

### Cutoff sensitivity (2.4 / 2.6 / 2.8 Å)

The 16 v5-new candidates are re-labelled at each cutoff in
`audit/selection_sensitivity.csv`. Robustness classes:

- **robust**: same label at all three cutoffs (bin also stable if coads)
- **threshold-sensitive-label**: label flips between cutoffs (e.g.
  br_PdPd at 2.6 Å becomes h3 at 2.8 Å or atop_Pd at 2.4 Å)
- **threshold-sensitive-inclusion**: label flip crosses the physi ↔
  chemisorbed boundary
- **diagnostic**: label unstable but candidate was selected precisely
  because we want DFT to arbitrate

Full table: `audit/selection_sensitivity.csv`. Quick summary shipped in
`audit/selection_group_summary.csv` (n_sites_distinct column).

### Known limitations (from `methodology_issues.md`)

1. Fixed 2.6 Å is set for Pd–ligand geometry. For rare Pd–H (<1.9 Å) it
   still works, but a 2.6 Å ball also captures 2nd-shell O_lattice on
   oxide surfaces where surface O–Pd is 2.0 Å + 2nd-shell O at 2.3–2.7 Å.
   Some `atop_Pd` labels on S3/S3b oxide slabs would flip to `br_PdO` at
   2.8 Å; conversely a couple of `br_PdO` on S2 flip to `atop_Pd` at 2.4 Å.
   The sensitivity CSV surfaces which specific candidates flip.

2. Tilted CH₃O with O below and C above (typical bridging methoxide)
   has O anchor 2.0 Å from surface Pd and C at 3.4 Å; anchor=O rule
   handles this correctly. But if MLIP relaxes to *inverted* geometry
   (C down, O up — nonphysical but seen in a few dissociated attempts),
   the anchor rule mislabels it. These are rejected by `valid_CH3O`
   (C–O bond outside 1.30–1.55 Å) before site classification runs.

3. Symmetry-equivalent sites on reconstructed surfaces may split into
   distinct labels because the neighbor pattern differs at the atomic
   level even if the site is topologically equivalent. Conversely,
   different physical sites (e.g. atop Pd_cus vs atop Pd_terrace on S2)
   can share the label `atop_Pd`. **This is why the region tag is
   introduced in v5 (§5).**

4. `physi` is literally "0 substrate atoms within 2.6 Å of anchor." It
   does not mean "no chemistry"; it means the MLIP-relaxed geometry has
   no atom in bonding range. DFT often confirms physisorption for these
   (see S4 CO discussion in §7).

5. **Site labels are diversity descriptors, not final adsorption-site
   claims.** DFT must confirm.

---

## 5. Region classification

`surface_region()` in `select_top5_v5.py:110-114`:

```python
if site=='physi':                    return 'physi'
if n_pd>0 and n_o==0:                return 'metal'
if n_pd==0 and n_o>0:                return 'oxide'
return 'interface'                   # mixed Pd + O neighbors
```

Region is derived **entirely from the neighbor composition**, not from
absolute atom coordinates in the slab. Consequences:

1. On S2 (PdO(101)/Pd(100)), a CO adsorbed on top of an interface Pd
   atom will typically have n_Pd=1, n_O=1 nearby → **`interface`**. A CO
   on the terrace Pd(100) portion will have n_Pd=1, n_O=0 → **`metal`**.
   No absolute geometric partition of the slab is used.

2. `br_PdO` is **generally interface**, but on S3/S3b (bulk PdO(100),
   O-rich) the whole top layer is a mix of Pd and O → almost any bridge
   qualifies as "br_PdO" even though it is chemically an oxide site.
   This is a limitation of the label; the region tag mitigates it but
   does not fully solve it.

3. Oxide top site (`atop_O`) vs metal terrace (`atop_Pd`) are cleanly
   separated by the anchor rule — no ambiguity.

**S2/CH3O idx=217** — MLIP-relaxed: O anchor within 2.6 Å of one Pd and
one O (n_Pd=1, n_O=1 → `br_PdO`, region=`interface`). Selected as
S2's interface representative.

**S2/CH3O idx=496** — MLIP-relaxed: O anchor within 2.6 Å of only one O
(n_Pd=0, n_O=1 → `atop_O`, region=`oxide`). Selected as S2's pure-oxide
representative — this is the only atop_O candidate in the S2 CH₃O pool.

**S2/coads idx=8079** — combo `atop_Pd + atop_Pd`, d_reactive=2.86 Å,
distinct from the dominant `br_PdPd+br_PdPd` cluster at ~3.3 Å. Chemically
different because both adsorbates sit atop coordinatively-unsaturated Pd
sites at the interface rather than sharing a bridge with the terrace Pd.

---

## 6. Spatial / structural dedup

### v4 (in `select_top5_v4.py:96-108, 132-145`)

- E-cluster tolerance: `DELTA_E_DUP = 0.03 eV`
- xy tolerance: `XY_DUP = 1.5 Å`, computed as **PBC-aware
  minimum-image** distance in the in-plane fractional basis
  (`inv @ dx → wrap to [-0.5, 0.5] → back to Cartesian`).
- Duplicate := (|ΔE| < 0.03 eV) **AND** (xy MIC < 1.5 Å) between anchor
  atoms. Site label, region, orientation, RMSD **not consulted.**
- coads dedup: BOTH anchors (C_CO and O_CH3) must satisfy xy MIC <
  1.5 Å simultaneously for the pair to be a duplicate. Site swap (CO/CH3O
  interchange) is **not** treated as equivalent.

### v5 (`select_top5_v5.py:167-186`)

- Same E and xy MIC computation as v4.
- **Additionally requires `site_class` match** before declaring a
  duplicate for singles (`select_top5_v5.py:206-207`).
- Coads: xy MIC on both anchors (same as v4).
- No RMSD, no explicit orientation vector, no local-environment
  fingerprint. Structural fingerprint field from MLIP screener
  (`unique_*.json:fingerprint`) is not consumed by v5.

### Why v4 lost site diversity, and how v5 recovers it

v4 built its dedup entirely on (ΔE, xy) pairs. Symmetry-equivalent
bridges on Pd(100) sit at *different* xy positions (0.5 Å apart in the
in-plane basis) yet correspond to the same chemistry. v4 kept all of
them because they passed the xy>1.5 Å criterion. Result: S1 CH₃O v4
picks are `br_PdPd 5`; S1 CO is `atop_Pd 1 + br_PdPd 4`.

v5 fixes this by:

1. **Adding site_class as a mandatory diversity dimension** — the
   selector prioritises candidates that occupy a new `site_class` or
   `region` bin, subject to an energy cap.
2. **Enforcing a distance-bin dimension for coads** — even if the same
   `br_PdPd+br_PdPd` combo has three near-degenerate candidates, only
   one per (site combo × distance bin) is taken.

Concretely, S1 CH₃O v5 = `br_PdPd 4 + atop_Pd 1` (idx=283 added at
ΔE=0.002 eV to introduce atop_Pd). S1 CO v5 = `atop_Pd 1 + br_PdPd 3 +
4f_4Pd0O 1` (idx=43 added at ΔE=0.046 eV to introduce the 4-fold hollow).

---

## 7. Single-ads shortlist algorithm

Pseudocode of `select_top5_v5.select_single()` (verbatim structure from
lines 189-227):

```
1. pool = characterize_single(sid, ads)    # filters + descriptors + sort by E
2. picks = []
3. for r in pool (sorted by E ascending):
4.     if r's anchor is xy-dup + site-match of any existing pick: skip
5.     if picks == []:                              # slot 0
6.         picks.append(r); reason = "global E min"
7.         continue
8.     if len(picks) >= 5: break
9.     new_site   = (r.site   not in seen_sites)
10.    new_region = (r.region not in seen_regions)
11.    dE = r.E - picks[0].E
12.    if (new_site or new_region) and dE <= 0.5:
13.        accept; reason = "new site/region"
14.    else:
15.        record SKIP reason (diversity but dE>0.5)
16. # 2nd pass: fill remaining slots by lowest E respecting dedup
17. for r in pool:
18.    if len(picks) >= 5: break
19.    if r duplicates any pick: skip
20.    accept; reason = "fill (ΔE=…)"
21. return picks, log
```

**Fixed constants**:
- Target `N_TARGET = 5` (`select_top5_v5.py:28`).
- Energy cap for diversity picks: `dE <= 0.5 eV`.
- xy dedup threshold: `XY_CLUSTER = 1.5 Å` (`select_top5_v5.py:38`).
- E dedup tolerance: `DELTA_E_DUP` was 0.03 eV in v4; v5 uses only
  `site_class + xy_MIC` for dedup, no explicit E-cluster tolerance.
- Global minimum: **always** kept (slot 0).
- Site bin priority: whatever appears first in ascending-E order.
- Region bin priority: same.
- No hard limit on # candidates per site cluster (governed only by xy
  and site_class dedup).
- Fallback: 2nd pass fills by lowest-E after dedup, regardless of site
  diversity — so if the pool is 5 copies of br_PdPd, v5 will produce 5
  br_PdPd (this is what happens on S1 CH₃O beyond slot 1).

### S4/CO worked example

- Pool: 26 valid candidates, `E_min = −756.579 eV`, `E_max = −753.440 eV`
  → ΔE_span = 3.14 eV.
- All 5 lowest-E candidates classify as `physi` (anchor >2.6 Å from any
  Pd or O). Their E values cluster within 0.05 eV of the global min.
- Chemisorbed candidates exist further down the list:
  - `idx 36`: br_PdO, ΔE = 2.54 eV
  - `idx 52`: br_PdO, ΔE = 2.69 eV
  - `idx 91`: atop_Pd, ΔE = 2.69 eV
- These are **rejected by the 0.5 eV cap** in `select_single`. v5 log
  reports "SKIP diversity candidate ΔE=2.54>0.5 eV" (see
  `report.md:89-91`).
- Interpretation: the MLIP finds physisorption to be ~2.5 eV lower than
  any chemisorbed alternative on stoichiometric PdO₂(110). This is
  consistent with the geometric fact that the surface has no
  undercoordinated Pd cus site available for CO chemisorption (all top-
  layer Pd are 6-coordinated). We accept 5 physi picks as v5's S4/CO
  slate.
- **Preliminary signal**: this argues Case D (PdO₂ = DMC-inactive
  reference) from Phase 1 alone. **DFT must confirm**: does PBE-D3+VASPsol
  keep CO at 3 Å from the surface (physisorbed) or does it push CO down
  onto Pd cus and form a chemisorbed state that MLIP missed? A DFT E_bind
  of −0.1 to −0.3 eV would confirm physisorption; anything below −0.6 eV
  would flag an MLIP failure worth investigating.

---

## 8. Coadsorption distance definition

`d_reactive` in the v5 code:

```python
d_reactive = a.get_distance(co_c, me_o, mic=True)   # select_top5_v5.py:157
```

- Between **the C of CO*** and **the O of CH₃O***.
- ASE's `get_distance(..., mic=True)` = periodic minimum-image
  convention, full 3D (not xy-only).
- Measured on the **MLIP-relaxed** geometry, not the initial guess.
- Initial coords are also in `MLIP_phase2/relaxed_SetA.traj` (the traj
  contains only final frames; initial frames are in
  `coads_guide/SetA.traj`).

### Distance bins used in v5 (`select_top5_v5.py:123-127`)

| bin | range (Å) | physical meaning | Phase 1 role | expected DFT outcome |
|---|---|---|---|---|
| `product-like` | <2.1 | C–O bond forming; effectively CH₃OCO* nascent state | drop — not a valid coads endpoint | — (none in pool for any surface) |
| `reactive-close` | 2.1–3.0 | genuine TS1 endpoint distance; CO and CH₃O both anchored on adjacent sites | required for T2.5 NEB endpoint | expected d in DFT: 2.7–3.4 Å, may drift slightly |
| `reactive-loose` | 3.0–4.0 | broader reactive window; matches workplan §9 SetA cutoff 2.1–4.0 Å | complementary reactive endpoint | expected d in DFT: 3.0–4.5 Å |
| `separated` | 4.0–5.0 | too far for concerted TS1, but not quite thermodynamic reference | secondary check for site-availability effects | often collapses to thermo bin in DFT |
| `thermodynamic` | ≥5.0 | pair reference where CO and CH₃O are essentially non-interacting | reference for interaction energy | d in DFT: ~5–7 Å depending on slab size |

Cutoff sourcing: workplan §9 (`docs/DMC_Pd_workplan.md:66`) prescribes
"reactive atom distance (C_CO···O_OCH₃) 2.1–4.0 Å" and "thermodynamic
coadsorption reference ≥ 5.0 Å". The 2.1 Å lower bound = C–O single
bond length ~1.4 Å + a 0.7 Å reactive buffer (i.e. distances shorter
than 2.1 Å indicate the C–O bond is already forming, so the structure
is no longer a coads endpoint). The split at 3.0 Å (close vs loose) is
project convention introduced in v5 to separate TS1-imminent (short)
from ready-but-not-imminent (mid-range) endpoints; there is no
literature-fixed value at 3.0 Å.

---

## 9. Why coads needs distance-bin diversity, not just energy

1. **Global E min ≠ reactive endpoint**. The lowest-E coads structure
   is very often the one where CO and CH₃O sit on well-separated sites
   (each near its individual optimum) with `d_reactive ≥ 5 Å`. E.g.
   S3 coads global min (idx 2408): d_reactive = 5.40 Å = `thermodynamic`.
2. If we only took the lowest-E structures, **none** of the S3/S3b picks
   would be within the reactive window — v4 in fact produced all 5
   S3/S3b picks in the `thermodynamic` bin (5.18–5.40 Å for S3b, similar
   for S3). This is precisely what workplan §T2.1 needs to sample: the
   *reactive* range 2.1–4.0 Å, not the thermodynamic tail.
3. `reactive-close/loose` candidates are the direct initial guesses for
   T2.5 NEB (TS1 = CO* + CH₃O* → CH₃OCO*). Without them, we have no way
   to start the barrier calculation.
4. `thermodynamic` candidates are needed to compute the coadsorption
   interaction energy against isolated coads on the same surface — the
   `ΔE_coads_interaction` in the workplan §T1.18 descriptor map.
5. **`product-like` (<2.1 Å)** would be Phase 2 (post-TS1) or dissociation
   artifacts. Currently no candidates fall in this bin — this is
   expected: LBFGS relaxation at Phase 1 does not spontaneously form the
   C–O bond because there is no dedicated reactive pathway. If Phase 2
   NEB or dimer finds any, they will show up as new picks in the T2
   pipeline.
6. **Site combination × d_reactive**: two coads with the same combo
   (`br_PdPd+br_PdPd`) but different d_reactive (say 3.1 Å vs 5.4 Å)
   are chemically distinct — one is a TS1 initial guess, the other is a
   thermodynamic reference. Selecting by combo alone would miss this.
7. **Orientation, alignment, steric overlap** — currently **NOT** used
   in the selection code. The MLIP dedup fingerprint has intra-adsorbate
   distances only. Post-DFT, we should verify orientation of the C=O
   axis relative to the O–CH₃ vector, since a reactive pair with
   misaligned axes may still not couple in DFT.

### What is / is not consumed by the current selector

| descriptor | selector uses? | Phase 1 role | plan for DFT check |
|---|---|---|---|
| C_CO–O_CH3O distance | **yes** (bins) | primary reactive descriptor | recompute post-DFT |
| CO molecular axis / tilt | **no** | not recorded | measure in CONTCAR |
| CH₃O orientation (C–O tilt vs surface normal) | **no** | not recorded | measure in CONTCAR |
| adsorbate height (anchor z − z_top) | recorded but not used in selection | descriptor for sanity checks | flag if <1.5 Å (compressed) or >3 Å (physi) |
| steric overlap check | dedup fingerprint rejects fragmented / collapsed (`refilter_phase2_geometry.geometry_valid`) | pre-selection filter | — |
| site combination | **yes** (dedup + diversity) | see §4 | confirm post-DFT |
| xy anchor separation | **yes** (dedup) | see §6 | — |
| surface reconstruction | **no** — MLIP screener does not detect it | none | inspect CONTCAR |
| CO₂-like collapse | **yes** (`valid_CO/coads` rejects 2nd O within 2.0 Å) | pre-selection filter | if DFT produces C–O_lattice bond →2 Å, side-path |

---

## 10. Per-group selection results

See **`audit/selection_candidate_rationale.csv`** — 86 rows, one per
DFT candidate, with `origin` (v4-existing / v5-new), `priority`,
`E_MLIP`, `E_bind_MLIP`, `dE_from_global`, `site`, `site_combo`,
`region`, `d_reactive`, `d_bin`, human-readable `selection_reason`, and
a `dft_hypothesis` column stating what each DFT calculation is expected
to confirm or refute.

Selection-reason vocabulary (assigned per candidate, not "diversity"):

- `global MLIP minimum`
- `lowest-E atop_Pd / br_PdPd / … representative`
- `PdO/Pd interface representative`
- `pure oxide top-O site representative`
- `4-fold hollow diversity`
- `reactive-close endpoint candidate`
- `reactive-loose endpoint candidate`
- `thermodynamic coadsorption reference`
- `low-energy physisorbed alternative`
- `high-energy diagnostic reactive candidate` (S3b/coads idx=2051)
- `top-5 MLIP E + xy dedup (v4-baseline)` — for v4-existing picks that
  were made before the diversity rule was in place

`dft_hypothesis` gives a concrete falsifiable prediction, e.g. "d_reactive
in DFT stays below 4 Å after relaxation" or "S4/CO physi remains
non-bonded with |E_bind| < 0.3 eV".

---

## 11. v4 vs v5

| aspect | v4 (`select_top5_v4.py`) | v5 (`select_top5_v5.py`) |
|---|---|---|
| ranking value | `E_MLIP` | `E_MLIP` (identical inside a group) |
| E-cluster tolerance | 0.03 eV | not used explicitly (deferred to site_class dedup) |
| site classification used | **no** | **yes** — anchor + 2.6 Å neighbor rule |
| region classification | **no** | **yes** — metal / oxide / interface / physi |
| xy dedup | PBC-aware MIC, 1.5 Å | same |
| structural dedup | not used | not used |
| distance bins (coads) | **not used** | 5 bins per §8 |
| energy cap for diversity | none (E-cluster + xy only) | 0.5 eV (single) / 0.8 eV (coads) |
| fallback | keep filling by lowest-E if dedup keeps rejecting | 2nd pass, lowest-E after dedup |
| candidate count | 5 per (sid × ads) | same |
| known limitation | site diversity was accidental | site diversity enforced, but ΔE cap can exclude chemically valuable outliers (mitigated by MUST-diagnostic tag) |

**What v4 missed and v5 fixed** — concrete cases in the current shortlist:

- **S1/CO**: v4 = atop_Pd 1 + br_PdPd 4 (i.e. one accidental atop pick,
  four bridge duplicates). v5 = atop_Pd 1 + br_PdPd 3 + **4f_4Pd0O 1**
  (idx=43 added at ΔE=0.046 eV — a 4-fold hollow site distinct from
  bridge).
- **S1/CH3O**: v4 = 5 br_PdPd. v5 = br_PdPd 4 + **atop_Pd 1** (idx=283
  added at ΔE=0.002 eV — nearly free).
- **S2/coads**: v4 = 5 br_PdPd+br_PdPd at d_reactive ≈ 3.3 Å.
  v5 = adds **reactive-close** (idx=8079, atop_Pd+atop_Pd at 2.86 Å) and
  **thermodynamic** (idx=3633, 5.19 Å) representatives.
- **S3/coads**: v4 = 5 thermodynamic (all d_reactive ≥ 5.36 Å) with **no
  reactive-close/loose**. v5 = adds reactive-close idx=3481 (2.99 Å) and
  reactive-loose idx=5161 (3.25 Å) as MUST — without these the T2.5 NEB
  cannot start on S3.
- **S3b/coads**: v4 = 5 thermodynamic. v5 = adds reactive-close idx=2051
  (2.90 Å, ΔE=0.55 eV, **diagnostic**) and reactive-loose idx=2754
  (3.66 Å, ΔE=0.12 eV).
- **S4/CO**: v4 = 5 physi. v5 = still 5 physi (see §7 — no chemisorbed
  candidate within ΔE cap).
- **S4/coads**: **intentionally excluded** in both v4 and v5 (researcher
  decision).

---

## 12. Selection robustness / sensitivity

See **`audit/selection_sensitivity.csv`** for the row-level table. For
each of the 16 v5-new candidates the site label is re-computed at
neighbor cutoffs 2.4 / 2.6 / 2.8 Å, and (for coads) the distance-bin
label is re-computed at reactive-close boundary 2.8 / 3.0 / 3.2 Å and
reactive-loose upper bound 3.8 / 4.0 / 4.2 Å.

Categorisation:

- `robust`: same site label at all 3 cutoffs AND (for coads) same
  distance bin at all bin-boundary settings
- `threshold-sensitive-label`: site label changes across cutoffs but
  candidate remains in the same broad chemistry (e.g. atop_Pd ↔ br_PdPd
  as cutoff widens)
- `threshold-sensitive-bin`: distance bin flips at bin-boundary changes
  (e.g. 2.99 Å candidate crosses the reactive-close/loose divide when
  boundary moves from 3.0 to 2.8 Å)
- `diagnostic`: candidate was included precisely because DFT needs to
  arbitrate, and label instability underlines that

Aggregate counts appear in the row `robustness` column. Single-ads
`energy_cap` at 0.3 vs 0.5 eV: candidates with `dE_from_global > 0.3`
would be dropped at the tighter cap; those with `dE > 0.5` are already
excluded. `xy_dedup` at 1.0 or 2.0 Å changes the number of candidates
that survive dedup but has minimal effect on the 16-approved subset (all
16 were beyond the 1.5 Å tolerance from any v4 pick — see
`xy_dist_to_nearest_v4` in `proposed_additions_only.csv`, min value 0.21
Å for S3/coads 3481 which is very close to a v4 pick geometrically but
in a different d_reactive bin).

---

## 13. MLIP selection limits and the role of DFT

MLIP (MACE+D3) selection is **necessary but not sufficient**. Explicit
open questions that DFT must resolve:

1. Absolute `E_bind_MLIP` values are not the final answer. MLIP has a
   ~50 meV systematic bias on adsorption energies (project MLIP
   benchmark). The final descriptor map in T1.19 must use DFT values.
2. Site labels are diversity descriptors, not final claims. A candidate
   labelled `br_PdO` at MLIP-relaxed geometry may relax to `atop_Pd`
   after DFT allows the surface to reconstruct.
3. `reactive-close` at MLIP may drift to `separated` under DFT if the
   real Born-Oppenheimer PES has a wider well than MACE thinks.
4. Conversely, some `thermodynamic` MLIP picks may re-arrange to
   `reactive-close` under DFT if the MLIP missed a coupling term.
5. MLIP energy ordering within a group can flip under DFT for candidates
   within ΔE < ~0.1 eV. Slot 0 = MLIP global min is not guaranteed to be
   the DFT global min.
6. CO₂-like or dissociation calls made by `valid_*` are geometric
   heuristics. DFT may show that a "CO₂-like" candidate actually kept
   the CO bond and has a strong side-path signal (or vice versa).

**Per-candidate DFT checks required** (populated in
`selection_candidate_rationale.csv:dft_hypothesis`):

| descriptor | pass criterion |
|---|---|
| convergence | OUTCAR `reached required accuracy`; max force < 0.03 eV/Å |
| final E_bind vs MLIP | |ΔE_DFT − ΔE_MLIP| < 0.2 eV within group (loose) |
| final adsorption site | recompute anchor + 2.6 Å neighbor rule on CONTCAR; compare to metadata |
| final anchor–surface distance | must match site (chemisorbed: 1.8–2.5 Å; physi: >2.8 Å) |
| final `d_reactive` (coads) | should fall in the same bin ±0.5 Å; if crosses a bin boundary, flag |
| adsorbate dissociation | check C–H count, C–O bond length, second O within 2.0 Å of C |
| CO₂-like conversion | if 2 O atoms within 2.0 Å of C → side-path (T2.7) candidate |
| collapse to same minimum | different candidates converging to identical CONTCAR ⇒ dedup post-DFT |
| dipole moment | for interpretation of ΔG_MeOH(U) formula |

---

## 14. Bottom-line questions

1. **Is the current shortlist selected on binding energy alone?**
   No. It is selected on `E_MLIP` ranking + site/region diversity + (for
   coads) distance-bin diversity, with energy caps of 0.5 eV (singles) /
   0.8 eV (coads). Binding energy is used only for visualisation and
   cross-surface interpretation.

2. **Where is `E_MLIP` used, where is `E_bind_MLIP` used?**
   `E_MLIP` = selector input (within-group ranking is invariant to any
   constant shift). `E_bind_MLIP` = violin plot y-axis and descriptor
   map input (cross-surface comparison requires reference subtraction).

3. **Why site diversity?**
   Because near-degenerate MLIP energies do not distinguish physically
   distinct sites. Workplan §T1.15 requires "거리 bin별 대표 구조 포함";
   v4 lost this requirement because it dedup'd on xy only.

4. **Why distance-bin diversity for coads?**
   Because the lowest-E coads structures are almost always
   thermodynamic (well-separated), so pure energy ranking misses the
   reactive endpoints needed for T2.5 TS1 NEB.

5. **Why not just take the lowest E?**
   Because a single minimum is not a descriptor map — we cannot resolve
   the sensitivity of the reaction path to site, orientation, or
   pairing distance from one point in phase space.

6. **Why not include arbitrarily high-E site-diverse candidates?**
   Because MLIP has ~50 meV systematic error and DFT reordering
   typically happens within ~0.1–0.3 eV. A ΔE > 0.5 eV candidate is more
   likely to reflect a MLIP artifact (or a genuinely unfavourable state)
   than a genuine physical minimum worth spending DFT wallclock on. The
   only exceptions are diagnostic candidates whose inclusion answers a
   specific scientific question (S3b/coads 2051).

7. **What uncertainty does v4-70 + v5-add-16 = 86 candidates reduce?**
   - v4 (70) resolves the lowest-E cluster on each surface, and gives
     multiple xy-distinct picks for redundancy.
   - v5-add (16) fills three gaps:
     (a) reactive-close/loose coads on S3/S3b (T2.5 NEB precondition),
     (b) distinct site combos on S1/S2 coads (interface/oxide),
     (c) alternative single-ads sites on S1/S2/S3/S3b/S4 (region
     diversity for descriptor map).
   - Overall, 86 is chosen to keep DFT cost linear with the number of
     independent scientific questions we care about (~10 hypotheses
     from §10), not the geometric parameter space (~10⁴ candidates).

8. **What decides the final 2–3 winners after DFT?**
   - **Convergence + physical sanity** (site remains, adsorbate intact,
     no unexpected CO₂ formation).
   - **Global DFT minimum** in each (sid, ads) group → primary reference.
   - **One reactive-close coads** per surface → TS1 NEB endpoint.
   - **Descriptor map placement** — surfaces are grouped into Case A–D;
     the 2–3 winners per Case get promoted to T1.20 → T2.
   - **Deduplication post-DFT** — if 3 v4 near-degenerate candidates
     collapse to the same CONTCAR under DFT, they count as 1.

9. **What's the biggest remaining risk in the current selection?**
   - **S4/CO physi cluster**. v5's ΔE cap correctly excluded chemisorbed
     candidates at ΔE > 2.5 eV, but if DFT shows one of those becomes
     a chemisorbed minimum after full ionic relaxation, we would have
     under-sampled S4/CO. Mitigation: re-run MLIP with a broader ΔE cap
     if the physi results are trivially weak.
   - **S3b/coads idx=2051** (ΔE = 0.55 eV) is a diagnostic candidate;
     if DFT rejects it, we lose the only S3b reactive-close endpoint.
   - **Site-label sensitivity at the 2.6 Å cutoff** — a few interface
     candidates flip labels at 2.8 Å; DFT will decide which side of the
     boundary they belong to.

10. **When does DFT refute MLIP-based selection?**
    - Any candidate whose DFT-relaxed site differs from its MLIP label
      by more than a "physi ↔ chemisorbed" transition → MLIP mis-classify.
    - Any candidate where `|E_bind_DFT − E_bind_MLIP| > 0.3 eV` → MLIP
      systematic bias breaks the ranking assumption.
    - Any coads pick where `|d_reactive_DFT − d_reactive_MLIP| > 0.5 Å`
      and the bin changes → MLIP geometry unreliable for T2.
    - Any candidate reported as `converged` by MLIP but not by DFT →
      MLIP LBFGS `fmax=0.05` was too loose for that geometry.

---

## Files consulted (as of audit date)

- `scripts/run_mace_phase1.py` (L54-199)
- `scripts/run_mace_phase2.py` (L82-190)
- `scripts/refilter_phase2_geometry.py` (L51-99, L102-190)
- `scripts/compute_mace_d3_references.py` (full)
- `scripts/select_top5_v4.py` (full)
- `scripts/select_top5_v5.py` (full)
- `scripts/refine_v5_additions.py` (full)
- `scripts/setup_v5add16.py` (full)
- `scripts/build_violin_v5add.py` (full)
- `scripts/audit_selection.py` (this audit)
- `docs/DMC_Pd_workplan.md` §T1.14–T1.20
- `docs/DMC_Pd_package_guideline.md` §T1.14–T1.15
- `calculations/G3_adsorption/DFT_shortlist_v3/summary.json`
- `calculations/G3_adsorption/DFT_shortlist_v5/*.csv, report.md, addendum.md`
- `calculations/T1_16_DFT_L2_v5add16/manifest.csv, combined_summary.csv`
- `calculations/G3_adsorption/mace_d3_references.json`
