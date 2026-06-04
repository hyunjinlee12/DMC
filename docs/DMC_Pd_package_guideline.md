# 연구단계별 패키지 사용 가이드라인 — Pd/PdO/PdO₂ DMC DFT

작성일: 2026-06-03
근거: 워크플랜(DMC_Pd_workplan.md) + 참고문헌 요약(DMC_Pd_references_summary.md)
목적: 각 작업(T-ID)에서 어떤 패키지를, 어떤 옵션으로, 어떤 순서로 쓰는지 명시.

## 0. 소프트웨어 스택 한눈에

| 역할 | 패키지 | 비고 |
|------|--------|------|
| 구조 생성·조작·I/O 허브 | ASE | 모든 단계의 공통 backbone (ase.io, ase.build, constraints) |
| Bulk/slab/표면 생성 | pymatgen (+ ASE) | SlabGenerator, Wulff, oxide termination 열거 |
| 흡착 구조 생성 | AutoAdsorbate (+ RDKit) | SMILES fragment → site 배치, co-adsorption custom wrapper |
| 사전 relax·ranking | MACE (foundation MLIP) | mace_mp ASE calculator, 필요 시 fine-tuning |
| DFT 판정 | VASP (+ VASPsol) | PBE+D3, asymmetric slab, implicit solvation |
| TS 탐색 | VTST tools (NEB/dimer) + ASE-NEB | CI-NEB IMAGES=5 |
| 열역학/자유에너지 | ASE thermochem + 자체 후처리 스크립트 | ZPE, CHE(U), ab initio thermodynamics |
| 시각화·검증 | ASE GUI / OVITO / VESTA | 구조 sanity check |

설치 권장(별도 env): ase, pymatgen, rdkit, autoadsorbate, mace-torch, numpy/scipy/pandas. DFT는 클러스터 VASP, MLIP는 GPU 노드 권장.

---

# PHASE 1

## T1.1–T1.4 — Bulk 최적화 (Pd / PdO / PdO₂)

패키지: pymatgen(또는 Materials Project 구조 가져오기) → ASE/VASP.

1. 초기 구조 확보: fcc Pd, tetragonal PdO, rutile/hydrophilite-like PdO₂ — pymatgen Structure 또는 MP API에서 가져와 ASE로 변환.
2. VASP relaxation (ISIF=3, IBRION=2, NSW=200), INCAR는 워크플랜 §P1-A 그대로.
   - Smearing: Pd ISMEAR=1/0.10, PdO·PdO₂ ISMEAR=0/0.05.
3. 수렴 테스트는 ASE/pymatgen로 ENCUT·k-mesh sweep 자동화, 실험 격자상수와 비교.

산출물: 최적 CONTCAR 3종 + bulk reference 에너지. 이후 모든 slab의 기준.

## T1.5–T1.9 — Clean slab 제작·검증

패키지: pymatgen SlabGenerator (절단·termination 열거) → ASE (vacuum·constraint·정리) → VASP relax.

- S1 Pd(100): `SlabGenerator(bulk_Pd, (1,0,0), ...)` → p(4×4), 4–5 layers. ASE `FixAtoms`로 bottom 2 layers 고정, `center(vacuum=10, axis=2)`로 vacuum 20 Å(양쪽 합) 확보 후 asymmetric로 조정.
- S2 1 ML PdO(101)/Pd(100): Pd(100)-(√5×√5)R27°-O 모델 기반 (참고문헌 B1). PdO(101) 단층을 Pd(100) 위에 epitaxial로 올리고 가능하면 2×2 supercell. 수동 구성 + ASE로 격자 정합·strain 확인.
- S3 bulk PdO(100)-PdO / O-rich PdO(100): pymatgen으로 PdO(100) termination 열거, O-rich top termination 선택.
- S4 PdO₂(110): stoichiometric rutile-like termination.

검증(T1.9): ASE로 rumpling·Pd–O bond length·coordination·slab dipole 계산. O-rich/PdO₂ termination 안정성은 ab initio thermodynamics(surface energy vs μ_O) — 자체 스크립트 + bulk reference.

공통 규칙: asymmetric, vacuum 20 Å, LDIPOL=.TRUE., IDIPOL=3, ISYM=0, INCAR는 워크플랜 §P1-B/§1.

팁: pymatgen이 만든 slab은 항상 ASE GUI/VESTA로 육안 확인 — termination·층수·진공이 의도대로인지.

## T1.10–T1.14 — 흡착 sampling (AutoAdsorbate + MLIP)

패키지: AutoAdsorbate(+RDKit) → MACE(mace_mp).

단일 흡착 (T1.11 CO*, T1.12 CH₃O*):

```python
from autoadsorbate import Surface, Fragment
surface = Surface(slab)  # clean slab의 site map

# CO* — top/bridge/hollow 모두: mode='all'
co = Fragment('Cl[C-]#[O+]', to_initialize=1)
co_structs = surface.get_populated_sites(
    co, mode='all', sample_rotation=False,
    conformers_per_site_cap=1, overlap_thr=1.25)

# CH3O* — 회전 자유도 큼: rotation + conformer 다수
och3 = Fragment('ClOC', to_initialize=20)
och3_structs = surface.get_populated_sites(
    och3, mode='all', sample_rotation=True,
    conformers_per_site_cap=3, overlap_thr=1.25)
```

반드시 `mode='all'` (heuristic은 top site만). 후보 개수 목표는 워크플랜 표(표면당 100–300).

Co-adsorption (T1.13): `get_populated_sites()` 단독으로 부족 → custom wrapper.

```python
def generate_coadsorption(surface, co_fragment, och3_fragment, out_path,
                          pair_dist_min=2.0, pair_dist_max=4.5,
                          co_conformers=(0,), och3_conformers=(0,0.25,0.5,0.75),
                          max_per_pair=2):
    # 1) clean site map에서 site_i에 CO, site_j에 CH3O 순차 배치
    # 2) site_i != site_j, site-center distance로 1차 필터
    # 3) C_CO–O_CH3O distance filter
    # 4) steric overlap filter (heavy<1.6Å, H-heavy<1.1–1.2Å reject)
    # 5) Set A(reactive)/Set B(thermo)/side-path로 분류
    ...
```

분류 cutoff(워크플랜 §P1-C): site-pair primary ≤4.5 Å, reactive atom C_CO···O_OCH₃ 2.1–4.0 Å, thermo ≥5.0 Å.

MLIP relaxation·ranking (T1.14):

```python
from mace.calculators import mace_mp
calc = mace_mp(model="medium", dispersion=True, default_dtype="float64")
# 각 후보 atoms.calc = calc 후 ASE BFGS/FIRE로 relax → energy로 정렬
```

MLIP는 ranking 전용(절대값 신뢰 X). 우선순위: ① 낮은 총에너지 ② C_CO–O_CH3O 2.0–3.5 Å ③ 서로 다른 기능성 site ④ S2는 interface pair 우선 ⑤ O-rich/PdO₂의 CO₂-like 붕괴는 side-path로 보관.

(선택) 정량 필요 시 표면당 5개 DFT single-point로 MACE fine-tuning.

## T1.15–T1.20 — DFT 흡착에너지·descriptor map

패키지: VASP(+VASPsol) → ASE/pandas 후처리.

- T1.15 shortlist: 거리 bin별 대표 포함, 표면·화학종별 2–4개(워크플랜 표).
- T1.16 Level 1: PBE-D3 vacuum relaxation (clean_slab_relax INCAR, EDIFFG=-0.03).
- T1.17 Level 2: VASPsol final energy/relax (LSOL=.TRUE., EB_K=32.6, TAU=0).
- T1.18 흡착에너지:
  - ΔG_CO* = G_slab+CO − G_slab − μ_CO
  - ΔG_CH3O*^rad = G_slab+CH3O − G_slab − G_CH3O
  - ΔG_CH3O*^MeOH(U) = G_slab+CH3O + ½G_H2 − G_slab − G_CH3OH − eU (CHE)
  - ZPE/열보정은 ASE Thermochemistry (adsorbate=harmonic, gas=ideal-gas).
- T1.19 descriptor map: pandas로 ΔG_CO* vs ΔG_CH3O*^MeOH(U) 산점도(matplotlib).
- T1.20: Case A–D 프레임으로 Phase 2 표면 선정.

**게이트 G3:** descriptor map 완성 + 후보 표면 확정.

---

# PHASE 2 — DMC Gibbs free-energy profile

## T2.1–T2.4 — reactive pair·intermediate sampling·endpoint 최적화

패키지: AutoAdsorbate custom wrapper → MACE → VASP(+VASPsol).

- T2.1 reactive pair: Set A 기준(C_CO···O_OCH₃ 2.1–4.0 Å, site-pair ≤4.5 Å)으로 좁혀 생성.
- T2.2 intermediates: CH₃OCO*(Fragment로 acetyl-유사 SMILES 구성), CH₃OCO*+CH₃O*, DMC* 후보 생성.
- T2.3 재분류: MLIP relax 후 C···O 거리 bin(<1.6 product / 1.6–2.3 TS-like / 2.3–4.0 reactive / ≥5.0 thermo)으로 DFT shortlist.
- T2.4 endpoint: PBE-D3+VASPsol 최적화.

## T2.5–T2.6 — TS1 / TS2

패키지: VTST(VASP CI-NEB) 또는 ASE-NEB, 보조로 dimer.

초기·최종 endpoint(T2.4)로 NEB image 생성(ASE interpolate 또는 IDPP). INCAR: IMAGES=5, SPRING=-5, LCLIMB=.TRUE., IBRION=3, POTIM=0, NSW=300, VASPsol 포함(워크플랜 §Phase2).
TS guess 거리: TS1 C···O 1.7–2.3 Å, TS2 1.8–2.4 Å. 수렴 후 frequency로 1개 허수진동 확인. 어려운 안장점은 dimer로 보완.

## T2.7 — side-path 검증

CO + O_lattice → CO₂ + V_O (C_CO···O_lattice ≤3.6 Å), methanol oxidation 경로. 동일 NEB/endpoint 방식.

## T2.8 — Gibbs free-energy profile

ASE thermochem으로 ΔG(U) 보정, CHE로 전위 의존 step 처리. matplotlib로 표면별 profile 중첩.
벤치마크 비교: Angew(Shi 2024) pure Pd TS1/TS2 = 1.08/0.85 eV, Pd₃Cu = 0.86/0.79 eV.

**게이트 G4:** 표면별 DMC profile + side-path 비교 → "anodic 산화가 DMC 촉진 vs 억제" 결론.

---

## 4. 단계→패키지 매핑 요약표

| 단계 | 핵심 패키지 | 출력 |
|------|-------------|------|
| T1.1–1.4 bulk | pymatgen/MP → VASP(PBE+D3) | 최적 bulk 3종 |
| T1.5–1.9 slab | pymatgen SlabGenerator → ASE → VASP | clean slab + 검증 |
| T1.10–1.13 sampling | AutoAdsorbate(+RDKit) | 후보 풀 |
| T1.14 ranking | MACE mace_mp + ASE opt | 정렬된 후보 |
| T1.15–1.20 흡착 DFT | VASP+VASPsol → ASE/pandas | descriptor map |
| T2.1–2.3 P2 sampling | AutoAdsorbate → MACE | endpoint shortlist |
| T2.4 endpoint | VASP+VASPsol | relaxed intermediates |
| T2.5–2.7 TS·side-path | VTST CI-NEB / ASE-NEB / dimer | barrier |
| T2.8 profile | ASE thermochem + matplotlib | DMC Gibbs profile |

## 5. 실무 주의사항

- ASE를 단일 I/O 표준으로: pymatgen↔ASE 변환 시 AseAtomsAdaptor, 좌표·셀 손실 없는지 확인.
- AutoAdsorbate mode: CO*/CH₃O*는 항상 mode='all'. co-adsorption은 내장 함수 말고 custom wrapper.
- MLIP는 ranking 전용: 절대 에너지로 결론 내지 말 것. 거리 bin 대표 구조까지 DFT 교차검증. 정량 필요 시 fine-tuning.
- VASPsol: vacuum relax(Level 1)에서 수렴 잡고 solvation(Level 2) 적용 — 처음부터 LSOL 켜면 수렴 불안정.
- asymmetric slab dipole: LDIPOL/IDIPOL=3 필수, vacuum 20 Å 유지, bottom 고정 일관성.
- ISYM=0: 흡착·결함 있는 slab은 대칭 끄기.
- 재현성: INCAR/KPOINTS/POTCAR·패키지 버전(특히 MACE 모델 체크포인트)을 결과와 함께 기록.

## 출처

- MACE foundation model (mace_mp ASE calculator) — https://mace-docs.readthedocs.io/en/latest/guide/foundation_models.html
- MACE fine-tuning tutorial — https://arxiv.org/html/2506.21935v2
- Pd(100) √5×√5 surface oxide = PdO(101)/Pd(100) — https://arxiv.org/abs/cond-mat/0304107
- Shi et al. Angew 2024 (TS barrier 벤치마크) — https://onlinelibrary.wiley.com/doi/10.1002/anie.202401311
- (재확인) AutoAdsorbate / ACS Catal 10.1021/acscatal.5c06553 — https://pubs.acs.org/doi/10.1021/acscatal.5c06553
- INCAR/VASPsol/cutoff 파라미터: 프로젝트 대화록 pd_dmc_conversation_export.md
