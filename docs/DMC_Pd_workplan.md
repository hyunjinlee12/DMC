# Pd/PdO/PdO₂ 표면 DMC Formation DFT 연구 — 프로젝트 워크플랜

작성일: 2026-06-03
근거: pd_dmc_conversation_export.md (대화록 §1–§10, 부록 A–C)
형식: 작업 분해(WBS) 중심 — 절대 일정 없음. 각 작업은 ID·선행조건·산출물·완료조건으로 정의.

## 0. 프로젝트 개요

목표 시스템: Pd nanocube {100} facet이 anodic methanol carbonylation 조건에서 산화·재구성될 때, CO*와 CH₃O* 흡착 균형이 어떻게 변하고 이 변화가 DMC carbonylation 경로를 촉진하는지(또는 methanol/CO oxidation side-path로 밀어내는지) DFT로 검증.

반응: 2 CH₃OH + CO − 2e⁻ → (CH₃O)₂CO + 2H⁺

최종 후보 표면 (확정, 부록 B):

| ID | 표면 모델 | 산화 상태 | 역할 |
|----|-----------|-----------|------|
| S1 | Pd(100) | Pd⁰ | metallic nanocube facet baseline |
| S2 | 1 ML PdO(101)/Pd(100) | Pd²⁺-rich surface oxide / Pd–PdO interface | 부분 산화 bifunctional interface |
| S3 | bulk PdO(100)-PdO / O-rich PdO(100) | bulk-like Pd²⁺ oxide | thermodynamic / O-rich limit |
| S4 | PdO₂(110) | Pd⁴⁺ reference | high-valence oxide reference |

Surface model hierarchy: S1 → S2 → S3 → S4 (산화 진행 축)

핵심 워크플로우: PBE+D3 bulk 최적화 → asymmetric slab 제작(vacuum 20 Å) → foundation MLIP + heuristic(AutoAdsorbate) sampling → MLIP ranking → 표면·화학종별 2–3개 DFT 후보 → PBE-D3(+VASPsol) 에너지 경향성 → DMC Gibbs free-energy profile.

## 1. 전체 계산 레벨 (전 작업 공통 규칙)

- 모든 slab은 asymmetric (adsorption top side only, bottom clean/fixed).
- Vacuum 20 Å, bottom 2 layers(또는 하단 30–40% 원자) 고정.
- Bulk optimization: PBE+D3 (IVDW=12).
- Dipole correction: LDIPOL=.TRUE., IDIPOL=3.
- Adsorption 후보 생성/pre-ranking: foundation MLIP + heuristic adsorption sampling.
- 최종 에너지 경향성: VASP PBE+D3, 필요 시 VASPsol (LSOL=.TRUE., EB_K=32.6, TAU=0).
- ENCUT=520, PREC=Accurate, LASPH=.TRUE., ADDGRID=.TRUE., ISPIN=2.

### 1-1. Bulk 구조 레퍼런스

| Material | Phase | Space Group | MP-ID | ICSD | k-mesh (adopted) |
|----------|-------|-------------|-------|------|-------------------|
| Pd | fcc | Fm-3m (#225) | mp-2 | — | 12×12×12 |
| PdO | tetragonal | P4₂/mmc (#131) | mp-1336 | — | 8×8×6 |
| PdO₂ | rutile | P4₂/mnm (#136) | mp-1018886 | 647283 | 6×6×8 |

### 1-2. 물질별 INCAR 차이

| Tag | Pd (금속) | PdO · PdO₂ (산화물) | 비고 |
|-----|-----------|---------------------|------|
| ISMEAR | 1 | 0 | 금속 = Methfessel-Paxton; 산화물 = Gaussian |
| SIGMA | 0.10 | 0.05 | 산화물은 band gap → 좁은 smearing |
| ISIF | 3 (bulk) / 2 (slab) | 3 (bulk) / 2 (slab) | bulk: full cell relax; slab: ionic only |
| IVDW | 12 | 12 | D3-BJ (전 물질 공통) |
| LDIPOL | — (bulk) / .TRUE. (slab) | — (bulk) / .TRUE. (slab) | asymmetric slab 전용, IDIPOL=3 |

---

# PHASE 1 — Surface model validation + Adsorption map

목적: 각 표면을 검증하고 CO*·CH₃O* (+ coadsorption) 흡착 에너지 맵을 만들어 Phase 2로 넘길 표면을 선정.

## P1-A. Bulk 구조 최적화

| ID | 작업 | 선행 | 산출물 / 완료조건 |
|----|------|------|-------------------|
| T1.1 | fcc Pd bulk 최적화 (PBE+D3, ISIF=3, ISMEAR=1/SIGMA=0.10) | — | 최적 격자상수, CONTCAR |
| T1.2 | tetragonal PdO bulk 최적화 (ISMEAR=0/SIGMA=0.05) | — | 최적 구조, CONTCAR |
| T1.3 | rutile/hydrophilite-like PdO₂ bulk 최적화 (ISMEAR=0/SIGMA=0.05) | — | 최적 구조, CONTCAR |
| T1.4 | k-point / ENCUT 수렴 테스트, 실험 격자상수 대비 검증 | T1.1–T1.3 | 수렴 리포트, bulk reference 에너지 |

**게이트 G1:** 3개 bulk 구조가 수렴·검증됨 → slab 제작 착수 가능.

## P1-B. Clean slab 제작 및 검증

| ID | 작업 | 선행 | 산출물 / 완료조건 |
|----|------|------|-------------------|
| T1.5 | S1 Pd(100) slab — p(4×4), 4–5 layers, bottom 2 layers fixed | T1.1 | relaxed clean slab |
| T1.6 | S2 1 ML PdO(101)/Pd(100) — (√5×√5)R27°-O 기반, 가능 시 2×2 supercell | T1.1, T1.2 | relaxed interface slab |
| T1.7 | S3 bulk PdO(100)-PdO + O-rich PdO(100) termination | T1.2 | relaxed slab(들), O-rich top termination |
| T1.8 | S4 PdO₂(110) stoichiometric (rutile-like/hydrophilite) | T1.3 | relaxed clean slab |
| T1.9 | Clean slab 검증: surface rumpling, Pd–O bond length·coordination, slab dipole, O-rich termination stability window, PdO₂(110) stoichiometric termination stability | T1.5–T1.8 | 검증 표/리포트 |

**게이트 G2:** 4개 표면 clean slab 검증 완료 → adsorption sampling 착수 가능.

## P1-C. Adsorption sampling (heuristic + MLIP)

AutoAdsorbate 사용. 단일 흡착은 heuristic로 바로 생성하되 CO*/CH₃O*는 top/bridge/hollow를 모두 봐야 하므로 `mode='all'`. Co-adsorption은 `get_populated_sites()`만으로 부족 → clean site map을 먼저 만들고 CO·CH₃O를 site-pair에 순차 배치하는 custom wrapper 사용.

| ID | 작업 | 선행 | 산출물 / 완료조건 |
|----|------|------|-------------------|
| T1.10 | clean slab별 site map 생성 (custom wrapper 기반) | T1.9 | site map per surface |
| T1.11 | CO* 후보 생성 — `Fragment('Cl[C-]#[O+]')`, mode='all', overlap_thr=1.25 | T1.10 | 표면당 100–300 구조 |
| T1.12 | CH₃O* 후보 생성 — `Fragment('ClOC')`, mode='all', sample_rotation=True, conformers≤3 | T1.10 | 표면당 100–300 구조 |
| T1.13 | CO*+CH₃O* co-adsorption 생성 (Set A reactive / Set B thermodynamic) | T1.11, T1.12 | 표면당 150–500 구조 |
| T1.14 | foundation MLIP relaxation + ranking | T1.11–T1.13 | 에너지순 정렬 후보 풀 |

Co-adsorption 생성 절차: clean site map → CO를 site_i 배치 → CH₃O를 site_j 배치 → overlap filter → C_CO–O_CH3O distance filter → reactive/thermodynamic/side-path 분류.

**Cutoff 기준 (§9):**

| 목적 | cutoff |
|------|--------|
| CO*+CH₃O* site-pair primary | ≤ 4.5 Å |
| site-pair loose | 4.5–5.5 Å |
| reactive atom distance (C_CO···O_OCH₃) | 2.1–4.0 Å |
| TS guess (C···O) | 1.7–2.3 Å |
| thermodynamic coadsorption reference | reactive atom distance ≥ 5.0 Å |
| steric reject (heavy-heavy / H-heavy) | < 1.6 Å / < 1.1–1.2 Å |

MLIP ranking 우선순위: ① coadsorption total energy 낮은 구조 ② C_CO–O_CH3O = 2.0–3.5 Å ③ CO·CH₃O가 서로 다른 기능성 site ④ S2에서는 interface pair 우선 ⑤ O-rich PdO/PdO₂에서 CO₂-like로 무너진 구조는 side-path로 보관.

## P1-D. DFT 흡착 에너지 계산 및 경향성 분석

DFT shortlist 개수 (표면별):

| 항목 | S1 Pd(100) | S2 PdO(101)/Pd(100) | S3 PdO(100) | S4 PdO₂(110) |
|------|-----------|---------------------|-------------|--------------|
| CO* DFT 후보 | 3 | 3–4 | 3–4 | 3–4 |
| CH₃O* DFT 후보 | 3 | 3–4 | 3–4 | 3–4 |
| CO*+CH₃O* reactive pair | 3–5 | 5–8 | 3–5 | 3–5 |
| CO₂ side-path 후보 | optional | 2–3 | 3–5 | 3–5 |

| ID | 작업 | 선행 | 산출물 / 완료조건 |
|----|------|------|-------------------|
| T1.15 | DFT shortlist 선정 (거리 bin별 대표 구조 포함) | T1.14 | 표면·화학종별 shortlist |
| T1.16 | Level 1: PBE-D3 vacuum relaxation | T1.15 | relaxed DFT 구조 + 에너지 |
| T1.17 | Level 2: PBE-D3+VASPsol final energy/relaxation | T1.16 | solvated 에너지 |
| T1.18 | 흡착 에너지 계산 (ΔG_CO*, ΔG_CH3O*: radical 기준 + MeOH(U) 기준 둘 다) | T1.17 | adsorption energy table |
| T1.19 | Descriptor map 작성: ΔG_CO* vs ΔG_CH3O*^MeOH(U) | T1.18 | descriptor map |
| T1.20 | Phase 2 후보 표면 선정 (Case A–D 프레임으로 해석) | T1.19 | 선정 근거 리포트 |

흡착 에너지 정의:
- ΔG_CO* = G_slab+CO − G_slab − μ_CO
- ΔG_CH3O*^rad = G_slab+CH3O − G_slab − G_CH3O
- ΔG_CH3O*^MeOH(U) = G_slab+CH3O + ½G_H2 − G_slab − G_CH3OH − eU

해석 프레임: Case A. Pd(100) 최유리 / Case B. PdO(101)/Pd(100) 최유리 / Case C. bulk·O-rich PdO(100) DMC-inactive / Case D. PdO₂(110) DMC-inactive.

**게이트 G3 (Phase 1 종료):** descriptor map 완성 + Phase 2 후보 표면 확정.

Phase 1 산출물 체크리스트: ① optimized bulk (Pd, PdO, PdO₂) ② clean slabs ③ 표면별 CO* top 2–3 DFT 구조 ④ 표면별 CH₃O* top 2–3 DFT 구조 ⑤ CO*+CH₃O* preliminary 구조 ⑥ adsorption energy table ⑦ descriptor map ⑧ Phase 2 후보 선정.

---

# PHASE 2 — DMC Gibbs free-energy profile

목적: Phase 1을 통과한 표면에서 DMC 형성 경로의 자유에너지 profile과 TS barrier를 계산하고 side-path와 비교.

DMC elementary steps:
1. \* + CO → CO\*
2. CO\* + CH₃OH + \* → CO\* + CH₃O\* + H⁺ + e⁻
3. CO\* + CH₃O\* → CH₃OCO\*
4. CH₃OCO\* + CH₃OH + \* → CH₃OCO\* + CH₃O\* + H⁺ + e⁻
5. CH₃OCO\* + CH₃O\* → DMC\*
6. DMC\* → DMC + \*

| ID | 작업 | 선행 | 산출물 / 완료조건 |
|----|------|------|-------------------|
| T2.1 | reactive CO*+CH₃O* pair sampling (C_CO···O_OCH₃ 2.1–4.0 Å, site-pair ≤4.5 Å) | G3 | reactive endpoint 후보 |
| T2.2 | CH₃OCO*, CH₃OCO*+CH₃O*, DMC* endpoint sampling | T2.1 | intermediate 후보 풀 |
| T2.3 | MLIP relaxation 후 C···O 거리 bin 재분류, DFT shortlist 선정 | T2.1, T2.2 | endpoint shortlist |
| T2.4 | Endpoint PBE-D3+VASPsol 최적화 | T2.3 | relaxed intermediates + 에너지 |
| T2.5 | TS1 (CO*+CH₃O*→CH₃OCO*) NEB/dimer 계산 | T2.4 | TS1 구조 + barrier |
| T2.6 | TS2 (CH₃OCO*+CH₃O*→DMC*) NEB/dimer 계산 | T2.4 | TS2 구조 + barrier |
| T2.7 | side-path 검증: CO + O_lattice → CO₂ + V_O, methanol oxidation | T2.4 | side-path 에너지/barrier |
| T2.8 | DMC Gibbs free-energy profile 작성 + 표면 간 비교 | T2.5–T2.7 | 최종 profile, 결론 리포트 |

NEB 설정: IMAGES=5, SPRING=−5, LCLIMB=.TRUE., IBRION=3, POTIM=0, EDIFFG=−0.05, VASPsol 포함.
비교 기준: 업로드 Angew 논문의 pure Pd TS1/TS2 = 1.08/0.85 eV, Pd₃Cu = 0.86/0.79 eV를 참조 벤치마크로 사용.

**게이트 G4 (프로젝트 종료):** 표면별 DMC profile + side-path 비교로 "anodic 산화가 DMC 촉진 vs 억제"에 대한 결론 도출.

---

## 3. 작업 의존성 요약 (Critical Path)

```
T1.1/T1.2/T1.3 (bulk)
 └─ T1.4 ─[G1]─ T1.5–T1.8 (clean slab)
     └─ T1.9 ─[G2]─ T1.10 (site map)
         └─ T1.11/T1.12 ─ T1.13 ─ T1.14 (MLIP)
             └─ T1.15–T1.20 ─[G3]
                 └─ T2.1 ─ T2.2 ─ T2.3 ─ T2.4
                     ├─ T2.5/T2.6 (TS)
                     ├─ T2.7 (side-path)
                     └─ T2.8 ─[G4]
```

병렬 가능: bulk(T1.1–1.3)는 서로 독립. clean slab(T1.5–1.8)은 해당 bulk 완료 후 표면별 병렬. CO*(T1.11)와 CH₃O*(T1.12)는 병렬. TS1(T2.5)·TS2(T2.6)·side-path(T2.7)는 endpoint 확보 후 병렬.

## 4. 의사결정 게이트 (Go/No-Go)

| 게이트 | 위치 | 통과 조건 | No-Go 시 |
|--------|------|-----------|----------|
| G1 | bulk → slab | 3개 bulk 수렴 + 실험 격자 검증 | 수렴 파라미터 재조정 |
| G2 | slab → sampling | 4개 clean slab 물리량 검증 통과 | termination/슬래브 두께 재검토 |
| G3 | Phase 1 → Phase 2 | descriptor map 완성 + 후보 표면 선정 | DMC-inactive 표면은 Phase 2 제외(side-path만) |
| G4 | 프로젝트 종료 | 표면별 DMC profile + 결론 | TS 재계산 / 추가 중간체 |

## 5. 리스크 및 완화

| 리스크 | 영향 | 완화 |
|--------|------|------|
| 개별 최저 흡착 구조 ≠ 반응 최저 구조 | DMC 경로 오판 | CO*·CH₃O* reactive coadsorption pair 별도 sampling |
| heuristic이 top site만 생성 | bridge/hollow 누락 | CO*/CH₃O*에 mode='all' 강제 |
| O-rich PdO/PdO₂에서 CO₂-like 붕괴 | reactive 구조 손실 | 붕괴 구조를 side-path 후보로 보존·분류 |
| MLIP 정확도 한계 | ranking 오류 | 거리 bin별 대표 구조까지 DFT로 교차검증 |
| asymmetric slab dipole 오차 | 흡착 에너지 편향 | LDIPOL/IDIPOL=3, vacuum 20 Å 유지 |
| PdO₂(110) termination 모호성 | reference 신뢰도 | stoichiometric termination stability 사전 검증(T1.9) |

## 6. 한 줄 요약 (부록 A)

Pd(100)을 metallic nanocube baseline으로, 1 ML PdO(101)/Pd(100)을 부분 산화 Pd/PdO interface로, bulk/O-rich PdO(100)을 thermodynamic Pd²⁺ oxide/O-rich limit로, PdO₂(110)을 Pd⁴⁺ high-valence reference로 사용한다. AutoAdsorbate/heuristic sampling + foundation MLIP로 CO*, CH₃O*, CO*+CH₃O* 후보를 넓게 생성하고, 최종 adsorption energy와 DMC Gibbs profile은 PBE-D3(+VASPsol) DFT로 판정한다.
