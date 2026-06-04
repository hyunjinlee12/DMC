# 학생용 가이드 — Pd/PdO/PdO₂ DMC Formation DFT 프로젝트

작성일: 2026-06-03 · 지도: 이태훈
먼저 읽을 것: 이 가이드 → 워크플랜 → 패키지 가이드라인 → 참고문헌 요약

이 프로젝트는 Pd nanocube {100} 표면이 anodic 조건에서 산화될 때 CO*·CH₃O* 흡착 균형이 어떻게 변하고, 그 변화가 DMC 형성을 촉진하는지/억제하는지를 DFT로 검증한다. 작업 전체는 두 단계다 — Phase 1(표면 검증 + 흡착 맵), Phase 2(DMC 자유에너지 profile).

## 1. 3개 문서에 무엇이 있나

### ① DMC_Pd_workplan.md — 워크플랜 (무엇을, 어떤 순서로)
- 최종 4개 표면(S1 Pd(100) → S2 1 ML PdO(101)/Pd(100) → S3 bulk/O-rich PdO(100) → S4 PdO₂(110))과 그 역할.
- Phase 1/2를 작업 ID(T1.1 … T2.8)로 분해. 각 작업의 선행조건·산출물·완료조건.
- 의사결정 게이트 G1–G4, critical path 의존성, 리스크·완화표.
- 언제 보나: "지금 뭘 해야 하고 다음은 뭔지" 확인할 때. 작업의 뼈대.

### ② DMC_Pd_package_guideline.md — 패키지 가이드라인 (어떻게)
- 각 작업 ID에서 어떤 패키지를 어떤 옵션으로 쓰는지. ASE(허브) / pymatgen(slab) / AutoAdsorbate(흡착 생성) / MACE(MLIP ranking) / VASP+VASPsol(판정) / VTST NEB(TS).
- 실행 코드 스니펫(CO*/CH₃O* sampling, co-adsorption wrapper, MLIP relax), INCAR 옵션, cutoff, 실무 주의사항.
- 언제 보나: 실제로 계산을 돌리기 직전. 명령·옵션 레퍼런스.

### ③ DMC_Pd_references_summary.md — 참고문헌 요약 (왜)
- 앵커 논문(Shi 2024 Angew, Pd₃Cu)과 TS barrier 벤치마크, Pd 표면 산화·PdO/PdO₂ 구조 문헌, DMC 메커니즘 연구 흐름.
- 각 표면 모델을 왜 그렇게 잡았는지의 근거(예: S2 = strained PdO(101)/Pd(100) 실측 문헌).
- 언제 보나: 모델 선택의 이유가 헷갈릴 때, 결과를 문헌과 비교할 때, 발표·논문 작성 시.

한 줄: 워크플랜=지도, 패키지 가이드=공구 설명서, 참고문헌=배경지식.

## 2. 나(지도교수)와 논의해야 하는 체크포인트

각 단계 끝난 직후, 다음 단계로 넘어가기 전에 아래 자료를 들고 오기. 게이트(G)를 통과해야 다음으로 진행.

### 체크포인트 A — DFT 구조 최적화 후 (bulk, T1.1–1.4 / 게이트 G1)
- 가져올 것: Pd·PdO·PdO₂ 최적 격자상수, 실험값 대비 오차, ENCUT·k-mesh 수렴 그래프.
- 함께 볼 것: bulk가 신뢰할 만한가? smearing(Pd는 ISMEAR=1, 산화물은 0) 적절했나? → 통과해야 slab 제작 착수.

### 체크포인트 B — Surface(clean slab) 구조 제작 후 (T1.5–1.9 / 게이트 G2)
- 가져올 것: 4개 표면 clean slab 그림, surface rumpling·Pd–O bond length·coordination·slab dipole, O-rich/PdO₂ termination 안정성.
- 함께 볼 것: S2 PdO(101)/Pd(100) 정합·strain이 합리적인가? termination 선택이 타당한가? 층수·진공(20 Å)·하단 고정 일관성. → 통과해야 흡착 sampling 착수.

### 체크포인트 C — Adsorbate initial population 생성 후 (T1.10–1.13)
- 가져올 것: CO*·CH₃O* 후보 개수와 site 분포, co-adsorption Set A(reactive)/Set B(thermo)/side-path 분류 통계, 적용한 cutoff.
- 함께 볼 것: mode='all'로 top/bridge/hollow 다 봤나? co-adsorption이 reactive pair를 충분히 포함하나? S2는 interface pair가 들어왔나? → 분포가 빈약하면 다시 생성.

### 체크포인트 D — MLIP 계산 + 2–3 구조 screening 후 (T1.14–1.15)
- 가져올 것: MLIP 에너지 순위표, 표면·화학종별로 고른 2–3개 DFT 후보와 선정 이유(거리 bin별 대표 포함).
- 함께 볼 것: ranking이 물리적으로 말이 되나? 최저에너지만 고르지 않고 reactive 구조를 챙겼나? (MLIP는 ranking 전용 — 절대값 신뢰 금지) → 합의된 shortlist만 DFT로.

### 체크포인트 E — 최종 DFT 계산 후 (T1.16–1.20 / 게이트 G3, 그리고 Phase 2 T2.x / G4)
- Phase 1 (G3): adsorption energy table, descriptor map(ΔG_CO* vs ΔG_CH3O*^MeOH(U)), Case A–D 해석 → Phase 2로 넘길 표면 선정.
- Phase 2 (G4): endpoint 최적화·TS1/TS2 barrier·side-path 비교·DMC Gibbs profile → 최종 결론. Shi 2024 벤치마크(pure Pd 1.08/0.85 eV, Pd₃Cu 0.86/0.79 eV)와 비교.

## 3. 체크포인트 ↔ 워크플랜 대응

| 체크포인트 | 작업 ID | 게이트 | 핵심 산출물 |
|-----------|---------|--------|-------------|
| A. DFT 구조 최적화 후 | T1.1–1.4 | G1 | 최적 bulk 3종 + 수렴 리포트 |
| B. Surface 제작 후 | T1.5–1.9 | G2 | clean slab 4종 + 물리량 검증 |
| C. Adsorbate population 후 | T1.10–1.13 | — | 흡착 후보 풀 + 분류 통계 |
| D. MLIP + screening 후 | T1.14–1.15 | — | 순위표 + DFT shortlist(2–3) |
| E. 최종 DFT 후 | T1.16–1.20 / T2.x | G3 / G4 | descriptor map / DMC profile |

## 4. 매 논의 때 지켜줄 것

- 막혔으면 혼자 오래 붙잡지 말고 해당 체크포인트에서 바로 공유 (특히 termination 선택, co-adsorption cutoff, MLIP ranking 판단).
- 게이트를 건너뛰고 다음 단계로 가지 말 것 — G2 전에 sampling, G3 전에 Phase 2 시작 금지.
- 결과는 항상 구조 그림 + 수치 표로. 구조는 ASE GUI/VESTA로 육안 확인 후 공유.
- 사용한 INCAR/KPOINTS/패키지 버전(특히 MACE 체크포인트)을 결과와 함께 기록 — 재현성.
