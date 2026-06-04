# 참고 문헌 요약 — Pd/PdO/PdO₂ 표면 DMC Formation DFT 연구

작성일: 2026-06-03
근거: 대화록(pd_dmc_conversation_export.md) + 웹 검증(2026-06 기준)
표기: 웹에서 직접 확인 / 프로젝트 노트(대화록)에 기술된 내용 — 인용 전 원문 재확인 권장.

## A. 앵커 논문 — Pd 전극 DMC 전기합성

### A1. Shi et al., Angew. Chem. Int. Ed. 2024 — Pd₃Cu 전기촉매

- 제목: Stabilization of Pd⁰ by Cu Alloying: Theory-Guided Design of Pd₃Cu Electrocatalyst for Anodic Methanol Carbonylation
- DOI: 10.1002/anie.202401311 · PubMed 38606491
- 문제의식: CO + methanol → DMC anodic carbonylation은 높은 작동 전위와 낮은 DMC 선택성(심한 methanol self-oxidation 동반)이 한계.
- 핵심 주장: Cu를 Pd에 도핑하면 산화 환경에서 저원자가 Pd(Pd⁰)를 안정화하고, 동시에 DMC 형성의 전체 에너지 장벽을 낮춤.

DFT 결과(대화록 정리):
- Pd⁰에서는 CO* adsorption이 유리하지만 Pd⁴⁺/PdO₂-like 표면에서는 CO* adsorption free energy가 양(+)으로 이동.
- CH₃O* adsorption은 Pd⁰·Pd⁴⁺ 모두 가능 → Pd 산화 시 DMC carbonylation보다 methanol oxidation side reaction 우세 가능.
- TS barrier: pure Pd = TS1/TS2 1.08 / 0.85 eV, Pd₃Cu = 0.86 / 0.79 eV로 감소.

실험: 2-step thermal shock 합성 Pd₃Cu가 DMC onset을 0.7 V vs Ag/AgCl까지 낮추고, 1.0 V에서 DMC FE 93% 달성.

본 프로젝트에서의 역할: TS barrier 벤치마크 + "anodic 산화가 CO* 안정성을 무너뜨려 side-path를 키운다"는 핵심 가설의 출발점.

### A2. 관련 DMC 전기합성 문헌 (보조)

- Guo et al., Angew. Chem. — CO₂로부터 DMC tandem 전기합성. CO₂ 경로 비교용.
- Mi et al., Small 2025 — intermediate stabilization·coupling 촉진으로 DMC 전기합성 (DOI 10.1002/smll.202501780). intermediate binding 전략 비교용.

## B. Pd 표면 산화 및 PdO / PdO₂ 표면 구조

### B1. Pd(100) surface oxide = strained PdO(101)/Pd(100)

- 핵심 문헌: "The Pd(100)-(√5×√5)R27°-O surface oxide revisited" (arXiv cond-mat/0304107) 및 후속 LEED+DFT+STM 정밀화.
- 결론: 종래의 rumpled PdO(001) 모델은 HRCLS/STM/DFT와 불일치. 실제로는 Pd(100) 위에 올라간 strained PdO(101) 단층이 맞으며, tensor-LEED I(V)와 DFT 좌표가 잘 일치(이후 PdO(101) 층의 수평 shift만 보정).
- 본 프로젝트 함의: 표면 S2 "1 ML PdO(101)/Pd(100)" 모델의 구조적 정당성 — Pd(100) 산화 초기 surface oxide의 표준 모델.

### B2. CO oxidation: Pd(100) vs PdO(101)-(√5×√5)R27°

- First-principles kinetic phase diagram + bistability 연구.
- 함의: working catalyst는 metal/oxide 사이를 전위·조건에 따라 동적으로 전환(ab initio thermodynamics + kMC). clean Pd만 계산하면 부족 → S1→S2→S3→S4 hierarchy 필요성의 근거.

### B3. PdO / Pd/PdO redox dynamics 관련 (대화록 정리)

| 문헌/방향 | 의미 |
|-----------|------|
| Epitaxial PdO(100) during Pd(100) oxidation | PdO(101)/Pd(100)만이 아니라 PdO(100)/Pd(100)도 후보로 검토 필요 |
| operando TEM + NAP-XPS + DFT (Pd/PdO redox) | strained PdO와 Pd/PdO interface가 반응성에 중요 |
| CO titration / DRIFTS / XRD / STEM + DFT | 실제 working catalyst는 topmost PdOₓ layer를 형성할 수 있음 |
| CO chemisorption DRIFTS on PdO | fully oxidized PdO는 CO adsorption에 약할 수 있음 (→ Case C/D의 DMC-inactive 가설 근거) |

위 항목들은 대화록 표(§4)에 정리된 방향성으로, 실제 논문 인용 시 서지정보를 재확인할 것.

### B4. PdO₂ / Pd⁴⁺ 표면 모델 (대화록 정리)

PdO₂는 dominant facet로 단정하지 말고 세 목적으로 분리 모델링:
1. DMC 문헌 benchmark — Angew 논문의 PdO₂-like reconstructed Pd(211) p(2×1)
2. bulk Pd⁴⁺ oxide reference — rutile/hydrophilite-like PdO₂(110) (← 최종 채택 S4)
3. Pd nanocube 현실성 — Pd(100)-derived PdO₂-like/PdO₂₋ₓ overlayer

SSW-NN + VASP DFT 등으로 high-valence Pd가 CO* adsorption을 약화시킬 수 있음이 보고됨.

## C. 표면 DMC formation 메커니즘 및 연구 흐름 (대화록 정리)

전체 반응: 2 CH₃OH + CO − 2e⁻ → (CH₃O)₂CO + 2H⁺
Elementary steps: CO+*→CO* → CH₃OH→CH₃O*+H⁺+e⁻ → CO*+CH₃O*→CH₃OCO* → CH₃OCO*+CH₃O*→DMC+*

연구 흐름 타임라인:

| 시기 | 초점 | 핵심 메시지 |
|------|------|-------------|
| ~2017 | in situ FTIR + 금속 전극 비교 | Au, Pd, Pt, Ag 등에서 DMC 형성 관찰 |
| 2019 | 금속 표면 DFT descriptor screening | *CO, *OCH₃ 흡착 자유에너지 + product reaction energy로 가능성 평가 |
| 2019 | Pd-B dopant | intermediate binding 조절로 C–O coupling barrier 저감 |
| 2024 | anodic Pd oxidation/reconstruction | Pd 산화 시 CO* adsorption 불리 → DMC보다 DMM/MF/CO₂ 우세 가능 |
| 2025–2026 | Pd-Au, Pd/PdO, Pd-Br 등 | Pd⁰ 유지, Pd²⁺/halide/interface 활용, 고전류·장시간 운전 방향 |

핵심 해석 요약:
- Pd⁰ surface: CO* stabilization 유리
- PdO/Pd²⁺ surface: CH₃O* stabilization 유리 가능
- PdO₂/Pd⁴⁺ surface: CO* adsorption 약화 가능
- Pd/PdO interface: CO*와 CH₃O*를 동시에 안정화하는 bifunctional site 가능성

## D. 계산 방법론 참고

### D1. AutoAdsorbate — heuristic 흡착 구조 생성

- 관련 워크플로우 논문: ACS Catal. 2025, DOI 10.1021/acscatal.5c06553 (foundation MLIP 기반 흡착 사이트 sampling/안정 구조 예측). 서지/저자 인용 시 재확인.
- 핵심 object: `from autoadsorbate import Surface, Fragment`
  - `Surface(slab)` — 노출 surface atom·가능한 adsorption site 탐색
  - `Fragment('*SMILES')` — 흡착 fragment 생성 (CO* = `Cl[C-]#[O+]`, CH₃O* = `ClOC`)
  - `get_populated_sites(...)` — fragment를 site에 배치
- 주의: `mode='heuristic'`는 connectivity==1(top site)만 배치 → CO*/CH₃O*는 top/bridge/hollow를 모두 봐야 하므로 `mode='all'` 필수. co-adsorption은 `get_populated_sites()`만으로 부족 → site map 선생성 후 custom wrapper로 site-pair 순차 배치.

### D2. Foundation MLIP — pre-screening/ranking

MACE-MP-0 (MPtrj 학습 foundation model)은 ASE calculator로 즉시 사용 가능(`from mace.calculators import mace_mp`). out-of-distribution 일반화가 좋고 heterogeneous catalysis(CO oxidation, CO₂→methanol 등)에서 TS 위치를 대략 맞춤 — 단, 절대 에너지는 부정확하므로 ranking·pre-screening 용도로만 사용하고 최종 판정은 DFT로.

정량화가 필요하면: 각 에너지 프로파일에서 소수(예: 5개) DFT single-point로 fine-tuning하면 DFT reference와 정량적으로 일치(MACE-MP-0b3+D3가 O/OH adsorption scaling 재현). MACE-OMol25 등 최신 foundation potential도 후보.

본 프로젝트 적용: 넓게 생성한 CO*/CH₃O*/co-adsorption 후보를 MLIP로 relax·ranking → 표면·화학종별 2–3개만 DFT.

### D3. DFT 레벨 — VASP PBE+D3 (+VASPsol)

- Functional: PBE+D3 (GGA=PE, IVDW=12 = D3-BJ), ENCUT=520, PREC=Accurate, LASPH=.TRUE., ADDGRID=.TRUE., ISPIN=2.
- Asymmetric slab: top side만 adsorption, bottom 2 layers(또는 30–40%) fixed, vacuum 20 Å, dipole correction LDIPOL=.TRUE., IDIPOL=3.
- Smearing: Pd metal ISMEAR=1/SIGMA=0.10, PdO·PdO₂ ISMEAR=0/SIGMA=0.05.
- Implicit solvation: VASPsol LSOL=.TRUE., EB_K=32.6(물), TAU=0.
- TS: CI-NEB IMAGES=5, SPRING=-5, LCLIMB=.TRUE., IBRION=3, POTIM=0 (또는 dimer).

## E. 한 줄 종합

앵커는 Shi 2024 Angew(Pd₃Cu) — "Pd 산화가 CO* 안정성을 무너뜨려 DMC 대신 side-path를 키운다"는 가설과 TS barrier 벤치마크를 제공한다. √5×√5 surface oxide = strained PdO(101)/Pd(100) 문헌이 S2 모델의 구조적 근거이고, Pd/PdO redox·kinetic phase diagram 문헌이 S1→S4 산화 hierarchy의 정당성을 준다. 방법론은 AutoAdsorbate(구조 생성) → foundation MLIP(MACE, ranking) → VASP PBE-D3(+VASPsol)(판정)의 3단 파이프라인.

## 출처 (웹 확인)

- Shi et al., Angew. Chem. Int. Ed. 2024 — https://onlinelibrary.wiley.com/doi/10.1002/anie.202401311 · https://pubmed.ncbi.nlm.nih.gov/38606491/
- Pd(100)-(√5×√5)R27°-O surface oxide revisited — https://arxiv.org/abs/cond-mat/0304107
- MACE foundation models — https://mace-docs.readthedocs.io/en/latest/guide/foundation_models.html · https://pubs.aip.org/aip/jcp/article/163/18/184110/3372267/
- MACE fine-tuning tutorial — https://arxiv.org/html/2506.21935v2
- Mi et al., Small 2025 (DMC) — https://onlinelibrary.wiley.com/doi/10.1002/smll.202501780
- (재확인 필요) AutoAdsorbate / ACS Catal. 10.1021/acscatal.5c06553 — https://pubs.acs.org/doi/10.1021/acscatal.5c06553
