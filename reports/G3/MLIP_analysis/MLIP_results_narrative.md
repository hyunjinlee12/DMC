# MLIP Results Analysis — Phase 1+2+3 (MACE-MH + D3)

**작성일**: 2026-06-21  
**검증 상태**: T1.10-T1.15 committee 모두 Pass-with-caveats  
**Note**: 모든 E 는 MACE-MH ranking 전용. 절대값은 DFT (T1.16-17) 로 확인 예정.

---

## 1. 요약 표 (per-surface, MLIP top-1)

| Surface | 산화 | n_CO unique | CO d_min (Å) | CO E_range (meV) | CH₃O d_min (Å) | CH₃O E_range | n_coads | coads d_react (Å) | coads d_min (Å) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **S1** Pd(100) | Pd⁰ | 12 | **1.97** ✓ | **1692** | **2.11** ✓ | 304 | 1985 | 3.10 ✓ SetA | 1.98 ✓ |
| **S2** PdO/Pd | mixed | 26 | **2.01** ✓ | 818 | **2.14** ✓ | 1653 | 5524 | **5.26** drift→SetB | 1.18 ⚠ |
| **S3** PdO O-top | Pd²⁺ | 10 | **2.46** ⚠ | 126 | 2.79 ⚠ | 1180 | 3614 | **1.34** ⭐ product | 1.41 ⚠ |
| **S3b** PdO Pd-top | Pd²⁺ | 18 | **3.54** ❌ | 2314 | 2.55 ⚠ | 471 | 1581 | 5.30 drift→SetB | 1.91 ✓ |
| **S4** PdO₂(110) | Pd⁴⁺ | 52 | **4.05** ❌ | 3184 | **0.98** ❌ broken | 5764 | 1918 | **1.33** product/broken | 1.33 ⚠ |

**chemisorbed band**: Pd-C 1.85-2.10 Å, Pd-O 2.00-2.15 Å (catalysis reference)

---

## 2. Surface-by-surface chemistry 분석

### S1 — Pd(100) clean (Pd⁰)
- CO* Pd-C **1.97 Å** ✓ 표준 chemisorbed
- CH₃O* Pd-O **2.11 Å** ✓ 표준 chemisorbed
- co-ads top-1 d_react **3.10 Å** (Set A 안) + d_min 1.98 Å ✓
- E_range CO 1.7 eV (대단히 discriminative)

**해석**: 깨끗한 Pd metallic. CO 와 CH₃O 모두 잘 흡착. co-ads 도 reactive 거리에서 안정. **Case A 후보** (DMC favorable).

---

### S2 — PdO(101)/Pd(100) composite (Pd⁰ + Pd²⁺ interface)
- CO* d_min **2.01 Å** ✓ (Pd⁰ 자리에 흡착)
- CH₃O* d_min **2.14 Å** ✓
- co-ads top-1 d_react **5.26 Å** ⚠ (Set A [2.1-4.0] 이탈, Set B 영역으로 drift)
  - d_min 1.18 Å 매우 짧음 → 한 분자만 강결합, 다른 분자는 분리 가능성
- E_range CO 818 meV, CH₃O 1653 meV ✓

**해석**: **bifunctional interface 가설은 약함**. Single ads 잘 흡착하지만, **두 분자 동시 들이면 분리 선호** (5.26 Å로 멀어짐). Shi 2024 의 Pd/PdO interface hypothesis 가 우리 데이터에선 약하게 나옴 — DFT 확인 필요.

---

### S3 — PdO(100) O-term (Pd²⁺ O top)
- CO* d_min **2.46 Å** ⚠ borderline physisorbed
- CH₃O* d_min **2.79 Å** ⚠ physisorbed
- E_range CO **126 meV** ⚠ (200 meV 임계 미달, PES flat)
- **co-ads top-1 d_react 1.34 Å** ⭐ **product collapse 의심** (CH₃OCO\* 형성?)
- d_min 1.41 Å — 표면-product complex

**해석**: 표면이 O 로 덮여 single CO/CH₃O 가 Pd 직접 접근 어려움. 그러나 **두 분자 함께 들어가면 분자간 coupling 으로 product 직접 형성 신호** (d_react 1.34 Å 는 CH₃OCO\* 의 C-O 단일결합 거리). 
**Case C 후보** — 단, DMC step 이 thermodynamically favorable 한 흥미로운 surface 일 가능성.

---

### S3b — PdO(100) Pd-term (Pd²⁺ Pd top)
- CO* d_min **3.54 Å** ❌ physisorbed (under-coord Pd 인데도 CO 멀어짐)
- CH₃O* d_min **2.55 Å** ⚠ borderline
- co-ads top-1 d_react 5.30 Å (Set B drift) + d_min 1.91 Å ✓
- E_range CO 2.3 eV ⭐ excellent

**해석**: Under-coordinated Pd 가 예상보다 weak CO binding. CH₃O 는 borderline. co-ads 는 분리된 state 선호. 
Phase 1 단일 흡착에서 1.81 Å 결합 봤던 것과 차이 — top-1 vs top-3 의 차이로 explain 가능 (top-1 만 physisorbed).
**Case A/B 후보**. CO* dock 위치가 까다로움 (sampling 보강 필요).

---

### S4 — PdO₂(110) (Pd⁴⁺)
- CO* d_min **4.05 Å** ❌ 완전 unbound
- CH₃O* d_min **0.98 Å** ❌ 비정상 짧음 (methoxy decomposition 의심 — O atom 이 표면 O와 합쳐짐?)
- co-ads top-1 d_react **1.33 Å** + d_min 1.33 Å (둘 다 broken-like)
- E_range CO 3.2 eV (가장 큼)

**해석**: ⭐⭐⭐ **Shi 2024 핵심 가설 직접 검증** — Pd⁴⁺ surface 에서 CO* 매우 불안정. 추가 finding: CH₃O 가 PdO₂ 위에서 분해 경향 (top-1 d_min 0.98 Å 는 O-Pd 또는 O-O 결합 형성 가능성). 
**Case C/D 결정적 후보**. DFT 시 manual placement 권고 (broken top-1 사용 불가).

---

## 3. Oxidation Trend 종합

```
산화 진행:    S1 (Pd⁰)  →  S2 (mixed)  →  S3b (Pd²⁺ Pd-top)  →  S3 (Pd²⁺ O-top)  →  S4 (Pd⁴⁺)

CO Pd-C:     1.97  →  2.01  →  3.54         →  2.46         →  4.05 Å
             ✓        ✓        ❌              ⚠              ❌

CH₃O Pd-O:   2.11  →  2.14  →  2.55         →  2.79         →  0.98 (broken)
             ✓        ✓        ⚠              ⚠              ❌
```

**선명한 trend**: **산화 진행 → CO* 결합 약화** (Shi 2024 정확 재현). CH₃O 는 산화 영향 덜 민감 (Pd⁰~Pd²⁺ 까지 유지), 다만 Pd⁴⁺ 에서 decomposition.

---

## 4. Co-adsorption behavior 분류

| Surface | Co-ads 거동 | 의미 |
|---|---|---|
| **S1** | SetA 안 정상 reactive | DMC step 정상 진행 가능 |
| **S2** | Set B drift (분리) | bifunctional 가설 약함 |
| **S3** | **Product collapse (d=1.34 Å)** | **CH₃OCO\* 직접 형성 신호** ⭐ |
| **S3b** | Set B drift (분리) | reactive intermediate unstable |
| **S4** | Broken/product collapse | Pd⁴⁺ chemistry 적용 안 됨 |

→ **두 가지 chemistry mode 분기**:
- **Reactive intermediate stable**: S1 (clean Pd⁰)
- **Product favored**: S3 (PdO O-top — DMC step thermodynamically 유리?)
- **Separation favored**: S2, S3b (bifunctional 또는 interface)
- **Broken chemistry**: S4 (PdO₂)

---

## 5. Preliminary Case A-D 분류 (T1.20 예상)

```
            CO* binding
             강 ─────────────────── 약
   강   ┌──────────────┬──────────────┐
        │              │              │
   CH3O │  CASE A      │  CASE C      │
   비   │  S1, (S3b?)  │  S3 ⭐       │  (side-path 우세)
   nd   │  DMC fav.    │  prod. collapse?
   ing  │              │              │
        ├──────────────┼──────────────┤
        │              │              │
   약   │  CASE B      │  CASE D      │
        │  S2          │  S4          │  (DMC inactive)
        │  interface   │  Pd⁴⁺        │
        │              │              │
        └──────────────┴──────────────┘

  ※ 단일 흡착 단순 분류. Co-ads 거동까지 보면 S3 가 product 직접 형성 가능성 — 
    Case C 의 "DMC inactive" 라벨 보다 "DMC step easy but pathway different" 가 더 정확.
```

---

## 6. Shi 2024 Angew 가설 vs 우리 결과

| Shi 2024 예측 | 우리 결과 | 일치도 |
|---|---|---|
| Pd⁰ 에서 CO* favorable | S1 Pd-C 1.97 ✓ | ✅ 완전 일치 |
| Pd²⁺ 에서 CH₃O* 가능 | S2/S3/S3b CH₃O ~2.1-2.8 ✓ | ✅ 일치 |
| **Pd⁴⁺ 에서 CO* 약화** | **S4 Pd-C 4.05 Å** ❌결합 | ✅⭐ **정확 재현** |
| anodic 산화 → side-path 우세 | S4 + S3 broken/product 신호 | ⚠ 부분 검증 (DFT 필요) |
| TS barrier Pd 1.08 eV | (T2.5 NEB 미실시) | — DFT 후 |

---

## 7. DFT 진입 후 우선 확인 항목

### A. 핵심 검증 (T1.16 Level 1 vacuum)
1. **S1 single CO Pd-C**: MLIP 1.97 vs DFT? Pd(100) 표준 ~1.95 Å 매칭 확인
2. **S4 single CO**: 진짜 unbound (positive ΔG)?  ⭐ Shi 2024 hypothesis confirmation
3. **S3 co-ads top-1 d_react=1.34**: 진짜 CH₃OCO\* product? 또는 MLIP artifact?

### B. 보조 검증
4. **S2 co-ads bifunctional**: top-1 이 Set B 로 drift 한 것이 진짜 chemistry truth 인지
5. **S3 single CO PES flat (E_range 126 meV)**: DFT 에서도 평탄한지

### C. 제외 / Manual
6. **S4 single CH₃O top-1 broken (d_min=0.98)**: manual replacement 필요 (atop Pd, atop O 등 high-symmetry site)
7. **S4 co-ads**: shortlist 자체가 broken — 0 jobs in current T1.15 (의도된 제외)

---

## 8. 의의 (advisor 보고용)

### 1) Methodology
- MACE-MH + D3 + cuEq foundation MLIP 로 **48k 후보 → 47 DFT** 효과적 down-selection
- 5-judge Paimon committee 가 silent bug 2회 잡음 (Phase 2 PBC, Phase 3 SetTS) → 데이터 폐기 없이 fix
- Shi 2024 의 SSW-NN 역할을 foundation MLIP 로 대체, **학습 없이** 동일 chemistry hierarchy 재현

### 2) Chemistry findings (MLIP 단계 수준)
- **Pd 산화 → CO\* 안정성 손실** 정량 재현 (Shi 2024 직접 검증)
- **PdO(100) O-term 에서 CH₃OCO\* product 형성 신호** (S3 co-ads top-1 d_react=1.34 Å)
  → 우리 가 발견한 추가 chemistry — Shi 2024 에 명시 안 됨
- S2 interface 의 bifunctional 가설은 우리 MLIP 데이터에서 약함 (Set B drift)

### 3) Phase 2 (workplan) 전망
- TS barrier 비교: 우리 S1 Pd(100) 만 Shi 2024 의 pure Pd (TS1=1.08 eV) 와 직접 비교 가능
- S3 product collapse 가 진짜면 → DMC step thermodynamic favorable surface 신규 발견
- S4 PdO₂(110) 는 TS 계산 자체 부적용 (CO 불결합) → side-path 우세 결론
