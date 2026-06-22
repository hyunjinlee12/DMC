# G1 + G2 Committee Re-audit Summary

**감사 일자**: 2026-06-20
**감사 방법**: Paimon-extended 5-judge blind committee (methods, physics, statistics, silent-error, malicious)
**감사 범위**: G1 bulk (Pd, PdO, PdO₂) + G2 slab (S1, S2, S3, S3b, S4)
**최종 판정**: 두 게이트 모두 **Pass-with-caveats**

---

## G1 Bulk 결과

| Judge | Decision | Score | 핵심 |
|---|---|---:|---|
| methods | Pass | 9 | INCAR/POTCAR/k-mesh/D3 모두 룰 준수 |
| physics | Concern | 8 | PdO₂ lattice +3% (PBE+D3 한계 근처), 잔여 magmom |
| statistics | Pass | 9 | ENCUT 0.17-1.16 meV/atom, k-mesh sweep 완벽 |
| silent-error | Concern | 7 | k-mesh 표 0.2-0.6 meV 차이 (benign), POTCAR 의 valence 클레임은 judge hallucination (memory 에 그런 클레임 없음) |
| malicious | Pass | – | 모든 claim 0.01 meV 정밀도로 reproduce |

### 주요 Lattice 결과 (vs 실험)
- **Pd a = 3.8907 Å (+0.02%)** — EXCELLENT
- PdO a/c: +0.35% / +1.31% (PBE+D3 정상 범위)
- PdO₂ a/c: +1.26% / +2.39% (rutile c축은 PBE+D3 한계)

### Caveats (G1)
1. `G1_report.md` k-mesh table 의 일부 수치가 실제 OUTCAR 와 0.2-0.6 meV 차이 (benign, 수렴 결론은 안전).
2. PdO₂ lattice ±3% 는 PBE+D3 정확도 한계. 흡착 에너지 계산에서 systematic 차이 있을 수 있음 — descriptor map 해석 시 명시.

---

## G2 Slab 결과

| Judge | Decision | Score | 핵심 |
|---|---|---:|---|
| **methods** | **REJECT** ← **FALSE ALARM** | 1 | "4/5 surface 미수렴" — fixed atom force 까지 합산한 오류 |
| physics | GO | 9 | 모든 구조 chemistry 정상 (Pd-O 2.00-2.05 Å, rumpling, dipole) |
| **statistics** | **GO** | **10** | **Free atom 만 분리, 모두 EDIFFG 만족** |
| **silent-error** | **REJECT** ← **FALSE ALARM** | 3 | 같은 fixed atom force 혼동 + termination 라벨 모호성 (true) |
| malicious | Pass | – | STATUS.md E 값에 100 meV 단위 차이 (data extraction 버그) |

### Chair override: **Pass-with-caveats**
methods/silent-error 의 REJECT 는 **VASP convention 오해** (EDIFFG 는 free atom 만 검사) 에서 비롯된 false alarm. 자체 Bash 검증 + statistics judge 결과 일치:

**Free atom max force (실제 VASP 검사 대상):**
- S1: 0.0194 eV/Å ✓
- S2: 0.0261 ✓
- S3: 0.0232 ✓
- S3b: 0.0275 ✓
- S4: 0.0246 ✓

모두 EDIFFG=-0.03 만족, **5 surface 모두 정상 수렴**.

### True findings (받아들임)

#### 1. STATUS.md E 값 데이터 추출 버그
이전: F (free energy, sigma 포함) vs E_sigma→0 혼용 → 100 meV 단위 차이
수정: E_sigma→0 표준으로 통일

| Slab | 이전 | 수정 | 차이 |
|---|---:|---:|---:|
| S1 | -434.408 (F) | -434.380 (E₀) | 28 meV |
| S2 | -618.565 (E₀) | -618.565 | 0 |
| S3 | -724.152 (F) | -724.103 (E₀) | 49 meV |
| S3b | -570.770 | -570.772 | 2 meV |
| S4 | -788.493 | -788.493 | 0 |

#### 2. Termination 라벨 명확화
- S3 "O-term" → 실제 top layer 8 Pd + 16 O (O-rich mixed)
- S3b "PdO-term" → 실제 top layer 순수 Pd (Pd-rich)

→ STATUS.md 의 표는 수정됨. `calculations/G2_slab/{S3_PdO100, S3b_PdO100_PdOterm}` 디렉토리 이름은 backward compatibility 위해 유지.

#### 3. S3 dipole 1.88 e·Å (physics)
다른 slab (0.02-0.26) 대비 큼. IDIPOL=3 보정됨. 흡착 에너지 계산 시 systematic offset 모니터링 권고.

---

## 결론

**두 게이트 모두 통과 확정.**

- G1 통과 (3 bulk converged, lattice 검증)
- G2 통과 (5 slab converged, T1.9 validation)
- Downstream G3 결과 (Phase 1/2/3 ranking + DFT shortlist) **재계산 불필요** — 기존 결과 그대로 사용 가능.
- STATUS.md 갱신 완료 (E 값 표준화 + termination 라벨 명확화).

### Lessons learned (memory 업데이트)
- VASP EDIFFG 는 **free atom (Selective Dynamics T T T) 만 검사** — fixed atom 무시.
- Judge 에게 DFT relax 검증 dispatch 시 위 logic 을 brief 에 명시.
- STATUS.md / 보고서 의 E 값은 항상 E_sigma→0 (E₀) 사용.
- `memory/feedback_vasp_convergence_audit.md` 에 기록.

---

## 다음 단계

T1.15 DFT shortlist (47 jobs) → **DFT 진입 user 승인 게이트** (memory rule: 비싼 작업).
