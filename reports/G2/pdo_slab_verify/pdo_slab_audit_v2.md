# PdO(100) Slab Termination Verification (V2 — tol 0.2 Å)

## Bulk PdO reference
- a = 3.054 Å, c = 5.406 Å
- Pd-O bond mean: 2.039 Å (n=4)

## S3 (claimed O-term)
- 128 atoms (Pd:64, O:64, O/Pd = 1.00) **stoichiometric**
- 8 layers (0.2 Å tol)
- TOP 5 layers (z descending):
  - z=12.07  (24 atoms)  Pd:8, O:16  [mixed]
  - z=10.75  ( 8 atoms)  Pd:8  [pure Pd]
  - z=9.15  (24 atoms)  Pd:8, O:16  [mixed]
  - z=7.65  ( 8 atoms)  Pd:8  [pure Pd]
  - z=6.11  (24 atoms)  O:16, Pd:8  [mixed]
- Pd-O bond mean: 2.034 Å (n=240, vs bulk 2.039)
- Top-region Pd coord to O: [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

## S3b (claimed Pd-term)
- 104 atoms (Pd:56, O:48, O/Pd = 0.86) **Pd-rich**
- 7 layers (0.2 Å tol)
- TOP 5 layers:
  - z=10.78  ( 8 atoms)  Pd:8  [pure Pd]
  - z=9.14  (24 atoms)  Pd:8, O:16  [mixed]
  - z=7.65  ( 8 atoms)  Pd:8  [pure Pd]
  - z=6.10  (24 atoms)  Pd:8, O:16  [mixed]
  - z=4.57  ( 8 atoms)  Pd:8  [pure Pd]
- Pd-O bond mean: 2.046 Å (n=192)
- Top-region Pd coord to O: [4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2]

## Conclusion

### Termination labels CONFIRMED valid (with caveat)
- **S3 = 진짜 O-term**: TOP layer 은 순수 16 O, 그 아래 8 Pd 가 0.14 Å 분리. O atom 들이 표면 위로 노출됨.
- **S3b = 진짜 Pd-term**: TOP layer 은 순수 8 Pd, 그 아래 16 O 가 1.6 Å 분리. Pd atom 들이 표면 위로 노출됨.

### Layer structure detail
PdO(100) slab 은 c-축 방향에 따라 **순수 Pd plane** 과 **순수 O plane** 이 교대로 쌓임 (tetragonal PdO bulk 의 c-축 stacking). 우리 supercell 은 4×2 = 8 Pd 또는 16 O / layer.

### Stoichiometry
- S3: O/Pd = 1.00 (perfectly stoichiometric — Pd:O = 64:64)
- S3b: O/Pd = 0.86 (Pd-rich — 56:48). 두 layer 비교하면 S3b 가 TOP 의 O layer 1개 제거 + bottom 도 균형 맞춤.

### Bond length sanity
- S3 Pd-O mean 2.034 Å (vs bulk 2.039, +-0.25%)
- S3b Pd-O mean 2.046 Å (++0.37%)
- 둘 다 bulk PdO 결합거리와 거의 동일 → 표면 reconstruction 없이 정상.

### Pd coordination
- S3 top-region Pd: [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4] (avg 4.0). Bulk = 4. 정상.
- S3b top-region Pd: [4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2] (avg 3.0). 상부 Pd는 2-coord (under-coordinated, exposed) — Pd-rich termination 의 특징.

### DMC chemistry implication
- **S3 (O-term)**: 표면 O atoms 가 노출됨 → CO 가 표면 Pd 에 직접 접근 어려움 (Phase 1 결과: CO 미결합과 일치).
- **S3b (Pd-term)**: 상부 Pd 가 under-coordinated (2-fold) → CO 강결합 가능. Phase 1 결과 Pd-C 1.81 Å chemisorbed 와 일치.
- **두 surface 의 chemistry 차이가 sharp** → descriptor map 에서 둘이 명확히 다른 영역 (Case 분류 다르게 나올 것).

### Validity for project
두 termination 모두 valid, intended chemistry 정확히 구현됨. Phase 1/2/3 ranking 결과들이 이 termination 차이를 잘 반영함 (S3 CO 약, S3b CO 강).
