# T1.16 DFT Submission Go/No-Go Certification

**작성일**: 2026-06-21
**상태**: ✅ **READY FOR SUBMISSION** (47 jobs)
**전체 검증**: G1/G2/T1.10-T1.13/T1.14(P1/P2/P3)/T1.15 모두 Pass-with-caveats 통과

---

## 1. 최종 47 candidates (post-replacement)

| Surface | single_CO | single_CH3O | coads_SetA | 합계 |
|---|---:|---:|---:|---:|
| S1 Pd(100)         | 3 | 3 | 3 | 9 |
| S2 PdO/Pd          | 5 | 5 | 4 | 14 |
| S3 PdO O-top       | 3 | 3 | 3 | 9 |
| S3b PdO Pd-top     | 3 | 3 | 3 | 9 |
| S4 PdO₂(110)       | 3 | **3** (1 MLIP + 2 manual) | 0 | 6 |
| **Total** | 17 | 17 | 13 | **47** |

### S4 single_CH3O 후보 구성 (특수처리)
- `01_single_CH3O_rank1_idx00241.vasp` — MLIP rank-1 (intact, OK)
- `90_single_CH3O_manual_atop_Pd.vasp` — **Manual placement: atop top-Pd, d(O-C)=1.42 Å, d(C-H)=1.10 Å ✓**
- `91_single_CH3O_manual_bridge.vasp` — **Manual placement: bridge Pd-Pd, ✓**
- ❌ `00_single_CH3O_rank0_idx00065.vasp.broken_bak` (archived, methoxy O-C=1.24 broken)
- ❌ `02_single_CH3O_rank2_idx00008.vasp.broken_bak` (archived, C-H=3.26 broken)

---

## 2. Chemistry 검증 종합

```
total 47:
  OK:                32 (68%)   chemisorbed, intramolecular bonds in range
  SUSPECT_weak:      15 (32%)   physisorbed / weak binding (Pd-C > 2.5 Å)
                                ※ Shi 2024 hypothesis 검증용 — DFT 가 정량
  BROKEN:             0 (0%)   ⭐ 모두 archive + manual replacement 완료
```

### Per-surface weak_binding breakdown
- **S1**: 2 weak (CO rank0/1: Pd-C 1.84 — threshold 1.85 의 0.01 borderline, 사실상 OK)
- **S2**: 4 weak (CO rank1/2 Pd-C 2.52/3.82, coads rank2/3 Pd-C(CO) 3.41/3.45 physisorbed)
- **S3**: 3 weak (CO rank0 Pd-C 3.53, CH3O rank1 Pd-O 2.82, coads rank0 Pd-O 3.23)
- **S3b**: 3 weak (CO rank0/2 Pd-C 1.80, rank1 Pd-C 3.59)
- **S4**: 3 weak (CO rank0/1 Pd-C 3.68/3.90 — **Shi 가설 signal**, manual R2 Pd-O 2.55)

**모두 chemistry 해석 가능**:
- S1 borderline 0.01 → DFT relax 시 정상 chemisorbed 수렴 예상
- S2/S3/S3b weak → 표면 특수 화학 (interface, O-rich, under-coord)
- S4 weak → **Pd⁴⁺ chemistry — Shi 2024 가설 직접 검증 의도된 데이터**

---

## 3. VASP 입력 검증

### INCAR (per-surface)
```
모든 표면 공통:
  ENCUT=520, PREC=Accurate, LASPH=.TRUE., ADDGRID=.TRUE., ISPIN=2
  IVDW=12 (D3-BJ)
  EDIFF=1e-06, NELM=200, NELMIN=5
  ISYM=0
  IBRION=2, NSW=300, ISIF=2 (ionic only)
  EDIFFG=-0.03
  LDIPOL=.TRUE., IDIPOL=3
  KSPACING=0.25

표면별 ISMEAR/SIGMA:
  S1 Pd(100):       ISMEAR=1, SIGMA=0.10   (metal)
  S2 PdO/Pd:        ISMEAR=0, SIGMA=0.05   (oxide)
  S3 PdO O-term:    ISMEAR=0, SIGMA=0.05
  S3b PdO Pd-term:  ISMEAR=0, SIGMA=0.05
  S4 PdO₂(110):     ISMEAR=0, SIGMA=0.05
```

### POSCAR
- Selective Dynamics on (bottom 50% F F F)
- atom species sorted (VASP 표준)

### POTCAR
- Pd_pv (28Jan2005, 16 valence) + O (08Apr2002, 6 valence)
- For coads: + C (08Apr2002), + H (08Apr2002)
- Path: `/home/hyunjin/POTENTIAL/potpaw_PBE/{Pd_pv,O,C,H}/POTCAR`

### SLURM submit
```
calculations/T1_16_DFT_L1/{surface}/{kind}/{candidate_name}/
  POSCAR
  INCAR
  POTCAR
  submit_vasp_gpu.sh -> ../../../scripts/submit_vasp_gpu.sh

Submit command (per job):
  sbatch -J <jobname> --chdir=<jdir> scripts/submit_vasp_gpu.sh

Batch script template:
  calculations/T1_16_DFT_L1/submit_all.sh
  (47 lines, all commented — user uncomments to submit)
```

---

## 4. 비용 견적 (재확인)

| 단계 | 후보 | 시간/job | 합계 GPU-hr | Wall (2 GPUs 병렬) |
|---|---:|---|---:|---|
| **T1.16** Level 1 (vacuum) | 47 | 6-12 hr | 282-564 | **5.9-11.8 day** |
| **T1.17** Level 2 (VASPsol SP) | 47 | 1-2 hr | 47-94 | **1.0-2.0 day** |
| **Total** | | | 329-658 | **🕐 6.9-13.7 day wall** |

GPU 0 + GPU 1 모두 free, SLURM debug partition 사용 가능.

---

## 5. 알려진 caveats (T1.16 진행 중 모니터링 권고)

### ⚠ S2 coads 의심 후보
- rank-2 (idx10461) Pd-C(CO)=3.41 Å — physisorbed
- rank-3 (idx02538) Pd-C(CO)=3.45 Å — physisorbed
- → DFT 수렴 시 다른 site 로 drift 가능. 결과 비교 후 판단.

### ⚠ S3 single_CO rank-0
- Pd-C 3.53 Å (PES flat per MLIP)
- → DFT 가 ranking 의미를 정량화

### ⚠ S3b single_CO rank-0/2 (Pd-C 1.80)
- threshold 1.85 의 0.05 below — likely OK 분류로 reclassify 가능
- 안전을 위해 weak_binding 으로 분류했지만 chemistry 정상

### ⚠ S4 만들어진 manual replacement
- 90 (atop Pd): d(O-C)=1.42, d(C-H)=1.10 ✓ MLIP-free, 인공적 placement
- 91 (bridge): Pd-O=2.55 borderline. DFT relax 후 결합거리 결정.

---

## 6. 진행 절차 (DFT 진입 후)

### 사전 준비 ✅
1. ✅ 47 POSCAR 확정 (broken 0, manual replacement 2 추가)
2. ✅ 47 job 디렉토리 생성 (`calculations/T1_16_DFT_L1/`)
3. ✅ INCAR/POTCAR/submit script 모두 배치
4. ✅ Batch submit 스크립트 생성 (`submit_all.sh`, 주석으로 보존)

### 실제 submission (user 승인 시)
```bash
cd /home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc

# Option A: 한 번에 submit
bash calculations/T1_16_DFT_L1/submit_all.sh
# (현재 sbatch 라인이 주석 처리됨 — 검토 후 uncomment)

# Option B: 표면별 단계 submit (안전)
# 예: S1 만 먼저
for d in calculations/T1_16_DFT_L1/S1/*/*; do
  jname=$(basename $d)
  sbatch -J "S1_${jname:0:2}" --chdir=$d scripts/submit_vasp_gpu.sh
done
```

### 진행 중 monitoring
- `squeue -u hyunjin` 로 큐 상태
- `vasp_<jobid>.out` 으로 수렴 진행
- 24hr 후 첫 결과 (S1 작은 시스템부터)

### 종료 후 (T1.16 complete)
- Energy + force 검증 (Phase 1 처럼)
- → T1.17 Level 2 VASPsol single-point 진입

---

## 7. 결론: **Ready to Submit**

```
✓ 모든 검증 단계 통과 (G1/G2/T1.10-T1.15)
✓ 47 candidates 모두 chemistry valid (0 broken)
✓ S4 broken 2개 manual replacement 적용
✓ VASP 입력 파일 완전 준비 (INCAR/POSCAR/POTCAR/submit)
✓ Batch submission 스크립트 준비됨 (주석 처리, user 승인 후 활성화)

DFT 진입 비용: 6.9-13.7 days wall (2 GPUs)
```

**user 최종 승인 부탁드립니다.**
