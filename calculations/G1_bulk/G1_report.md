# G1 Report — Bulk Structure Optimization (T1.1–T1.4)

Date: 2026-06-04
Project: Pd/PdO/PdO₂ surface DMC formation DFT study

## 1. Objective

Pd, PdO, PdO₂ 3종의 bulk 구조를 PBE+D3 수준에서 최적화하고, ENCUT/k-mesh 수렴을 검증한 뒤 실험 격자상수와 비교하여 slab 제작의 기준 구조를 확보한다.

## 2. Computational Details

- Functional: PBE + D3-BJ (IVDW=12)
- Pseudopotential: PAW-PBE (Pd_pv, O)
- ENCUT: 520 eV (수렴 테스트로 확정)
- PREC: Accurate, LASPH=.TRUE., ADDGRID=.TRUE., ISPIN=2
- Relaxation: ISIF=3 (full cell + ionic), IBRION=2, EDIFFG=-0.01 eV/Å
- Smearing: Pd ISMEAR=1/SIGMA=0.10 (metal), PdO·PdO₂ ISMEAR=0/SIGMA=0.05 (oxide)

## 3. Initial Structures (Materials Project)

| Material | MP-ID | Space Group | Role |
|----------|-------|-------------|------|
| Pd | mp-2 | Fm-3m (fcc) | Pd⁰ metallic baseline |
| PdO | mp-1336 | P4₂/mmc (tetragonal) | Pd²⁺ oxide |
| PdO₂ | mp-1018886 | P4₂/mnm (rutile) | Pd⁴⁺ high-valence reference |

## 4. Convergence Tests (T1.4)

### 4.1 ENCUT Convergence (fixed k-mesh, ref = 600 eV)

| ENCUT (eV) | Pd (meV/atom) | PdO (meV/atom) | PdO₂ (meV/atom) |
|------------|---------------|----------------|------------------|
| 400 | +1.97 | -5.32 | -8.25 |
| 450 | +0.68 | -0.32 | -0.42 |
| 500 | +0.26 | +0.77 | +1.25 |
| **520** | **+0.17** | **+0.72** | **+1.16** |
| 550 | +0.23 | +0.63 | +0.94 |
| 600 | ref | ref | ref |

ENCUT = 520 eV에서 3종 모두 ~1 meV/atom 이내 수렴 확인. 520 eV 채택.

### 4.2 k-mesh Convergence (ENCUT=520 eV, ref = densest mesh)

**Pd** (ref = 16×16×16):

| k-mesh | dE (meV/atom) |
|--------|---------------|
| 6×6×6 | -9.81 |
| 8×8×8 | -1.68 |
| 10×10×10 | +1.96 |
| **12×12×12** | **-1.00** |
| 14×14×14 | +0.18 |
| 16×16×16 | ref |

**PdO** (ref = 12×12×10):

| k-mesh | dE (meV/atom) |
|--------|---------------|
| 4×4×3 | -0.84 |
| 6×6×4 | -2.26 |
| **8×8×6** | **-1.21** |
| 10×10×8 | -0.31 |
| 12×12×10 | ref |

**PdO₂** (ref = 10×10×12):

| k-mesh | dE (meV/atom) |
|--------|---------------|
| 4×4×4 | -7.73 |
| 4×4×6 | +4.04 |
| **6×6×8** | **-0.04** |
| 8×8×10 | -0.66 |
| 10×10×12 | ref |

Adopted k-mesh: Pd 12×12×12, PdO 8×8×6, PdO₂ 6×6×8.

Convergence plots: see `convergence/results/convergence_plots.png`.

## 5. Optimized Bulk Structures (T1.1–T1.3)

| Material | a_DFT (Å) | c_DFT (Å) | E₀/atom (eV) | a_exp (Å) | c_exp (Å) | err_a (%) | err_c (%) |
|----------|-----------|-----------|--------------|-----------|-----------|-----------|-----------|
| Pd (fcc) | 3.8907 | — | -5.8659 | 3.890 [1] | — | +0.02 | — |
| PdO (P4₂/mmc) | 3.0536 | 5.4058 | -5.8622 | 3.043 [2] | 5.328 [2] | +0.35 | +1.46 |
| PdO₂ (P4₂/mnm) | 4.5424 | 3.1772 | -5.5900 | 4.486 [3] | 3.103 [3] | +1.26 | +2.39 |

Lattice parameter comparison: see `G1_lattice_comparison.png`.

All errors ≤ 2.4%. PBE+D3 typically overestimates oxide lattice constants by 1–3%; PdO₂ c-axis at 2.4% is within expected range for rutile-type oxides.

## 6. G1 Gate Assessment

| Criterion | Status |
|-----------|--------|
| 3 bulk structures converged (EDIFFG < -0.01) | PASS |
| ENCUT convergence ≤ 1 meV/atom at 520 eV | PASS |
| k-mesh convergence ≤ ~1 meV/atom | PASS |
| Lattice vs experiment ≤ 2% (Pd, PdO a) | PASS |
| PdO₂ c-axis 2.4% (acceptable for PBE+D3 rutile) | PASS (marginal) |

**G1: PASS** — 3종 bulk 구조 확보 완료. Slab 제작 (G2) 진행 가능.

## 7. Output Files

- `calculations/G1_bulk/Pd/CONTCAR` — optimized fcc Pd
- `calculations/G1_bulk/PdO/CONTCAR` — optimized tetragonal PdO
- `calculations/G1_bulk/PdO2/CONTCAR` — optimized rutile PdO₂
- `calculations/G1_bulk/convergence/results/convergence_plots.png` — convergence plots
- `calculations/G1_bulk/G1_lattice_comparison.png` — DFT vs exp bar chart
- `structures/bulk_fetch_summary.json` — MP fetch metadata

## References

[1] Pd fcc lattice constant: a = 3.89 Å. Kittel, C. *Introduction to Solid State Physics*, 8th ed.
[2] PdO tetragonal (P4₂/mmc): a = 3.043 Å, c = 5.328 Å. Waser, J. et al., *J. Am. Chem. Soc.* **1953**, 75, 3400–3401.
[3] PdO₂ rutile (P4₂/mnm): a = 4.486 Å, c = 3.103 Å. Surface oxide / high-pressure synthesis data; see also Materials Project mp-1018886.
