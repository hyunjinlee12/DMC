"""Benchmark MACE-MH medium with/without cueq, float32/64.

Settings tested:
  baseline : medium float64, cueq=False
  cueq64   : medium float64, cueq=True       ← relax 권장
  cueq32   : medium float32, cueq=True       ← MD only (이용혁 교수 조언, 비교용)

Structure samples:
  small: S1 Pd(100) + CO*  (~82 atoms)
  large: S3 PdO(100) + CO+CH3O co-ads (~135 atoms)
Run: 30-step LBFGS each, on GPU 1.
"""
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '1'  # GPU 0 점거되어 GPU 1 사용

import time, warnings
warnings.filterwarnings('ignore')
import torch
from mace.calculators import mace_mp
from ase.io import read
from ase.optimize import LBFGS
import numpy as np


def reattach(cand, slab, n_ads):
    """Re-attach trimmed substrate."""
    ads = cand[-n_ads:]
    full = slab + ads
    full.set_pbc(slab.pbc); full.set_cell(slab.cell)
    return full


def bench(label, atoms, model, dtype, cueq, steps=30):
    a = atoms.copy()
    t0 = time.time()
    calc = mace_mp(model=model, device='cuda', default_dtype=dtype, enable_cueq=cueq, dispersion=False)
    t_load = time.time() - t0
    a.calc = calc
    # warmup
    _ = a.get_potential_energy()
    # timed
    t0 = time.time()
    E0 = a.get_potential_energy()
    F0 = a.get_forces()
    t_sp = time.time() - t0
    opt = LBFGS(a, logfile=None)
    t0 = time.time()
    opt.run(fmax=0.001, steps=steps)
    t_relax = time.time() - t0
    E_final = a.get_potential_energy()
    print(f'{label:<30} load={t_load:5.1f}s  SP={t_sp*1000:6.0f}ms  relax_{steps}step={t_relax:5.1f}s  ({t_relax/steps*1000:.0f} ms/step)')
    return {'E0': E0, 'F0_max': float(np.abs(F0).max()), 'E_final': E_final}


# Load structures
slab_S1 = read('calculations/G2_slab/S1_Pd100/CONTCAR')
slab_S3 = read('calculations/G2_slab/S3_PdO100/CONTCAR')

co_S1 = read('calculations/G3_adsorption/S1_Pd100/CO/candidates.traj', index=0)  # has full slab
coa_S3 = read('calculations/G3_adsorption/S3_PdO100/coads_guide/SetA.traj', index=0)
coa_S3 = reattach(coa_S3, slab_S3, 7)

small = co_S1
large = coa_S3
print(f"\nSmall sample : S1 Pd(100) + CO*   ({len(small)} atoms)")
print(f"Large sample : S3 PdO(100) + co-ads ({len(large)} atoms)")
print(f"GPU         : {torch.cuda.get_device_name(0)} (CUDA_VISIBLE_DEVICES=1)\n")

MODEL = 'medium'  # MACE-MP (MH) medium

# Run benchmark matrix
# (cueq 비활성: mace-torch 0.3.16 + cuequivariance API 호환 이슈)
configs = [
    ('medium_float64', 'float64', False),
    ('medium_float32', 'float32', False),
]

results = {}
print('=== Small (S1 CO*, 82 atoms) ===')
for label, dtype, cueq in configs:
    results[('small', label)] = bench(label, small, MODEL, dtype, cueq)

print()
print('=== Large (S3 co-ads, 135 atoms) ===')
for label, dtype, cueq in configs:
    results[('large', label)] = bench(label, large, MODEL, dtype, cueq)

print()
print('=== Energy / force accuracy comparison ===')
ref_small = results[('small', 'medium_float64')]
ref_large = results[('large', 'medium_float64')]
for label, dtype, cueq in configs[1:]:
    rs = results[('small', label)]
    rl = results[('large', label)]
    print(f'{label:<30} ΔE_small={(rs["E0"]-ref_small["E0"])*1000:7.3f} meV  ΔE_large={(rl["E0"]-ref_large["E0"])*1000:7.3f} meV')
