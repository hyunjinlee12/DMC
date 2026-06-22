"""Systematic comparison: MACE-MH+D3 vs SevenNet-Omni+D3.

Tests on representative structures from all 5 surfaces (CO* + CH3O*).
Measures: relax rate (ms/step), final E (eV), final d_min ads-substrate (Å).

USAGE — run in two passes due to e3nn version conflict:
  $ # Pass 1: MACE (needs e3nn 0.4.4)
  $ conda run -n pddmc pip install --quiet "e3nn==0.4.4"
  $ conda run -n pddmc python scripts/mlip_compare_benchmark.py --mace
  $
  $ # Pass 2: SevenNet (needs e3nn 0.6.0)
  $ conda run -n pddmc pip install --quiet "e3nn==0.6.0"
  $ conda run -n pddmc python scripts/mlip_compare_benchmark.py --sevenet
  $
  $ # Pass 3: plot from both
  $ conda run -n pddmc python scripts/mlip_compare_benchmark.py --plot

Output: /home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/mlip_compare/
"""
import argparse, json, os, time, warnings
from pathlib import Path

os.environ['CUDA_VISIBLE_DEVICES'] = '1'
warnings.filterwarnings('ignore')

import numpy as np
from ase.io import read
from ase.optimize import LBFGS

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'   # noD3 also OK; we use raw heuristic candidates
G2 = ROOT / 'calculations/G2_slab'
OUT = ROOT / 'reports/mlip_compare'
OUT.mkdir(parents=True, exist_ok=True)

SURFACES = {
    'S1':  ('S1_Pd100',            'Pd(100)'),
    'S2':  ('S2_PdO101_Pd100',     'PdO(101)/Pd(100)'),
    'S3':  ('S3_PdO100',           'PdO(100) O-term'),
    'S3b': ('S3b_PdO100_PdOterm',  'PdO(100) PdO-term'),
    'S4':  ('S4_PdO2_110',         'PdO2(110)'),
}
N_STEPS = 80   # enough to reach steady state for binding & timing


def get_test_atoms(sid, sdir, ads):
    """Pick the FIRST candidate (deterministic) from each (surface, adsorbate) pool.
    For CO*/CH3O*, candidates.traj already contains full slab."""
    fpath = G3 / sdir / ads / 'candidates.traj'
    if not fpath.exists():
        # fall back to noD3 backup
        fpath = ROOT / 'calculations/G3_adsorption_noD3' / sdir / ads / 'candidates.traj'
    return read(fpath, index=0)


def n_ads_for(ads):
    return 2 if ads == 'CO' else 5


def relax_and_measure(atoms, calc, n_ads, n_steps=N_STEPS):
    a = atoms.copy()
    a.calc = calc
    _ = a.get_potential_energy()   # warmup, kernel JIT
    t0 = time.time()
    opt = LBFGS(a, logfile=None)
    opt.run(fmax=0.05, steps=n_steps)
    t_relax = time.time() - t0
    E = float(a.get_potential_energy())
    sub = list(range(len(a) - n_ads))
    d_all = a.get_all_distances(mic=True)
    d_min = float(min(d_all[len(a) - n_ads, j] for j in sub))
    return {
        'E_relaxed': E,
        'd_min_relaxed': d_min,
        'rate_ms_step': t_relax / max(opt.nsteps, 1) * 1000,
        'n_atoms': len(a),
        'n_steps_taken': int(opt.nsteps),
        'converged': bool(opt.converged()),
    }


def run_mace():
    from mace.calculators import mace_mp
    print('Loading MACE-MH + D3 ...')
    calc = mace_mp(model='mh-1', head='oc20_usemppbe',
                   default_dtype='float64', enable_cueq=True, device='cuda',
                   dispersion=True, damping='bj', dispersion_xc='pbe')
    results = {}
    for sid, (sdir, _) in SURFACES.items():
        for ads in ['CO', 'CH3O']:
            atoms = get_test_atoms(sid, sdir, ads)
            print(f'  [{sid} {ads}] n={len(atoms)}')
            r = relax_and_measure(atoms, calc, n_ads_for(ads))
            print(f'    E={r["E_relaxed"]:.3f} eV  d_min={r["d_min_relaxed"]:.2f} Å  rate={r["rate_ms_step"]:.0f} ms/step')
            results[f'{sid}_{ads}'] = r
    json.dump(results, open(OUT / 'mace_results.json', 'w'), indent=2)
    print(f'Saved {OUT}/mace_results.json')


def run_sevenet():
    from sevenn.calculator import SevenNetCalculator
    from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
    from ase.calculators.mixing import SumCalculator
    import torch
    print('Loading SevenNet-Omni + D3 ...')
    base = SevenNetCalculator(model='7net-omni', modal='oc20',
                              device='cuda', enable_cueq=True, enable_flash=False)
    d3 = TorchDFTD3Calculator(device='cuda', damping='bj', xc='pbe', dtype=torch.float64)
    calc = SumCalculator([base, d3])
    results = {}
    for sid, (sdir, _) in SURFACES.items():
        for ads in ['CO', 'CH3O']:
            atoms = get_test_atoms(sid, sdir, ads)
            print(f'  [{sid} {ads}] n={len(atoms)}')
            r = relax_and_measure(atoms, calc, n_ads_for(ads))
            print(f'    E={r["E_relaxed"]:.3f} eV  d_min={r["d_min_relaxed"]:.2f} Å  rate={r["rate_ms_step"]:.0f} ms/step')
            results[f'{sid}_{ads}'] = r
    json.dump(results, open(OUT / 'sevenet_results.json', 'w'), indent=2)
    print(f'Saved {OUT}/sevenet_results.json')


def make_plots():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({'font.size': 13})

    mace = json.load(open(OUT / 'mace_results.json'))
    seve = json.load(open(OUT / 'sevenet_results.json'))

    keys = list(mace.keys())   # S1_CO, S1_CH3O, S2_CO, ...
    n_atoms = [mace[k]['n_atoms'] for k in keys]
    sids = [k.split('_')[0] for k in keys]
    ads = [k.split('_', 1)[1] for k in keys]

    rate_m = [mace[k]['rate_ms_step'] for k in keys]
    rate_s = [seve[k]['rate_ms_step'] for k in keys]
    dmin_m = [mace[k]['d_min_relaxed'] for k in keys]
    dmin_s = [seve[k]['d_min_relaxed'] for k in keys]
    E_m = [mace[k]['E_relaxed'] for k in keys]
    E_s = [seve[k]['E_relaxed'] for k in keys]

    # =====  Figure 1: speed scaling  =====
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.scatter(n_atoms, rate_m, s=100, color='#1f4e79', label='MACE-MH + D3 (cueq)', marker='o', edgecolor='black')
    ax.scatter(n_atoms, rate_s, s=100, color='#e76f51', label='SevenNet-Omni + D3 (cueq)', marker='s', edgecolor='black')
    for x, y, lbl in zip(n_atoms, rate_m, keys):
        ax.annotate(lbl.replace('_', '\n'), (x, y), fontsize=9, alpha=0.6, ha='center', va='bottom')
    ax.set_xlabel('n_atoms (slab + adsorbate)')
    ax.set_ylabel('Relax rate / (ms/step)')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT / 'fig1_speed_scaling.png', dpi=140, bbox_inches='tight')
    plt.close()

    # =====  Figure 2: binding distance (chemisorption check) =====
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(keys))
    w = 0.35
    ax.bar(x - w/2, dmin_m, w, label='MACE-MH + D3', color='#1f4e79', edgecolor='black')
    ax.bar(x + w/2, dmin_s, w, label='SevenNet-Omni + D3', color='#e76f51', edgecolor='black')
    ax.axhspan(1.7, 2.3, color='green', alpha=0.18, label='Pd-C chemisorbed (1.7–2.3 Å)')
    ax.axhline(3.0, color='red', linestyle='--', alpha=0.4, label='physisorption threshold (~3 Å)')
    ax.set_xticks(x); ax.set_xticklabels([k.replace('_', '\n') for k in keys], fontsize=10)
    ax.set_ylabel(r'$d_{\min}$ adsorbate-substrate / Å')
    ax.legend(loc='upper left', fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(OUT / 'fig2_binding_distance.png', dpi=140, bbox_inches='tight')
    plt.close()

    # =====  Figure 3: side-by-side rate/d_min  =====
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    ax.scatter(n_atoms, rate_m, s=120, color='#1f4e79', label='MACE-MH + D3', marker='o', edgecolor='black')
    ax.scatter(n_atoms, rate_s, s=120, color='#e76f51', label='SevenNet-Omni + D3', marker='s', edgecolor='black')
    # connect with lines for visual scaling
    order = np.argsort(n_atoms)
    ax.plot(np.array(n_atoms)[order], np.array(rate_m)[order], '--', color='#1f4e79', alpha=0.5)
    ax.plot(np.array(n_atoms)[order], np.array(rate_s)[order], '--', color='#e76f51', alpha=0.5)
    ax.set_xlabel('n_atoms')
    ax.set_ylabel('Relax rate / (ms/step)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Speed scaling')

    ax = axes[1]
    ax.bar(x - w/2, dmin_m, w, label='MACE-MH + D3', color='#1f4e79', edgecolor='black')
    ax.bar(x + w/2, dmin_s, w, label='SevenNet-Omni + D3', color='#e76f51', edgecolor='black')
    ax.axhspan(1.7, 2.3, color='green', alpha=0.18, label='chemisorbed band')
    ax.axhline(3.0, color='red', linestyle='--', alpha=0.4, label='physisorption (~3 Å)')
    ax.set_xticks(x); ax.set_xticklabels([k.replace('_', '\n') for k in keys], fontsize=10)
    ax.set_ylabel(r'$d_{\min}$ / Å')
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_title('Binding (chemisorption) check')
    plt.tight_layout()
    plt.savefig(OUT / 'fig3_combined.png', dpi=140, bbox_inches='tight')
    plt.close()

    # Combined summary CSV
    import csv
    with open(OUT / 'summary.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['key','n_atoms','MACE_E_eV','MACE_dmin_A','MACE_rate_ms','SevenNet_E_eV','SevenNet_dmin_A','SevenNet_rate_ms','speedup_MACE_over_Sevenet'])
        for k, n in zip(keys, n_atoms):
            sp = seve[k]['rate_ms_step'] / mace[k]['rate_ms_step']
            w.writerow([k, n,
                       f'{mace[k]["E_relaxed"]:.3f}', f'{mace[k]["d_min_relaxed"]:.2f}', f'{mace[k]["rate_ms_step"]:.0f}',
                       f'{seve[k]["E_relaxed"]:.3f}', f'{seve[k]["d_min_relaxed"]:.2f}', f'{seve[k]["rate_ms_step"]:.0f}',
                       f'{sp:.2f}'])
    print(f'Saved 3 figures + summary.csv → {OUT}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mace', action='store_true')
    ap.add_argument('--sevenet', action='store_true')
    ap.add_argument('--plot', action='store_true')
    args = ap.parse_args()
    if args.mace: run_mace()
    if args.sevenet: run_sevenet()
    if args.plot: make_plots()
    if not any([args.mace, args.sevenet, args.plot]):
        print(__doc__)


if __name__ == '__main__':
    main()
