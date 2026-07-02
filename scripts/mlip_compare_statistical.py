"""Statistical MLIP comparison — 5 random samples per (surface, ads) × 2 MLIPs.

Two-pass execution due to e3nn version conflict:
  Pass 1 (--mace):    e3nn==0.4.4, MACE-MH + D3 + cueq
  Pass 2 (--sevenet): e3nn==0.6.0, SevenNet-Omni + D3
  Pass 3 (--plot):    compare distributions

Sampling: deterministic seed for reproducibility.
"""
import argparse, json, os, random, time, warnings
from pathlib import Path

os.environ['CUDA_VISIBLE_DEVICES'] = '1'
warnings.filterwarnings('ignore')

import numpy as np
from ase.io import read
from ase.optimize import LBFGS
from ase.constraints import FixAtoms

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
G2 = ROOT / 'calculations/G2_slab'
OUT = ROOT / 'reports/mlip_compare/statistical/N20'
OUT.mkdir(parents=True, exist_ok=True)

SURFACES = {
    'S1':  'S1_Pd100',
    'S2':  'S2_PdO101_Pd100',
    'S3':  'S3_PdO100',
    'S3b': 'S3b_PdO100_PdOterm',
    'S4':  'S4_PdO2_110',
}
N_SAMPLES = 20
SEED = 42


def sample_candidates(sid, sdir, ads):
    """Return N_SAMPLES deterministic random indices from candidates.traj.
    For CO: pool = candidates.traj (full slab + 2 atoms CO)
    For CH3O: candidates.traj (full slab + 5 atoms CH3O)
    """
    src = G3 / sdir / ads / 'candidates.traj'
    if not src.exists():
        src = ROOT / 'calculations/G3_adsorption_noD3' / sdir / ads / 'candidates.traj'
    cands = list(read(src, index=':'))
    rng = random.Random(SEED + hash(sid + ads) % 100)
    idxs = rng.sample(range(len(cands)), min(N_SAMPLES, len(cands)))
    return [(i, cands[i]) for i in idxs]


def fix_bottom_half(atoms, n_sub):
    z = atoms.positions[:n_sub, 2]
    z_med = np.median(z)
    fixed = [i for i in range(n_sub) if atoms.positions[i, 2] < z_med]
    atoms.set_constraint(FixAtoms(indices=fixed))


def relax_measure(atoms, calc, n_ads, n_steps=80):
    a = atoms.copy()
    a.calc = calc
    _ = a.get_potential_energy()
    t0 = time.time()
    opt = LBFGS(a, logfile=None)
    opt.run(fmax=0.05, steps=n_steps)
    t = time.time() - t0
    E = float(a.get_potential_energy())
    n = len(a)
    sub_inds = list(range(n - n_ads))
    d_all = a.get_all_distances(mic=True)
    d_min = min(d_all[n - n_ads, j] for j in sub_inds)
    return {'E': E, 'd_min': float(d_min), 'time_s': t, 'n_steps': int(opt.nsteps),
            'converged': bool(opt.converged())}


def run_mlip(mlip):
    if mlip == 'mace':
        from mace.calculators import mace_mp
        calc = mace_mp(model='mh-1', head='oc20_usemppbe',
                       default_dtype='float64', enable_cueq=True, device='cuda',
                       dispersion=True, damping='bj', dispersion_xc='pbe')
    elif mlip == 'sevenet':
        from sevenn.calculator import SevenNetCalculator
        from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
        from ase.calculators.mixing import SumCalculator
        import torch
        base = SevenNetCalculator(model='7net-omni', modal='oc20',
                                  device='cuda', enable_cueq=True, enable_flash=False)
        d3 = TorchDFTD3Calculator(device='cuda', damping='bj', xc='pbe', dtype=torch.float64)
        calc = SumCalculator([base, d3])
    else:
        raise ValueError(mlip)

    print(f'\n=== {mlip} ===')
    results = []
    for sid, sdir in SURFACES.items():
        slab = read(G2 / sdir / 'CONTCAR')
        n_sub = len(slab)
        for ads in ['CO', 'CH3O']:
            samples = sample_candidates(sid, sdir, ads)
            for k, (idx, atoms) in enumerate(samples):
                n_ads = 2 if ads == 'CO' else 5
                # Reattach if needed (single ads traj usually has full slab; safe-check)
                if len(atoms) != n_sub + n_ads:
                    ads_atoms = atoms[-n_ads:]
                    atoms = slab.copy()
                    atoms += ads_atoms
                fix_bottom_half(atoms, n_sub)
                r = relax_measure(atoms, calc, n_ads)
                r['surface'] = sid; r['ads'] = ads; r['orig_idx'] = idx; r['sample_k'] = k
                results.append(r)
                print(f'  {sid} {ads} k={k} idx={idx}  d_min={r["d_min"]:.2f}  E={r["E"]:.2f}  conv={r["converged"]}')
    json.dump(results, open(OUT / f'{mlip}_results.json', 'w'), indent=2)
    print(f'\nSaved {OUT/(mlip+"_results.json")}')


def plot_compare():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
    plt.rcParams['font.size'] = 13
    plt.rcParams['axes.labelsize'] = 15
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 11
    plt.rcParams['axes.linewidth'] = 1.2

    mace = json.load(open(OUT / 'mace_results.json'))
    seve = json.load(open(OUT / 'sevenet_results.json'))

    # Group by (surface, ads, sample_k)
    def key(r): return (r['surface'], r['ads'], r['sample_k'])
    mace_d = {key(r): r for r in mace}
    seve_d = {key(r): r for r in seve}

    surfaces_list = list(SURFACES.keys())
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

    for ax, ads in [(axes[0], 'CO'), (axes[1], 'CH3O')]:
        positions = np.arange(len(surfaces_list))
        mace_data = []
        seve_data = []
        for sid in surfaces_list:
            mace_data.append([mace_d.get((sid, ads, k), {}).get('d_min')
                              for k in range(N_SAMPLES) if (sid, ads, k) in mace_d])
            seve_data.append([seve_d.get((sid, ads, k), {}).get('d_min')
                              for k in range(N_SAMPLES) if (sid, ads, k) in seve_d])

        w = 0.35
        # MACE: blue circles (jittered to avoid overlap)
        rng = np.random.default_rng(1)
        for i, d in enumerate(mace_data):
            xs = i - w/2 + rng.normal(0, 0.03, size=len(d))
            ax.scatter(xs, d, s=40, color='#1f4e79', edgecolor='black',
                       linewidth=0.5, alpha=0.7, zorder=3)
            if len(d) > 1:
                ax.boxplot([d], positions=[i - w/2], widths=0.25, patch_artist=True,
                           boxprops=dict(facecolor='#5f7faa', alpha=0.4),
                           medianprops=dict(color='red', lw=2),
                           showfliers=False)
        # SevenNet: orange squares
        for i, d in enumerate(seve_data):
            xs = i + w/2 + rng.normal(0, 0.03, size=len(d))
            ax.scatter(xs, d, s=40, color='#e76f51', edgecolor='black', marker='s',
                       linewidth=0.5, alpha=0.7, zorder=3)
            if len(d) > 1:
                ax.boxplot([d], positions=[i + w/2], widths=0.25, patch_artist=True,
                           boxprops=dict(facecolor='#f4a261', alpha=0.4),
                           medianprops=dict(color='red', lw=2),
                           showfliers=False)
        ax.axhspan(1.85, 2.10 if ads == 'CO' else 2.15, alpha=0.15, color='green')
        ax.axhline(3.0, ls='--', color='red', alpha=0.4, lw=1.2)
        # Direct labels (no legend)
        ax.text(0.98, 0.04, 'MACE = blue circles\nSevenNet = orange squares\nGreen band = chemisorbed\nRed dashed = physisorption ~3 Å',
                transform=ax.transAxes, fontsize=9, va='bottom', ha='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='gray'))
        ax.set_xticks(positions); ax.set_xticklabels(surfaces_list)
        elem = 'C' if ads == 'CO' else 'O'
        ax.set_ylabel(r'$d_{\mathrm{min}}$ (Pd–' + elem + r') / Å')
        ax.text(0.04, 0.96, f'({"a" if ads == "CO" else "b"}) {ads}*',
                transform=ax.transAxes, fontsize=14, fontweight='bold', va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.5, len(surfaces_list)-0.5)
    plt.tight_layout()
    plt.savefig(OUT / 'F_statistical_compare.png', dpi=160, bbox_inches='tight')
    plt.close()
    print(f'Plot: {OUT/"F_statistical_compare.png"}')

    # Summary table
    print('\n=== Summary statistics ===')
    print(f"{'Surface':<5} {'ads':<6} {'MLIP':<10} {'n':<3} {'min':<6} {'max':<6} {'median':<8} {'range':<8}")
    for sid in surfaces_list:
        for ads in ['CO', 'CH3O']:
            for mlip_name, data in [('MACE', mace_d), ('SevenNet', seve_d)]:
                d = [data[(sid, ads, k)]['d_min'] for k in range(N_SAMPLES) if (sid, ads, k) in data]
                if not d: continue
                print(f"{sid:<5} {ads:<6} {mlip_name:<10} {len(d):<3} "
                      f"{min(d):<6.2f} {max(d):<6.2f} {np.median(d):<8.2f} {max(d)-min(d):<8.2f}")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--mace', action='store_true')
    p.add_argument('--sevenet', action='store_true')
    p.add_argument('--plot', action='store_true')
    args = p.parse_args()
    if args.mace: run_mlip('mace')
    if args.sevenet: run_mlip('sevenet')
    if args.plot: plot_compare()
    if not any([args.mace, args.sevenet, args.plot]):
        print(__doc__)
