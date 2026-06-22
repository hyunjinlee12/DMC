"""T1.14 Phase 3 — MACE-MH ranking for TS guess (SetTS) + thermo reference (SetB).

5 surfaces × SetTS (~1.2k/surface) + SetB subsample (~500/surface).

Why split from Phase 2:
- SetTS is the *transition-state initial guess* pool: C_CO ↔ O_CH3O = 1.7–2.3 Å.
  Used to seed T2.5 CI-NEB. Ranking by E_MACE is informative but treat with
  extra caution (TS region is poorly trained for MACE foundation models).
- SetB is the *thermodynamic reference* pool: distance ≥ 5.0 Å (non-interacting).
  Only a handful is needed to set the zero-point for descriptor map.

Pipeline identical to Phase 2 (reattach G2 slab, fix bottom 50%, LBFGS relax
with mh-1 + oc20_usemppbe + cueq + float64, dedup, rank).

SetB subsampling: random 100/surface (stratified by site_type if possible) —
plenty for descriptor map's non-interacting baseline.

Usage:
  conda run -n pddmc python scripts/run_mace_phase3.py --dry-run
  conda run -n pddmc python scripts/run_mace_phase3.py
  conda run -n pddmc python scripts/run_mace_phase3.py --surfaces S1
"""
import argparse, json, os, random, time, warnings
from pathlib import Path

# Under SLURM, GPU is assigned via --gres / cgroup — do NOT override.
# Outside SLURM (manual nohup), use GPU 1 (GPU 0 reserved for VASP).
if 'SLURM_JOB_ID' not in os.environ:
    os.environ['CUDA_VISIBLE_DEVICES'] = '1'
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.constraints import FixAtoms
from ase.optimize import LBFGS

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab'
G3 = ROOT / 'calculations/G3_adsorption'
RANDOM_SEED = 42
N_SETB_PER_SURFACE = 100

SURFACES = {
    'S1':  'S1_Pd100',
    'S2':  'S2_PdO101_Pd100',
    'S3':  'S3_PdO100',
    'S3b': 'S3b_PdO100_PdOterm',
    'S4':  'S4_PdO2_110',
}
N_ADS = 7
IDX_C_CO = 0
IDX_O_CH3O = 2


def reattach(cand, slab):
    full = slab + cand[-N_ADS:]
    full.set_cell(slab.cell); full.set_pbc(slab.pbc)
    return full


def constrain_bottom_half(atoms, n_sub):
    z = atoms.positions[:n_sub, 2]
    z_med = np.median(z)
    fixed = [i for i in range(n_sub) if atoms.positions[i, 2] < z_med]
    atoms.set_constraint(FixAtoms(indices=fixed))


def fingerprint(atoms):
    p = atoms.positions[-N_ADS:]
    dists = []
    for i in range(N_ADS):
        for j in range(i + 1, N_ADS):
            dists.append(round(float(np.linalg.norm(p[i] - p[j])), 2))
    return tuple(sorted(dists))


def metrics(atoms, n_sub):
    d_all = atoms.get_all_distances(mic=True)
    sub = np.arange(n_sub)
    ads = np.arange(n_sub, n_sub + N_ADS)
    return {'d_min': float(d_all[np.ix_(ads, sub)].min()),
            'd_reactive': float(d_all[n_sub + IDX_C_CO, n_sub + IDX_O_CH3O]),
            'd_z': float(atoms.positions[-N_ADS:, 2].min() - atoms.positions[:n_sub, 2].max())}


def load_calculator():
    from mace.calculators import mace_mp
    # D3 dispersion ENABLED per advisor 이태훈 (2026-06-16): 분자 다루므로 필요.
    # PBE+D3(BJ) — DFT INCAR 설정과 동일.
    return mace_mp(model='mh-1', head='oc20_usemppbe',
                   default_dtype='float64', enable_cueq=True,
                   device='cuda',
                   dispersion=True, damping='bj', dispersion_xc='pbe')


def relax_one(cand, slab, n_sub, calc):
    atoms = reattach(cand, slab)
    constrain_bottom_half(atoms, n_sub)
    atoms.calc = calc
    opt = LBFGS(atoms, logfile=None)
    converged = opt.run(fmax=0.05, steps=300)   # bumped 200→300 per committee S4/S3 convergence concern
    E = float(atoms.get_potential_energy())
    m = metrics(atoms, n_sub)
    return atoms, {
        'E': E, 'converged': bool(converged), 'n_steps': int(opt.nsteps),
        'd_min': m['d_min'], 'd_reactive': m['d_reactive'], 'd_z': m['d_z'],
        'fingerprint': fingerprint(atoms),
    }


def process_pool(sid, sdir, pool_name, src_traj, calc, log_fp,
                  subsample_n=None) -> dict:
    """pool_name in {'SetTS', 'SetB'}."""
    out_dir = G3 / sdir / 'MLIP_phase3'
    out_dir.mkdir(exist_ok=True, parents=True)
    summary_file = out_dir / f'summary_{pool_name}.json'
    if summary_file.exists():
        return json.load(open(summary_file))

    slab = read(G2 / sdir / 'CONTCAR')
    n_sub = len(slab)
    cands = list(read(src_traj, index=':'))
    if subsample_n is not None and len(cands) > subsample_n:
        rng = random.Random(RANDOM_SEED)
        cands = rng.sample(cands, subsample_n)
    print(f'[{sid} {pool_name}] {len(cands)} candidates')

    relaxed_traj = out_dir / f'relaxed_{pool_name}.traj'
    records = []
    tw = Trajectory(str(relaxed_traj), 'w')

    t_start = time.time()
    for k, c in enumerate(cands):
        atoms, rec = relax_one(c, slab, n_sub, calc)
        rec['idx'] = k
        records.append(rec)
        atoms.info.update({'idx': k, 'E_MACE': rec['E']})
        tw.write(atoms)
        if (k + 1) % 50 == 0 or k == 0:
            elapsed = time.time() - t_start
            eta = elapsed / (k + 1) * (len(cands) - k - 1)
            msg = f'[{sid} {pool_name}] {k+1}/{len(cands)}  ETA={eta/3600:.2f}hr'
            print(msg); log_fp.write(msg + '\n'); log_fp.flush()
    tw.close()
    t_total = time.time() - t_start

    records.sort(key=lambda r: r['E'])
    if records:
        E0 = records[0]['E']
        for r in records:
            r['dE_rel_meV'] = (r['E'] - E0) * 1000.0
    json.dump(records, open(out_dir / f'ranked_{pool_name}.json', 'w'), indent=1)

    # dedup
    seen = {}
    for r in records:
        key = (round(r['E'], 2), r['fingerprint'])
        if key not in seen:
            seen[key] = r
    unique = sorted(seen.values(), key=lambda r: r['E'])
    json.dump(unique, open(out_dir / f'unique_{pool_name}.json', 'w'), indent=1)

    # plot
    if unique:
        fig, ax = plt.subplots(figsize=(8, 5.5))
        ax.scatter([r['d_reactive'] for r in unique],
                   [r['dE_rel_meV'] for r in unique],
                   s=20, alpha=0.5, color='steelblue')
        ax.set_xlabel(r'$d(C_{CO} - O_{CH_3O})$  /  Å')
        ax.set_ylabel(r'$\Delta E_{MACE}$  /  meV')
        ax.set_title(f'{sid} {pool_name}  (unique n={len(unique)})')
        plt.tight_layout()
        plt.savefig(out_dir / f'rank_{pool_name}.png', dpi=130, bbox_inches='tight')
        plt.close()

    summary = {
        'surface': sid, 'pool': pool_name,
        'n_raw': len(cands), 'n_converged': sum(r['converged'] for r in records),
        'n_unique': len(unique),
        'time_seconds': t_total,
    }
    json.dump(summary, open(summary_file, 'w'), indent=2)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--surfaces', default='')
    ap.add_argument('--setb-n', type=int, default=N_SETB_PER_SURFACE,
                    help='SetB subsample size per surface (default 100)')
    args = ap.parse_args()

    selected = [s for s in SURFACES if (not args.surfaces or s in args.surfaces.split(','))]
    print(f'=== Phase 3 — surfaces: {selected}  SetB sub={args.setb_n}/surface ===')
    for sid in selected:
        sdir = SURFACES[sid]
        sts = G3 / sdir / 'coads_guide/SetTS.traj'
        stb = G3 / sdir / 'coads_guide/SetB.traj'
        slab = G2 / sdir / 'CONTCAR'
        n_ts = len(list(read(sts, index=':'))) if sts.exists() else 0
        n_b = len(list(read(stb, index=':'))) if stb.exists() else 0
        print(f'  {sid}  slab={slab.exists()} SetTS={n_ts} SetB={n_b}')

    if args.dry_run:
        print('dry-run done.'); return

    log_file = G3 / 'MLIP_phase3.log'
    log_fp = open(log_file, 'a')
    print('Loading MACE-MH calculator...')
    calc = load_calculator()

    summaries = []
    for sid in selected:
        sdir = SURFACES[sid]
        for pool, sub in [('SetTS', None), ('SetB', args.setb_n)]:
            src = G3 / sdir / f'coads_guide/{pool}.traj'
            if not src.exists():
                print(f'[{sid} {pool}] source missing, skip'); continue
            try:
                s = process_pool(sid, sdir, pool, src, calc, log_fp, subsample_n=sub)
                summaries.append(s)
            except Exception as e:
                msg = f'[{sid} {pool}] ERROR: {type(e).__name__}: {e}'
                print(msg); log_fp.write(msg + '\n'); log_fp.flush()

    json.dump(summaries, open(G3 / 'MLIP_phase3_summary.json', 'w'), indent=2)
    log_fp.close()
    print('=== Phase 3 done ===')
    for s in summaries:
        print(f"  {s['surface']:<4} {s['pool']:<6} raw={s['n_raw']:<4} "
              f"conv={s['n_converged']:<4} unique={s['n_unique']:<4} "
              f"time={s['time_seconds']/3600:.2f}hr")


if __name__ == '__main__':
    main()
