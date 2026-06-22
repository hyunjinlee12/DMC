"""T1.14 Phase 2 — MACE-MH ranking for co-adsorption (CO* + CH3O*) Set A reactive.

5 surfaces × ~3.8k–11.6k candidates per surface = ~37,956 relaxations.
Per-candidate cost ~2 s with cueq → expected ~22 hr total on GPU 1.

Pipeline:
  load SetA.traj candidate (substrate=trimmed, ads=last 7 [C,O,O,C,H,H,H])
  → reattach full G2 slab (read CONTCAR)
  → FixAtoms on bottom 50% of slab atoms
  → LBFGS relax (fmax=0.05 eV/Å, max 200 steps) with mh-1 + oc20_usemppbe + cueq + float64
  → record E_MACE, d_min, d_reactive (C_CO ↔ O_CH3O), n_steps, conv, fingerprint
  → save to calculations/G3_adsorption/{surface}/MLIP_phase2/

Post-relax dedup:
  group by (round(E_MACE, 2), fingerprint) — 10 meV × adsorbate distances
  keep 1 representative per group

Distance-bin stratified DFT shortlist:
  bins: [2.1, 2.5, 3.0, 3.5, 4.0] Å (reactive range)
  pick top-1 (lowest E) per bin per surface → up to 4 candidates per surface
  for S2 add extra 2 (interface-rich) → 6 candidates

NOTE: This is Phase 2 — DO NOT run until Phase 1 has been reviewed and approved
by the orchestrator. Designed for resumability (skip surfaces with existing
summary.json under MLIP_phase2/).

Usage:
  conda run -n pddmc python scripts/run_mace_phase2.py --dry-run    # validate inputs
  conda run -n pddmc python scripts/run_mace_phase2.py              # run all 5 surfaces
  conda run -n pddmc python scripts/run_mace_phase2.py --surfaces S1,S2   # subset
"""
import argparse, json, os, sys, time, warnings
from pathlib import Path
from collections import defaultdict

os.environ['CUDA_VISIBLE_DEVICES'] = '1'
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.constraints import FixAtoms
from ase.optimize import LBFGS

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab'
G3 = ROOT / 'calculations/G3_adsorption'

SURFACES = {
    'S1':  'S1_Pd100',
    'S2':  'S2_PdO101_Pd100',
    'S3':  'S3_PdO100',
    'S3b': 'S3b_PdO100_PdOterm',
    'S4':  'S4_PdO2_110',
}
N_ADS = 7   # CO (2) + CH3O (5) — order [C, O, O, C, H, H, H]
IDX_C_CO = 0     # carbon of CO (1st atom)
IDX_O_CH3O = 2   # methoxy O (3rd atom — the O bonded to CH3 C)


def reattach(cand: Atoms, slab: Atoms) -> Atoms:
    """Append last N_ADS atoms of cand (absolute xyz) to the full G2 slab."""
    ads = cand[-N_ADS:]
    full = slab + ads
    full.set_cell(slab.cell)
    full.set_pbc(slab.pbc)
    return full


def constrain_bottom_half(atoms: Atoms, n_sub: int) -> None:
    """Fix bottom 50% of substrate atoms by z."""
    sub_z = atoms.positions[:n_sub, 2]
    z_median = np.median(sub_z)
    fixed = [i for i in range(n_sub) if atoms.positions[i, 2] < z_median]
    atoms.set_constraint(FixAtoms(indices=fixed))


def fingerprint(atoms: Atoms) -> tuple:
    """Adsorbate-internal sorted distance tuple, rounded to 0.02 Å."""
    ads_pos = atoms.positions[-N_ADS:]
    n = N_ADS
    dists = []
    for i in range(n):
        for j in range(i + 1, n):
            dists.append(round(float(np.linalg.norm(ads_pos[i] - ads_pos[j])), 2))
    return tuple(sorted(dists))


def metrics(atoms: Atoms, n_sub: int) -> dict:
    """Compute ads-substrate distances and reactive C-O distance."""
    d_all = atoms.get_all_distances(mic=True)
    sub = np.arange(n_sub)
    ads = np.arange(n_sub, n_sub + N_ADS)
    d_min = float(d_all[np.ix_(ads, sub)].min())
    # closest substrate atom to each ads
    closest = {int(i - n_sub): (int(sub[d_all[i, sub].argmin()]),
                                 float(d_all[i, sub].min()))
               for i in ads}
    # reactive C_CO to O_CH3O
    d_react = float(d_all[n_sub + IDX_C_CO, n_sub + IDX_O_CH3O])
    return {'d_min': d_min, 'd_reactive': d_react,
            'd_z': float(atoms.positions[-N_ADS:, 2].min() - atoms.positions[:n_sub, 2].max())}


def load_calculator():
    from mace.calculators import mace_mp
    # D3 dispersion ENABLED per advisor 이태훈 (2026-06-16): 분자 다루므로 필요.
    # PBE+D3(BJ) — DFT INCAR 설정과 동일.
    return mace_mp(
        model='mh-1', head='oc20_usemppbe',
        default_dtype='float64', enable_cueq=True,
        device='cuda',
        dispersion=True, damping='bj', dispersion_xc='pbe',
    )


def process_surface(sid: str, sdir: str, calc, log_fp, max_n: int = None) -> dict:
    out_dir = G3 / sdir / 'MLIP_phase2'
    out_dir.mkdir(exist_ok=True, parents=True)
    summary_file = out_dir / 'summary.json'
    if summary_file.exists():
        print(f'[{sid}] summary.json exists, skipping (delete to redo)')
        return json.load(open(summary_file))

    slab = read(G2 / sdir / 'CONTCAR')
    n_sub = len(slab)

    seta_path = G3 / sdir / 'coads_guide/SetA.traj'
    if not seta_path.exists():
        print(f'[{sid}] SetA.traj missing — SKIP')
        return None
    cands = list(read(seta_path, index=':'))
    if max_n is not None:
        cands = cands[:max_n]
    print(f'[{sid}] {len(cands)} candidates  (slab={n_sub} atoms)')

    relaxed_traj = out_dir / 'relaxed_SetA.traj'
    records = []
    traj_writer = Trajectory(str(relaxed_traj), 'w')

    t_start = time.time()
    for k, cand in enumerate(cands):
        atoms = reattach(cand, slab)
        constrain_bottom_half(atoms, n_sub)
        atoms.calc = calc
        opt = LBFGS(atoms, logfile=None)
        converged = opt.run(fmax=0.05, steps=300)   # bumped 200→300 per committee S4/S3 convergence concern
        E = float(atoms.get_potential_energy())
        m = metrics(atoms, n_sub)
        fp = fingerprint(atoms)
        rec = {
            'idx': k, 'E': E, 'converged': bool(converged),
            'n_steps': int(opt.nsteps),
            'd_min': m['d_min'], 'd_reactive': m['d_reactive'], 'd_z': m['d_z'],
            'fingerprint': fp,
        }
        records.append(rec)
        atoms.info.update({'idx': k, 'E_MACE': E, 'd_min': m['d_min'],
                           'd_reactive': m['d_reactive']})
        traj_writer.write(atoms)
        if (k + 1) % 100 == 0 or k == 0:
            elapsed = time.time() - t_start
            rate = elapsed / (k + 1)
            eta = rate * (len(cands) - k - 1)
            msg = f'[{sid}] {k+1}/{len(cands)}  rate={rate:.2f}s/struct  ETA={eta/3600:.2f}hr'
            print(msg)
            log_fp.write(msg + '\n'); log_fp.flush()
    traj_writer.close()
    t_total = time.time() - t_start

    # rank by E ascending
    records.sort(key=lambda r: r['E'])
    E_min = records[0]['E']
    for r in records:
        r['dE_rel_meV'] = (r['E'] - E_min) * 1000.0

    json.dump(records, open(out_dir / 'ranked_SetA.json', 'w'), indent=1)

    # dedup
    seen = {}
    for r in records:
        key = (round(r['E'], 2), r['fingerprint'])
        if key not in seen:
            seen[key] = r
    unique = sorted(seen.values(), key=lambda r: r['E'])
    json.dump(unique, open(out_dir / 'unique_SetA.json', 'w'), indent=1)

    # distance-bin stratified shortlist
    bins = [(2.1, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.0)]
    shortlist = []
    for lo, hi in bins:
        in_bin = [r for r in unique if lo <= r['d_reactive'] < hi]
        if in_bin:
            shortlist.append(in_bin[0])
    json.dump(shortlist, open(out_dir / 'shortlist_SetA.json', 'w'), indent=1)

    # plot: E vs d_reactive scatter
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.scatter([r['d_reactive'] for r in unique],
               [r['dE_rel_meV'] for r in unique],
               s=20, alpha=0.5, color='steelblue', label=f'unique (n={len(unique)})')
    ax.scatter([r['d_reactive'] for r in shortlist],
               [r['dE_rel_meV'] for r in shortlist],
               s=80, color='red', edgecolor='black', zorder=5,
               label=f'shortlist (n={len(shortlist)})')
    ax.set_xlabel(r'$d(C_{CO} - O_{CH_3O})$  /  Å')
    ax.set_ylabel(r'$\Delta E_{MACE}$  /  meV')
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'shortlist_SetA.png', dpi=130, bbox_inches='tight')
    plt.close()

    summary = {
        'surface': sid, 'sdir': sdir,
        'n_raw': len(cands), 'n_records': len(records),
        'n_converged': sum(1 for r in records if r['converged']),
        'n_unique_after_dedup': len(unique),
        'n_shortlist': len(shortlist),
        'E_range_meV': records[-1]['dE_rel_meV'],
        'time_seconds': t_total,
        'time_per_struct': t_total / max(len(records), 1),
        'mace_config': {'model': 'mh-1', 'head': 'oc20_usemppbe',
                        'dtype': 'float64', 'cueq': True},
    }
    json.dump(summary, open(summary_file, 'w'), indent=2)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true',
                    help='validate inputs only, no relaxations')
    ap.add_argument('--surfaces', default='', help='comma-separated S1,S2,...')
    ap.add_argument('--max-n', type=int, default=None,
                    help='cap candidates per surface (debugging)')
    args = ap.parse_args()

    selected = [s for s in SURFACES if (not args.surfaces or s in args.surfaces.split(','))]
    print(f'=== Phase 2 (co-ads SetA) — surfaces: {selected} ===')
    for sid in selected:
        sdir = SURFACES[sid]
        seta = G3 / sdir / 'coads_guide/SetA.traj'
        slab = G2 / sdir / 'CONTCAR'
        ok = seta.exists() and slab.exists()
        try:
            n = len(list(read(seta, index=':'))) if ok else 0
        except Exception:
            n = -1
        print(f'  {sid}  slab={slab.exists()} SetA={seta.exists()} n_cands={n}')

    if args.dry_run:
        print('dry-run done.')
        return

    log_file = G3 / 'MLIP_phase2.log'
    log_fp = open(log_file, 'a')
    log_fp.write(f'\n===== Phase 2 start {time.strftime("%Y-%m-%d %H:%M:%S")} =====\n')
    log_fp.flush()
    print(f'Loading MACE-MH calculator (mh-1 + oc20_usemppbe + cueq)...')
    calc = load_calculator()

    global_summary = []
    for sid in selected:
        sdir = SURFACES[sid]
        try:
            s = process_surface(sid, sdir, calc, log_fp, max_n=args.max_n)
            if s: global_summary.append(s)
        except Exception as e:
            msg = f'[{sid}] ERROR: {type(e).__name__}: {e}'
            print(msg); log_fp.write(msg + '\n'); log_fp.flush()
            continue

    json.dump(global_summary, open(G3 / 'MLIP_phase2_summary.json', 'w'), indent=2)
    log_fp.close()
    print('=== Phase 2 done ===')
    for s in global_summary:
        print(f"  {s['surface']:<4} raw={s['n_raw']:<6} conv={s['n_converged']:<6} "
              f"unique={s['n_unique_after_dedup']:<5} short={s['n_shortlist']} "
              f"time={s['time_seconds']/3600:.1f}hr")


if __name__ == '__main__':
    main()
