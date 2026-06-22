#!/usr/bin/env python3
"""
T1.14: MACE-MH ranking of single-adsorbate candidates (CO*, CH3O*)
Ranks ~2,516 heuristic candidates and selects DFT shortlist (top-3 per surface per adsorbate).
RANKING ONLY — never trust absolute MACE energies; DFT is ground truth.
"""
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '1'  # GPU 0 occupied by VASP job 334

import json
import time
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from mace.calculators import mace_mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Project paths
BASE = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2_SLAB = BASE / 'calculations/G2_slab'
G3_ADS = BASE / 'calculations/G3_adsorption'
LOG_FILE = G3_ADS / 'MLIP_phase1.log'

SURFACES = ['S1_Pd100', 'S2_PdO101_Pd100', 'S3_PdO100', 'S3b_PdO100_PdOterm', 'S4_PdO2_110']
ADSORBATES = ['CO', 'CH3O']

# MACE calculator (project-confirmed settings + D3)
def get_mace_calc():
    """MACE-MH + oc20_usemppbe head + cueq + D3(BJ).
    D3 enabled per advisor 이태훈 (2026-06-16) — 분자 다루므로 필요.
    PBE+D3(BJ) — DFT INCAR 설정과 동일.
    """
    return mace_mp(
        model='mh-1',
        head='oc20_usemppbe',
        default_dtype='float64',
        enable_cueq=True,
        device='cuda',
        dispersion=True, damping='bj', dispersion_xc='pbe',
    )

def log(msg):
    """Thread-safe logging."""
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
    line = f"[{timestamp}] {msg}\n"
    print(line.strip())
    with open(LOG_FILE, 'a') as f:
        f.write(line)

def classify_site(atoms, ads_indices, surf_indices):
    """
    Best-effort site classification: atop (1 nearest), bridge (2), hollow (3+).
    Returns site_type string.
    """
    from ase.geometry import get_distances

    ads_pos = atoms.positions[ads_indices]
    surf_pos = atoms.positions[surf_indices]

    # Find min distances from each ads atom to substrate
    min_dists = []
    for ap in ads_pos:
        dists = np.linalg.norm(surf_pos - ap, axis=1)
        min_dists.append(dists.min())

    overall_min = min(min_dists)

    # Find substrate atoms within 0.5 Å buffer of closest ads-sub distance
    ads_center = ads_pos.mean(axis=0)
    dists_center = np.linalg.norm(surf_pos - ads_center, axis=1)
    threshold = overall_min + 0.5
    nearby = np.sum(dists_center < threshold)

    if nearby == 1:
        return 'atop'
    elif nearby == 2:
        return 'bridge'
    elif nearby >= 3:
        return 'hollow'
    else:
        return 'unknown'

def get_fingerprint(atoms, ads_indices):
    """
    Structural fingerprint: sorted tuple of rounded interatomic distances in adsorbate.
    Used for deduplication.
    """
    from ase.geometry import get_distances
    if len(ads_indices) == 1:
        return ()

    ads_pos = atoms.positions[ads_indices]
    dists = []
    for i in range(len(ads_indices)):
        for j in range(i+1, len(ads_indices)):
            d = np.linalg.norm(ads_pos[i] - ads_pos[j])
            dists.append(round(d, 2))
    return tuple(sorted(dists))

def relax_candidate(atoms_orig, ref_slab, calc, ads_type):
    """
    Relax a single candidate with MACE-MH.

    Returns: dict with {
        'atoms': relaxed Atoms,
        'E_MACE': float (eV),
        'converged': bool,
        'n_steps': int,
        'd_min_ads_sub': float (Å),
        'site_type': str,
        'fingerprint': tuple
    }
    """
    from ase.geometry import get_distances

    atoms = atoms_orig.copy()

    # Determine adsorbate size
    n_ads = 2 if ads_type == 'CO' else 5
    n_slab = len(ref_slab)

    # Sanity check: if substrate length mismatches G2 slab, reattach
    if len(atoms) != n_slab + n_ads:
        log(f"  WARNING: candidate length {len(atoms)} != {n_slab + n_ads}, reattaching")
        # Take adsorbate atoms (last n_ads) and attach to ref_slab
        ads_atoms = atoms[-n_ads:]
        atoms = ref_slab.copy()
        atoms.extend(ads_atoms)
        atoms.set_cell(ref_slab.cell)
        atoms.set_pbc(ref_slab.pbc)

    # Constrain bottom 2 layers (bottom 50% by z)
    z_coords = atoms.positions[:n_slab, 2]
    z_mid = (z_coords.max() + z_coords.min()) / 2
    bottom_mask = [i for i in range(n_slab) if atoms.positions[i, 2] < z_mid]
    atoms.set_constraint(FixAtoms(indices=bottom_mask))

    # Set calculator
    atoms.calc = calc

    # Relax with LBFGS
    opt = LBFGS(atoms, logfile=None, trajectory=None)
    try:
        opt.run(fmax=0.05, steps=200)
        converged = opt.converged()
        n_steps = opt.nsteps
    except Exception as e:
        log(f"  ERROR during relax: {e}")
        return None

    # Extract results
    E = atoms.get_potential_energy()

    # Min distance ads-substrate
    ads_indices = list(range(n_slab, n_slab + n_ads))
    surf_indices = list(range(n_slab))
    ads_pos = atoms.positions[ads_indices]
    surf_pos = atoms.positions[surf_indices]

    min_dist = 999.0
    for ap in ads_pos:
        dists = np.linalg.norm(surf_pos - ap, axis=1)
        min_dist = min(min_dist, dists.min())

    # Site type
    site = classify_site(atoms, ads_indices, surf_indices)

    # Fingerprint
    fp = get_fingerprint(atoms, ads_indices)

    return {
        'atoms': atoms,
        'E_MACE': E,
        'converged': converged,
        'n_steps': n_steps,
        'd_min_ads_sub': min_dist,
        'site_type': site,
        'fingerprint': fp
    }

def deduplicate(results):
    """
    Group by (round(E, 2), fingerprint) and keep one representative per group.
    Returns indices of unique structures.
    """
    groups = {}
    for i, r in enumerate(results):
        key = (round(r['E_MACE'], 2), r['fingerprint'])
        if key not in groups:
            groups[key] = []
        groups[key].append(i)

    # Pick representative (first in each group)
    unique_indices = [grp[0] for grp in groups.values()]
    return sorted(unique_indices)

def process_surface(surf, calc):
    """Process one surface: relax CO* and CH3O* candidates, rank, dedup, shortlist."""

    log(f"\n{'='*60}")
    log(f"Processing {surf}")
    log(f"{'='*60}")

    # Paths
    surf_dir = G3_ADS / surf
    mlip_dir = surf_dir / 'MLIP_phase1'
    mlip_dir.mkdir(exist_ok=True, parents=True)

    summary_path = mlip_dir / 'summary.json'
    if summary_path.exists():
        log(f"  {surf} already processed (summary.json exists), skipping")
        with open(summary_path) as f:
            return json.load(f)

    # Load G2 reference slab
    ref_slab_path = G2_SLAB / surf / 'CONTCAR'
    ref_slab = read(ref_slab_path)
    log(f"  Loaded G2 slab: {len(ref_slab)} atoms")

    summary = {'surface': surf, 'adsorbates': {}}

    for ads in ADSORBATES:
        t0 = time.time()
        log(f"\n  --- {ads} ---")

        # Load candidates
        cand_path = surf_dir / ads / 'candidates.traj'
        if not cand_path.exists():
            log(f"  WARNING: {cand_path} not found, skipping")
            continue

        candidates = read(cand_path, ':')
        n_raw = len(candidates)
        log(f"  Loaded {n_raw} candidates")

        # Relax all
        results = []
        for i, atoms in enumerate(candidates):
            if (i + 1) % 50 == 0 or i == n_raw - 1:
                log(f"    Relaxing {i+1}/{n_raw}")

            res = relax_candidate(atoms, ref_slab, calc, ads)
            if res is not None:
                res['idx'] = i
                results.append(res)

        n_conv = sum(1 for r in results if r['converged'])
        log(f"  Relaxed {len(results)}/{n_raw}, converged {n_conv}")

        if not results:
            log(f"  No successful relaxations for {ads}, skipping")
            continue

        # Sort by energy
        results.sort(key=lambda x: x['E_MACE'])

        # Relative energies
        E_min = results[0]['E_MACE']
        for r in results:
            r['dE_rel'] = (r['E_MACE'] - E_min) * 1000  # meV

        # Dedup
        unique_idx = deduplicate(results)
        unique_results = [results[i] for i in unique_idx]
        log(f"  Dedup: {len(results)} → {len(unique_results)} unique")

        # Top-3 for DFT shortlist
        top3_indices = [unique_results[i]['idx'] for i in range(min(3, len(unique_results)))]

        # Save relaxed structures
        relaxed_traj = mlip_dir / f'relaxed_{ads}.traj'
        relaxed_atoms = [r['atoms'] for r in results]
        for i, (atm, res) in enumerate(zip(relaxed_atoms, results)):
            atm.info['mace_E'] = res['E_MACE']
            atm.info['converged'] = res['converged']
            atm.info['n_steps'] = res['n_steps']
            atm.info['d_min'] = res['d_min_ads_sub']
            atm.info['site_type'] = res['site_type']
            atm.info['dE_rel_meV'] = res['dE_rel']
            atm.info['original_idx'] = res['idx']
        write(relaxed_traj, relaxed_atoms)
        log(f"  Saved {relaxed_traj}")

        # Save ranked JSON
        ranked_json = mlip_dir / f'ranked_{ads}.json'
        ranked_data = [{
            'idx': int(r['idx']),
            'E': float(r['E_MACE']),
            'dE_rel': float(r['dE_rel']),
            'converged': bool(r['converged']),
            'n_steps': int(r['n_steps']),
            'd_min': float(r['d_min_ads_sub']),
            'site_type': r['site_type'],
            'fingerprint': list(r['fingerprint'])
        } for r in results]
        with open(ranked_json, 'w') as f:
            json.dump(ranked_data, f, indent=2)
        log(f"  Saved {ranked_json}")

        # Save unique JSON
        unique_json = mlip_dir / f'unique_{ads}.json'
        unique_data = [{
            'idx': int(r['idx']),
            'E': float(r['E_MACE']),
            'dE_rel': float(r['dE_rel']),
            'converged': bool(r['converged']),
            'n_steps': int(r['n_steps']),
            'd_min': float(r['d_min_ads_sub']),
            'site_type': r['site_type'],
            'fingerprint': list(r['fingerprint'])
        } for r in unique_results]
        with open(unique_json, 'w') as f:
            json.dump(unique_data, f, indent=2)
        log(f"  Saved {unique_json}")

        # Plot top-10
        plot_path = mlip_dir / f'top10_{ads}.png'
        fig, ax = plt.subplots(figsize=(8, 5))

        # All unique structures
        x = [r['d_min_ads_sub'] for r in unique_results]
        y = [r['dE_rel'] for r in unique_results]
        ax.scatter(x, y, s=30, alpha=0.5, label='Unique candidates')

        # Highlight top-3
        if len(unique_results) >= 3:
            x_top = [unique_results[i]['d_min_ads_sub'] for i in range(3)]
            y_top = [unique_results[i]['dE_rel'] for i in range(3)]
            ax.scatter(x_top, y_top, s=100, c='red', marker='*', label='Top-3 DFT shortlist', zorder=5)

        ax.set_xlabel('d_min (ads-substrate) [Å]')
        ax.set_ylabel('ΔE_rel [meV]')
        ax.set_title(f'{surf} — {ads}* (N_unique={len(unique_results)})')
        ax.legend()
        ax.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150)
        plt.close()
        log(f"  Saved {plot_path}")

        # Summary stats
        E_range = (results[-1]['E_MACE'] - results[0]['E_MACE']) * 1000  # meV
        avg_steps = np.mean([r['n_steps'] for r in results])
        time_taken = time.time() - t0

        summary['adsorbates'][ads] = {
            'n_raw': n_raw,
            'n_converged': n_conv,
            'n_unique': len(unique_results),
            'n_top3': len(top3_indices),
            'top3_indices': top3_indices,
            'E_range_meV': round(E_range, 1),
            'avg_steps': round(avg_steps, 1),
            'time_sec': round(time_taken, 1),
            'top1_site_type': unique_results[0]['site_type'] if unique_results else None
        }

        log(f"  Stats: E_range={E_range:.1f} meV, avg_steps={avg_steps:.1f}, time={time_taken:.1f}s")

    # Save surface summary
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    log(f"\nSaved {summary_path}")

    return summary

def main():
    """Main pipeline."""

    log("="*80)
    log("T1.14: MACE-MH Phase 1 ranking — single adsorbates (CO*, CH3O*)")
    log("="*80)

    # Initialize MACE calculator
    log("\nInitializing MACE-MH calculator (mh-1 + oc20_usemppbe + cueq)...")
    t_init = time.time()
    calc = get_mace_calc()
    log(f"  Calculator ready ({time.time() - t_init:.1f}s)")

    # Process all surfaces
    t_total = time.time()
    all_summaries = []

    for surf in SURFACES:
        summary = process_surface(surf, calc)
        all_summaries.append(summary)

    # Global summary
    global_summary = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'total_time_sec': round(time.time() - t_total, 1),
        'surfaces': all_summaries
    }

    global_path = G3_ADS / 'MLIP_phase1_summary.json'
    with open(global_path, 'w') as f:
        json.dump(global_summary, f, indent=2)

    log("\n" + "="*80)
    log(f"COMPLETE. Total time: {global_summary['total_time_sec']:.1f}s")
    log(f"Global summary: {global_path}")
    log("="*80)

    # Quick report table
    log("\nPer-surface summary:")
    log(f"{'Surface':<20} {'CO_raw':>7} {'CO_conv':>8} {'CO_uniq':>8} {'CO_E_meV':>10} {'CH3O_raw':>9} {'CH3O_conv':>10} {'CH3O_uniq':>10} {'CH3O_E_meV':>12}")
    log("-"*120)

    for s in all_summaries:
        surf = s['surface']
        co = s['adsorbates'].get('CO', {})
        ch3o = s['adsorbates'].get('CH3O', {})
        log(f"{surf:<20} {co.get('n_raw', 0):>7} {co.get('n_converged', 0):>8} "
            f"{co.get('n_unique', 0):>8} {co.get('E_range_meV', 0):>10.1f} "
            f"{ch3o.get('n_raw', 0):>9} {ch3o.get('n_converged', 0):>10} "
            f"{ch3o.get('n_unique', 0):>10} {ch3o.get('E_range_meV', 0):>12.1f}")

    log("\nDFT shortlist (top-3 per surface per adsorbate):")
    for s in all_summaries:
        surf = s['surface']
        for ads in ['CO', 'CH3O']:
            if ads in s['adsorbates']:
                indices = s['adsorbates'][ads]['top3_indices']
                log(f"  {surf}/{ads}: {indices}")

if __name__ == '__main__':
    main()
