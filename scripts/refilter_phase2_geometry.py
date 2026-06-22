"""Refilter Phase 2 SetA results with post-relax geometry validator.

Catches silent PBC fragmentation errors:
  - mic=True distance metrics MASK split-across-cell adsorbates
  - This script uses DIRECT (non-MIC) inter-atomic distances to detect fragmentation

For each surface:
  1. Read relaxed_SetA.traj
  2. Apply geometry_valid(atoms) — drop fragmented / collapsed / broken-bond structures
  3. Re-rank surviving structures by E_MACE+D3 (already in info)
  4. Re-dedup with LOOSER fingerprint (10 meV E bin + 0.1 Å d_reactive bin)
  5. Re-pick distance-bin stratified shortlist
  6. Save outputs to MLIP_phase2_filtered/

Adsorbate atom order (last 7): [C, O, O, C, H, H, H]
  - idx 0-1: CO (C=O ~1.15 Å)
  - idx 2-3: O-C of methoxy (~1.40 Å)
  - idx 3-{4,5,6}: methyl C-H (~1.10 Å)
"""
import json, shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read
from ase.io.trajectory import Trajectory
from scipy.spatial.distance import pdist

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
G2 = ROOT / 'calculations/G2_slab'

SURFACES = {
    'S1':  'S1_Pd100',
    'S2':  'S2_PdO101_Pd100',
    'S3':  'S3_PdO100',
    'S3b': 'S3b_PdO100_PdOterm',
    'S4':  'S4_PdO2_110',
}
N_ADS = 7

# Geometry validator thresholds
ADS_DIRECT_MAX = 6.0     # Å — if any adsorbate pair > 6 Å in direct distance → fragmented
ADS_DIRECT_MIN = 0.8     # Å — atoms overlap → collapsed
CO_BOND_RANGE = (0.9, 1.5)     # C=O (idx 0-1)
METHOXY_OC_RANGE = (1.1, 1.8)  # O-C of methoxy (idx 2-3)
CH_RANGE = (0.85, 1.4)         # methyl C-H (idx 3-{4,5,6})


def geometry_valid(atoms):
    """Check adsorbate (last 7 atoms) is structurally intact under direct (no MIC) distances."""
    ads_pos = atoms.positions[-N_ADS:]
    # All pairwise direct distances
    d = pdist(ads_pos)
    if d.max() > ADS_DIRECT_MAX:
        return False, 'fragmented'
    if d.min() < ADS_DIRECT_MIN:
        return False, 'collapsed'
    # Intramolecular bond sanity (direct, since the molecule should be unwrapped within cell)
    d_co = np.linalg.norm(ads_pos[1] - ads_pos[0])
    if not (CO_BOND_RANGE[0] <= d_co <= CO_BOND_RANGE[1]):
        return False, f'CO_bond_broken d={d_co:.2f}'
    d_oc = np.linalg.norm(ads_pos[3] - ads_pos[2])
    if not (METHOXY_OC_RANGE[0] <= d_oc <= METHOXY_OC_RANGE[1]):
        return False, f'methoxy_OC_broken d={d_oc:.2f}'
    for h_idx in [4, 5, 6]:
        d_ch = np.linalg.norm(ads_pos[h_idx] - ads_pos[3])
        if not (CH_RANGE[0] <= d_ch <= CH_RANGE[1]):
            return False, f'CH_broken d={d_ch:.2f}'
    return True, 'ok'


def ads_substrate_dmin_direct(atoms, n_sub):
    """Direct (non-MIC) ads-substrate min distance — for chemisorption confirm."""
    # Use MIC for ads-substrate because slab is periodic; that's intended.
    # But ads-ads is what gets corrupted by PBC, so we use direct only for ads-ads.
    d_all = atoms.get_all_distances(mic=True)
    sub = np.arange(n_sub)
    ads = np.arange(n_sub, n_sub + N_ADS)
    return float(d_all[np.ix_(ads, sub)].min())


def d_reactive_direct(atoms):
    """C_CO to O_methoxy direct distance (non-MIC)."""
    p = atoms.positions[-N_ADS:]
    return float(np.linalg.norm(p[0] - p[2]))


def loose_fingerprint(atoms, E):
    """LOOSER fingerprint than original (0.02 Å) → 0.1 Å.
    Plus E bin (10 meV) + d_reactive bin (0.2 Å) for dedup key.
    """
    p = atoms.positions[-N_ADS:]
    d_intra = []
    for i in range(N_ADS):
        for j in range(i + 1, N_ADS):
            d_intra.append(round(float(np.linalg.norm(p[i] - p[j])), 1))   # 0.1 Å bin
    return tuple(sorted(d_intra))


def process_surface(sid, sdir):
    src_dir = G3 / sdir / 'MLIP_phase2'
    out_dir = G3 / sdir / 'MLIP_phase2_filtered'
    out_dir.mkdir(exist_ok=True, parents=True)

    if not (src_dir / 'relaxed_SetA.traj').exists():
        print(f'[{sid}] missing src — skip')
        return None

    slab = read(G2 / sdir / 'CONTCAR')
    n_sub = len(slab)

    print(f'[{sid}] Reading relaxed_SetA.traj ...')
    cands = list(read(src_dir / 'relaxed_SetA.traj', index=':'))
    ranked_orig = json.load(open(src_dir / 'ranked_SetA.json'))
    # Map from idx -> rec
    ranked_by_idx = {r['idx']: r for r in ranked_orig}
    n_in = len(cands)
    print(f'[{sid}] {n_in} input structures')

    # Apply geometry validator
    survivors = []
    drop_reasons = {'fragmented': 0, 'collapsed': 0,
                    'CO_bond_broken': 0, 'methoxy_OC_broken': 0, 'CH_broken': 0,
                    'unconverged': 0}
    for i, atoms in enumerate(cands):
        # Atomic info: idx is sequential in the saved traj — but may not match `ranked_SetA.json` idx.
        # The saved traj was written in candidate-load order (idx 0, 1, 2, ...).
        orig_idx = i
        rec_old = ranked_by_idx.get(orig_idx)
        if rec_old is None:
            # shouldn't happen if traj has same indexing
            continue
        if not rec_old.get('converged', False):
            drop_reasons['unconverged'] += 1
            continue
        ok, reason = geometry_valid(atoms)
        if not ok:
            short_reason = reason.split()[0]  # take main keyword
            drop_reasons[short_reason] = drop_reasons.get(short_reason, 0) + 1
            continue
        # Build new record with direct distances
        rec = {
            'idx': orig_idx,
            'E': rec_old['E'],
            'converged': True,
            'n_steps': rec_old.get('n_steps', -1),
            'd_min': ads_substrate_dmin_direct(atoms, n_sub),
            'd_reactive': d_reactive_direct(atoms),
            'fingerprint': loose_fingerprint(atoms, rec_old['E']),
        }
        survivors.append(rec)

    n_out = len(survivors)
    print(f'[{sid}] survivors {n_out}/{n_in} ({100*n_out/n_in:.0f}% kept)')
    print(f'         drops: {drop_reasons}')

    if n_out == 0:
        print(f'[{sid}] no survivors!'); return None

    # Re-rank by E
    survivors.sort(key=lambda r: r['E'])
    E_min = survivors[0]['E']
    for r in survivors:
        r['dE_rel_meV'] = (r['E'] - E_min) * 1000.0

    json.dump(survivors, open(out_dir / 'ranked_SetA.json', 'w'), indent=1)

    # Re-dedup with LOOSER fingerprint
    seen = {}
    for r in survivors:
        key = (round(r['E'], 2), round(r['d_reactive'], 1), r['fingerprint'])
        if key not in seen:
            seen[key] = r
    unique = sorted(seen.values(), key=lambda r: r['E'])
    json.dump(unique, open(out_dir / 'unique_SetA.json', 'w'), indent=1)
    print(f'[{sid}] unique after dedup: {len(unique)} (from {n_out})')

    # Re-pick distance-bin shortlist FROM filtered+unique pool
    # Only keep in canonical Set A band [2.1, 4.0)
    bins = [(2.1, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.0)]
    shortlist = []
    for lo, hi in bins:
        in_bin = [r for r in unique if lo <= r['d_reactive'] < hi]
        if in_bin:
            shortlist.append(in_bin[0])   # lowest E
    json.dump(shortlist, open(out_dir / 'shortlist_SetA.json', 'w'), indent=1)
    print(f'[{sid}] shortlist (in SetA band 2.1-4.0): {len(shortlist)}')

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.scatter([r['d_reactive'] for r in unique],
               [r['dE_rel_meV'] for r in unique],
               s=15, alpha=0.4, color='steelblue', label=f'filtered unique (n={len(unique)})')
    if shortlist:
        ax.scatter([r['d_reactive'] for r in shortlist],
                   [r['dE_rel_meV'] for r in shortlist],
                   s=120, color='red', edgecolor='black', zorder=5,
                   label=f'shortlist (n={len(shortlist)})')
    ax.axvspan(2.1, 4.0, alpha=0.10, color='green', label='Set A band')
    ax.set_xlabel(r'$d(C_{CO} - O_{CH_3O})$  /  Å (direct)')
    ax.set_ylabel(r'$\Delta E_{MACE+D3}$  /  meV')
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'shortlist_SetA.png', dpi=130, bbox_inches='tight')
    plt.close()

    summary = {
        'surface': sid,
        'sdir': sdir,
        'n_input': n_in,
        'n_unconverged_dropped': drop_reasons['unconverged'],
        'n_fragmented_dropped': drop_reasons.get('fragmented', 0),
        'n_collapsed_dropped': drop_reasons.get('collapsed', 0),
        'n_bond_broken_dropped': (drop_reasons.get('CO_bond_broken', 0)
                                   + drop_reasons.get('methoxy_OC_broken', 0)
                                   + drop_reasons.get('CH_broken', 0)),
        'n_survivors': n_out,
        'n_unique_after_dedup': len(unique),
        'n_shortlist_setA_band': len(shortlist),
        'E_range_meV': survivors[-1]['dE_rel_meV'] if len(survivors) > 1 else 0,
        'shortlist_d_reactive': [r['d_reactive'] for r in shortlist],
        'shortlist_E_meV': [r['dE_rel_meV'] for r in shortlist],
        'drop_reasons': drop_reasons,
    }
    json.dump(summary, open(out_dir / 'summary.json', 'w'), indent=2)
    return summary


def main():
    print('=== Phase 2 re-filter with geometry validator ===')
    print()
    all_summary = []
    for sid, sdir in SURFACES.items():
        s = process_surface(sid, sdir)
        if s: all_summary.append(s)
        print()

    json.dump(all_summary, open(G3 / 'MLIP_phase2_filtered_summary.json', 'w'), indent=2)
    print('=== Summary ===')
    print(f"{'Sur':<5} {'n_in':>6} {'survive':>8} {'%':>6} {'unique':>7} {'short':>6} {'E_rng_meV':>10}")
    for s in all_summary:
        pct = 100 * s['n_survivors'] / s['n_input']
        print(f"{s['surface']:<5} {s['n_input']:>6} {s['n_survivors']:>8} {pct:>5.1f}% "
              f"{s['n_unique_after_dedup']:>7} {s['n_shortlist_setA_band']:>6} {s['E_range_meV']:>10.0f}")


if __name__ == '__main__':
    main()
