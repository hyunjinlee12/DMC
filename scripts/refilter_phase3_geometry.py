"""Refilter Phase 3 SetTS + SetB results with post-relax geometry validator.

Same fix pattern as scripts/refilter_phase2_geometry.py (PBC fragmentation +
intramolecular bond sanity + looser fingerprint), adapted for two pools.

Per guideline (docs/DMC_Pd_workplan.md T2.5 + docs/DMC_Pd_package_guideline.md):
  - SetB: thermodynamic reference, d_reactive ≥ 5.0 Å maintained (filter violators)
  - SetTS: pool no longer needed (T2.5 uses ASE interpolate / IDPP on endpoints;
    MACE cannot find saddle points). We still refilter for diagnostic value.

For each surface × pool:
  1. Read relaxed_{pool}.traj
  2. Apply geometry_valid()  -- drop fragmented / collapsed / bond-broken / unconverged
  3. SetB additional check: d_reactive_direct >= 5.0 Å (else drop)
  4. Re-rank by E
  5. Re-dedup (looser fingerprint)
  6. Save to MLIP_phase3_filtered/
"""
import json, shutil
from pathlib import Path
import numpy as np
from ase.io import read
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
POOLS = ['SetTS', 'SetB']

# Validator thresholds (same as Phase 2 refilter)
ADS_DIRECT_MAX = 6.0
ADS_DIRECT_MIN = 0.8
CO_BOND_RANGE = (0.9, 1.5)
METHOXY_OC_RANGE = (1.1, 1.8)
CH_RANGE = (0.85, 1.4)
SETB_MIN_REACTIVE = 5.0    # Å — SetB must keep reactive distance >= 5


def geometry_valid(atoms):
    ads_pos = atoms.positions[-N_ADS:]
    d = pdist(ads_pos)
    if d.max() > ADS_DIRECT_MAX: return False, 'fragmented'
    if d.min() < ADS_DIRECT_MIN: return False, 'collapsed'
    d_co = np.linalg.norm(ads_pos[1] - ads_pos[0])
    if not (CO_BOND_RANGE[0] <= d_co <= CO_BOND_RANGE[1]):
        return False, 'CO_bond_broken'
    d_oc = np.linalg.norm(ads_pos[3] - ads_pos[2])
    if not (METHOXY_OC_RANGE[0] <= d_oc <= METHOXY_OC_RANGE[1]):
        return False, 'methoxy_OC_broken'
    for h_idx in [4, 5, 6]:
        d_ch = np.linalg.norm(ads_pos[h_idx] - ads_pos[3])
        if not (CH_RANGE[0] <= d_ch <= CH_RANGE[1]):
            return False, 'CH_broken'
    return True, 'ok'


def d_reactive_direct(atoms):
    p = atoms.positions[-N_ADS:]
    return float(np.linalg.norm(p[0] - p[2]))


def ads_substrate_dmin(atoms, n_sub):
    """MIC is fine here (slab is genuinely periodic)."""
    d_all = atoms.get_all_distances(mic=True)
    sub = np.arange(n_sub)
    ads = np.arange(n_sub, n_sub + N_ADS)
    return float(d_all[np.ix_(ads, sub)].min())


def loose_fingerprint(atoms):
    p = atoms.positions[-N_ADS:]
    dists = []
    for i in range(N_ADS):
        for j in range(i + 1, N_ADS):
            dists.append(round(float(np.linalg.norm(p[i] - p[j])), 1))
    return tuple(sorted(dists))


def process_pool(sid, sdir, pool):
    src_dir = G3 / sdir / 'MLIP_phase3'
    out_dir = G3 / sdir / 'MLIP_phase3_filtered'
    out_dir.mkdir(exist_ok=True, parents=True)

    traj_path = src_dir / f'relaxed_{pool}.traj'
    ranked_path = src_dir / f'ranked_{pool}.json'
    if not traj_path.exists() or not ranked_path.exists():
        print(f'[{sid} {pool}] missing src — skip')
        return None

    slab = read(G2 / sdir / 'CONTCAR')
    n_sub = len(slab)

    cands = list(read(traj_path, index=':'))
    ranked_orig = json.load(open(ranked_path))
    ranked_by_idx = {r['idx']: r for r in ranked_orig}

    survivors = []
    drops = {'fragmented': 0, 'collapsed': 0, 'CO_bond_broken': 0,
             'methoxy_OC_broken': 0, 'CH_broken': 0, 'unconverged': 0,
             'setb_too_close': 0}

    for i, atoms in enumerate(cands):
        rec_old = ranked_by_idx.get(i)
        if rec_old is None: continue
        if not rec_old.get('converged', False):
            drops['unconverged'] += 1
            continue
        ok, reason = geometry_valid(atoms)
        if not ok:
            drops[reason] = drops.get(reason, 0) + 1
            continue
        d_react = d_reactive_direct(atoms)
        if pool == 'SetB' and d_react < SETB_MIN_REACTIVE:
            drops['setb_too_close'] += 1
            continue
        rec = {
            'idx': i,
            'E': rec_old['E'],
            'converged': True,
            'n_steps': rec_old.get('n_steps', -1),
            'd_min': ads_substrate_dmin(atoms, n_sub),
            'd_reactive': d_react,
            'fingerprint': loose_fingerprint(atoms),
        }
        survivors.append(rec)

    n_in = len(cands)
    n_out = len(survivors)
    print(f'[{sid} {pool}] {n_out}/{n_in} survive ({100*n_out/n_in:.0f}%)  drops={drops}')

    if n_out == 0:
        return None

    survivors.sort(key=lambda r: r['E'])
    E_min = survivors[0]['E']
    for r in survivors:
        r['dE_rel_meV'] = (r['E'] - E_min) * 1000.0
    json.dump(survivors, open(out_dir / f'ranked_{pool}.json', 'w'), indent=1)

    seen = {}
    for r in survivors:
        key = (round(r['E'], 2), round(r['d_reactive'], 1), r['fingerprint'])
        if key not in seen:
            seen[key] = r
    unique = sorted(seen.values(), key=lambda r: r['E'])
    json.dump(unique, open(out_dir / f'unique_{pool}.json', 'w'), indent=1)

    summary = {
        'surface': sid, 'pool': pool,
        'n_input': n_in, 'n_survivors': n_out, 'n_unique': len(unique),
        'E_range_meV': survivors[-1]['dE_rel_meV'] if survivors else 0,
        'drops': drops,
    }
    json.dump(summary, open(out_dir / f'summary_{pool}.json', 'w'), indent=2)
    return summary


def main():
    print('=== Phase 3 refilter (PBC fragmentation + SetB ≥5 Å) ===\n')
    all_s = []
    for sid, sdir in SURFACES.items():
        for pool in POOLS:
            s = process_pool(sid, sdir, pool)
            if s: all_s.append(s)
        print()

    json.dump(all_s, open(G3 / 'MLIP_phase3_filtered_summary.json', 'w'), indent=2)
    print('=== Summary ===')
    print(f"{'Sur':<5} {'Pool':<6} {'n_in':>5} {'survive':>8} {'%':>6} {'unique':>7} {'E_meV':>7}")
    for s in all_s:
        pct = 100 * s['n_survivors'] / s['n_input']
        print(f"{s['surface']:<5} {s['pool']:<6} {s['n_input']:>5} {s['n_survivors']:>8} {pct:>5.1f}% "
              f"{s['n_unique']:>7} {s['E_range_meV']:>7.0f}")


if __name__ == '__main__':
    main()
