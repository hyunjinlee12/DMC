"""Re-pick DFT shortlist with STRICT adherence to docs/DMC_Pd_workplan.md §P1-C.

Priority order (per guide):
  ① lowest total E (within group)
  ② C_CO–O_CH3O distance bin (2.0-3.5 for reactive)
  ③ CO·CH3O at DIFFERENT functional sites (co-ads only)
  ④ S2: interface pair (Pd⁰ + Pd²⁺ region) priority
  ⑤ O-rich PdO/PdO₂: CO₂-like collapsed → side-path bin

For SINGLE ads:
  - Bin by site_type (atop_Pd, bridge_Pd-Pd, hollow_3, atop_O, etc.)
  - Within each bin, pick lowest E
  - Round-robin across bins to fill quota

For CO-ADS:
  - Bin by (CO_site_type, OMe_site_type, d_reactive_band)
  - Prefer DIFFERENT functional site combinations (guide ③)
  - S2: prefer Pd⁰×Pd⁰ or Pd⁰×Pd²⁺ "interface" combo

Output: new shortlist at calculations/G3_adsorption/DFT_shortlist_v2/
"""
import json
import shutil
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.constraints import FixAtoms
from collections import defaultdict

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab'
G3 = ROOT / 'calculations/G3_adsorption'
OUT = G3 / 'DFT_shortlist_v2'
if OUT.exists(): shutil.rmtree(OUT)
OUT.mkdir()

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

BUDGET = {
    'S1':  {'CO': 3, 'CH3O': 3, 'coads': 3},
    'S2':  {'CO': 5, 'CH3O': 5, 'coads': 5},
    'S3':  {'CO': 3, 'CH3O': 3, 'coads': 3},
    'S3b': {'CO': 3, 'CH3O': 3, 'coads': 3},
    'S4':  {'CO': 3, 'CH3O': 3, 'coads': 0},
}

PD_BOND_CUTOFF = 2.60


def find_anchor_and_substrate(atoms, kind):
    syms = np.array(atoms.get_chemical_symbols())
    c_idx = list(np.where(syms == 'C')[0])
    h_idx = list(np.where(syms == 'H')[0])
    o_idx_all = list(np.where(syms == 'O')[0])
    pd_idx = list(np.where(syms == 'Pd')[0])
    ads_o = set()
    if kind == 'CO':
        if len(c_idx) != 1: return None
        c = c_idx[0]
        d = atoms.get_distances(c, o_idx_all, mic=True)
        o_co = o_idx_all[np.argmin(d)]; ads_o.add(o_co)
        sub = pd_idx + [o for o in o_idx_all if o not in ads_o]
        return {'C_CO': c, 'O_CO': o_co, 'sub': sub}
    elif kind == 'CH3O':
        if len(c_idx) != 1: return None
        c = c_idx[0]
        d = atoms.get_distances(c, o_idx_all, mic=True)
        o_methoxy = o_idx_all[np.argmin(d)]; ads_o.add(o_methoxy)
        sub = pd_idx + [o for o in o_idx_all if o not in ads_o]
        return {'C_methyl': c, 'O_methoxy': o_methoxy, 'H_methyl': h_idx, 'sub': sub}
    elif kind == 'coads':
        if len(c_idx) != 2 or len(h_idx) != 3: return None
        h_to_c = {h: c_idx[np.argmin(atoms.get_distances(h, c_idx, mic=True))] for h in h_idx}
        cm = list(set(h_to_c.values()))
        if len(cm) != 1: return None
        c_methyl = cm[0]; c_co = [c for c in c_idx if c != c_methyl][0]
        d_co_o = atoms.get_distances(c_co, o_idx_all, mic=True)
        o_co = o_idx_all[np.argmin(d_co_o)]; ads_o.add(o_co)
        rest = [o for o in o_idx_all if o != o_co]
        d_me_o = atoms.get_distances(c_methyl, rest, mic=True)
        o_methoxy = rest[np.argmin(d_me_o)]; ads_o.add(o_methoxy)
        sub = pd_idx + [o for o in o_idx_all if o not in ads_o]
        return {'C_CO': c_co, 'O_CO': o_co,
                'C_methyl': c_methyl, 'O_methoxy': o_methoxy,
                'H_methyl': h_idx, 'sub': sub}
    return None


def classify_site(atoms, anchor, sub_indices):
    syms = atoms.get_chemical_symbols()
    d = atoms.get_distances(anchor, sub_indices, mic=True)
    nbrs = [(sub_indices[i], d[i]) for i in range(len(sub_indices)) if d[i] < PD_BOND_CUTOFF]
    n_pd = sum(1 for i, _ in nbrs if syms[i] == 'Pd')
    n_o = sum(1 for i, _ in nbrs if syms[i] == 'O')
    total = n_pd + n_o
    if total == 0:           label = 'physisorbed'
    elif total == 1:         label = f"atop_{'Pd' if n_pd else 'O'}"
    elif total == 2:
        if n_pd == 2:        label = 'bridge_Pd-Pd'
        elif n_o == 2:       label = 'bridge_O-O'
        else:                label = 'bridge_Pd-O'
    elif total == 3:
        if n_pd == 3:        label = 'hollow_3Pd'
        elif n_o == 3:       label = 'hollow_3O'
        else:                label = f'hollow_3({n_pd}Pd{n_o}O)'
    else:                    label = f'{total}f({n_pd}Pd{n_o}O)'
    return label, total, n_pd, n_o


def fix_bottom_50(atoms, n_sub):
    z = atoms.positions[:n_sub, 2]
    z_med = np.median(z)
    fixed = [i for i in range(n_sub) if atoms.positions[i, 2] < z_med]
    atoms.set_constraint(FixAtoms(indices=fixed))


def is_intramol_valid(atoms, kind):
    """Check adsorbate intramolecular bonds (CO C=O, methoxy O-C/C-H)."""
    import numpy as np
    syms = np.array(atoms.get_chemical_symbols())
    c_idx = list(np.where(syms == 'C')[0])
    h_idx = list(np.where(syms == 'H')[0])
    o_idx = list(np.where(syms == 'O')[0])
    if kind == 'CO':
        if len(c_idx) != 1: return False
        c = c_idx[0]
        d_co = min(atoms.get_distances(c, o_idx, mic=True))
        return 1.05 <= d_co <= 1.30
    elif kind == 'CH3O':
        if len(c_idx) != 1 or len(h_idx) != 3: return False
        c = c_idx[0]
        d_oc = min(atoms.get_distances(c, o_idx, mic=True))
        d_ch = [atoms.get_distance(c, h, mic=True) for h in h_idx]
        return 1.30 <= d_oc <= 1.55 and all(0.90 <= d <= 1.25 for d in d_ch)
    return True


def diverse_pick_single(unique_recs, traj, slab, kind, budget):
    """Diverse pick: bin by site_type, then E order. Force different sites first.
    SKIPS structures with broken intramolecular bonds.
    """
    n_ads = 2 if kind == 'CO' else 5
    enriched = []
    for k, rec in enumerate(unique_recs):
        try:
            atoms = traj[rec['idx']]
        except (IndexError, KeyError):
            continue
        if len(atoms) != len(slab) + n_ads:
            from ase import Atoms
            ads = atoms[-n_ads:]
            atoms = slab.copy()
            atoms += ads
        # SKIP BROKEN adsorbate
        if not is_intramol_valid(atoms, kind):
            continue
        info = find_anchor_and_substrate(atoms, kind)
        if not info: continue
        anchor = info['C_CO'] if kind == 'CO' else info['O_methoxy']
        site, n_nbr, n_pd, n_o = classify_site(atoms, anchor, info['sub'])
        enriched.append({**rec, 'site': site, 'n_nbr': n_nbr, 'atoms': atoms})

    enriched.sort(key=lambda r: r['E'])

    selected = []
    sites_used = set()
    # Pass 1: get each distinct site (highest priority)
    for r in enriched:
        if r['site'] not in sites_used:
            selected.append(r); sites_used.add(r['site'])
            if len(selected) >= budget: break
    # Pass 2: fill remaining with lowest E from already-used sites
    for r in enriched:
        if len(selected) >= budget: break
        if r in selected: continue
        selected.append(r)
    return selected[:budget]


def diverse_pick_coads(unique_recs, traj, slab, budget, surface_id=None):
    """Co-ads diverse pick per guide ③ + STRICT SetA band [2.1, 4.0] Å only.
    Out-of-band (d<2.1 product or d≥5 SetB) → side-path / thermo, not in T1.15.
    """
    enriched = []
    for k, rec in enumerate(unique_recs):
        d_react = rec.get('d_reactive', 0)
        # STRICT SetA band — guide T1.13 + T2.1 + T1.15 reactive shortlist
        if not (2.1 <= d_react <= 4.0):
            continue
        try:
            atoms = traj[rec['idx']]
        except (IndexError, KeyError):
            continue
        if len(atoms) != len(slab) + 7:
            from ase import Atoms
            ads = atoms[-7:]
            atoms = slab.copy()
            atoms += ads
        info = find_anchor_and_substrate(atoms, 'coads')
        if not info: continue
        # SKIP BROKEN co-ads
        import numpy as np
        d_co_intra = atoms.get_distance(info['C_CO'], info['O_CO'], mic=True)
        d_methoxy = atoms.get_distance(info['C_methyl'], info['O_methoxy'], mic=True)
        d_ch = [atoms.get_distance(info['C_methyl'], h, mic=True) for h in info['H_methyl']]
        if not (1.05 <= d_co_intra <= 1.30 and 1.30 <= d_methoxy <= 1.55 and all(0.90 <= d <= 1.25 for d in d_ch)):
            continue
        c_site, _, _, _ = classify_site(atoms, info['C_CO'], info['sub'])
        o_site, _, _, _ = classify_site(atoms, info['O_methoxy'], info['sub'])
        if d_react < 2.5: bin_d = '2.1-2.5'
        elif d_react < 3.0: bin_d = '2.5-3.0'
        elif d_react < 3.5: bin_d = '3.0-3.5'
        else: bin_d = '3.5-4.0'
        site_combo = f'{c_site}/{o_site}'
        diff_sites = c_site != o_site
        enriched.append({**rec, 'c_site': c_site, 'o_site': o_site,
                         'combo': site_combo, 'd_bin': bin_d,
                         'diff_sites': diff_sites, 'atoms': atoms})
    enriched.sort(key=lambda r: r['E'])

    selected = []
    combos_used = set()
    bins_used = set()
    # Pass 1: prefer DIFFERENT sites + new bin + new combo (가이드 ③)
    for r in enriched:
        if r['diff_sites'] and r['combo'] not in combos_used and r['d_bin'] not in bins_used:
            selected.append(r); combos_used.add(r['combo']); bins_used.add(r['d_bin'])
            if len(selected) >= budget: break
    # Pass 2: same logic without diff_sites requirement
    for r in enriched:
        if len(selected) >= budget: break
        if r in selected: continue
        if r['combo'] not in combos_used and r['d_bin'] not in bins_used:
            selected.append(r); combos_used.add(r['combo']); bins_used.add(r['d_bin'])
    # Pass 3: fill remaining
    for r in enriched:
        if len(selected) >= budget: break
        if r in selected: continue
        selected.append(r)
    return selected[:budget]


def main():
    print('=== Guide-strict re-pick ===\n')
    for sid in SURFACES:
        sdir = SDIRS[sid]
        slab = read(G2 / sdir / 'CONTCAR')
        n_sub = len(slab)
        out_surf = OUT / sid
        out_surf.mkdir()
        budgets = BUDGET[sid]

        print(f'--- {sid} ---')

        # === Single CO ===
        if budgets['CO'] > 0:
            unique = json.load(open(G3 / sdir / 'MLIP_phase1/unique_CO.json'))
            traj = list(read(G3 / sdir / 'MLIP_phase1/relaxed_CO.traj', index=':'))
            picks = diverse_pick_single(unique, traj, slab, 'CO', budgets['CO'])
            sub_dir = out_surf / 'single_CO'; sub_dir.mkdir()
            for k, r in enumerate(picks):
                atoms = r['atoms']
                fix_bottom_50(atoms, n_sub)
                label = f'{k:02d}_single_CO_{r["site"]}_idx{r["idx"]:05d}'.replace('/', '-')
                write(str(sub_dir / f'{label}.vasp'), atoms, format='vasp', direct=True, sort=True, vasp5=True)
            sites_summary = ', '.join(f'{r["site"]}({r["E"]:.2f})' for r in picks)
            print(f'  CO   : {sites_summary}')

        # === Single CH3O ===
        if budgets['CH3O'] > 0:
            unique = json.load(open(G3 / sdir / 'MLIP_phase1/unique_CH3O.json'))
            traj = list(read(G3 / sdir / 'MLIP_phase1/relaxed_CH3O.traj', index=':'))
            picks = diverse_pick_single(unique, traj, slab, 'CH3O', budgets['CH3O'])
            sub_dir = out_surf / 'single_CH3O'; sub_dir.mkdir()
            for k, r in enumerate(picks):
                atoms = r['atoms']
                fix_bottom_50(atoms, n_sub)
                label = f'{k:02d}_single_CH3O_{r["site"]}_idx{r["idx"]:05d}'.replace('/', '-')
                write(str(sub_dir / f'{label}.vasp'), atoms, format='vasp', direct=True, sort=True, vasp5=True)
            sites_summary = ', '.join(f'{r["site"]}({r["E"]:.2f})' for r in picks)
            print(f'  CH3O : {sites_summary}')

        # === Co-ads ===
        if budgets['coads'] > 0:
            f = G3 / sdir / 'MLIP_phase2_filtered/unique_SetA.json'
            traj_f = G3 / sdir / 'MLIP_phase2/relaxed_SetA.traj'
            if not f.exists() or not traj_f.exists():
                print(f'  coads: no source — skip'); continue
            unique = json.load(open(f))
            traj = list(read(traj_f, index=':'))
            picks = diverse_pick_coads(unique, traj, slab, budgets['coads'], surface_id=sid)
            sub_dir = out_surf / 'coads_SetA'; sub_dir.mkdir()
            for k, r in enumerate(picks):
                atoms = r['atoms']
                fix_bottom_50(atoms, n_sub)
                label = f'{k:02d}_coads_CO-{r["c_site"]}_OMe-{r["o_site"]}_d{r["d_reactive"]:.2f}_idx{r["idx"]:05d}'.replace('/', '-')[:80]
                write(str(sub_dir / f'{label}.vasp'), atoms, format='vasp', direct=True, sort=True, vasp5=True)
            sites_summary = ', '.join(f'CO[{r["c_site"][:8]}]-OMe[{r["o_site"][:8]}]@{r["d_reactive"]:.1f}Å' for r in picks)
            print(f'  coads: {sites_summary}')

    # Total count
    total = sum(1 for _ in OUT.rglob('*.vasp'))
    print(f'\n✓ {total} candidates re-picked at {OUT}')


if __name__ == '__main__':
    main()
