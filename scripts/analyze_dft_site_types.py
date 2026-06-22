"""Analyze adsorption site type for each of 47 DFT shortlist candidates.

Fixed: POSCAR sort=True changes atom order. Use species + adsorbate-identification
(not last-N) to correctly distinguish substrate from adsorbate.

For each POSCAR:
  - Identify adsorbate atoms by species + bond connectivity
  - All Pd atoms = substrate Pd (Pd_pv)
  - All O atoms minus methoxy O / CO O = substrate O
  - All C atoms = adsorbate (substrate has no C)
  - All H atoms = adsorbate (substrate has no H)
  - Find substrate atoms within 2.6 Å of anchor → site classification
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read
from collections import Counter

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
SHORTLIST = ROOT / 'calculations/G3_adsorption/DFT_shortlist'
OUT = ROOT / 'reports/G3/DFT_certification'
OUT.mkdir(parents=True, exist_ok=True)

PD_BOND_CUTOFF = 2.60   # Å
SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']


def find_substrate_and_adsorbate(atoms, kind):
    """Use species + bond connectivity to identify substrate vs adsorbate atoms.
    Returns:
      sub_indices: list of substrate atom indices (Pd + lattice O)
      ads_indices: dict of adsorbate atoms by role
    """
    syms = np.array(atoms.get_chemical_symbols())
    c_idx = list(np.where(syms == 'C')[0])
    h_idx = list(np.where(syms == 'H')[0])
    o_idx_all = list(np.where(syms == 'O')[0])
    pd_idx = list(np.where(syms == 'Pd')[0])

    ads_indices = {}
    ads_o_set = set()   # adsorbate O atoms (NOT lattice)

    if kind == 'single_CO':
        if len(c_idx) != 1: return None, None
        c = c_idx[0]
        dists = atoms.get_distances(c, o_idx_all, mic=True)
        o_co = o_idx_all[np.argmin(dists)]
        ads_o_set.add(o_co)
        ads_indices = {'C_CO': c, 'O_CO': o_co}

    elif kind == 'single_CH3O':
        if len(c_idx) != 1: return None, None
        c = c_idx[0]
        dists = atoms.get_distances(c, o_idx_all, mic=True)
        o_methoxy = o_idx_all[np.argmin(dists)]
        ads_o_set.add(o_methoxy)
        ads_indices = {'C_methyl': c, 'O_methoxy': o_methoxy, 'H_methyl': h_idx}

    elif kind == 'coads_SetA':
        if len(c_idx) != 2 or len(h_idx) != 3: return None, None
        h_to_c = {h: c_idx[np.argmin(atoms.get_distances(h, c_idx, mic=True))] for h in h_idx}
        c_methyl = list(set(h_to_c.values()))
        if len(c_methyl) != 1: return None, None
        c_methyl = c_methyl[0]
        c_co = [c for c in c_idx if c != c_methyl][0]
        d_co_o = atoms.get_distances(c_co, o_idx_all, mic=True)
        o_co = o_idx_all[np.argmin(d_co_o)]
        rest = [o for o in o_idx_all if o != o_co]
        d_me_o = atoms.get_distances(c_methyl, rest, mic=True)
        o_methoxy = rest[np.argmin(d_me_o)]
        ads_o_set.update([o_co, o_methoxy])
        ads_indices = {'C_CO': c_co, 'O_CO': o_co,
                       'C_methyl': c_methyl, 'O_methoxy': o_methoxy,
                       'H_methyl': h_idx}

    # substrate = all Pd + (all O minus adsorbate O)
    sub_o = [o for o in o_idx_all if o not in ads_o_set]
    sub_indices = pd_idx + sub_o
    return sub_indices, ads_indices


def classify_anchor(atoms, anchor_idx, sub_indices):
    """Classify site based on substrate neighbors within PD_BOND_CUTOFF."""
    syms = atoms.get_chemical_symbols()
    pos = atoms.positions
    dists = atoms.get_distances(anchor_idx, sub_indices, mic=True)
    nbrs = [(sub_indices[k], dists[k]) for k in range(len(sub_indices)) if dists[k] < PD_BOND_CUTOFF]
    nbrs.sort(key=lambda x: x[1])
    # z above substrate top (compute regardless)
    sub_z = [pos[i, 2] for i in sub_indices]
    z_top = max(sub_z) if sub_z else 0.0
    if not nbrs:
        return {'site': 'physisorbed (>2.6 Å)', 'n': 0, 'n_Pd': 0, 'n_O': 0,
                'nearest_d': None, 'nearest_species': None,
                'z_above': float(pos[anchor_idx, 2] - z_top)}
    n_pd = sum(1 for i, _ in nbrs if syms[i] == 'Pd')
    n_o = sum(1 for i, _ in nbrs if syms[i] == 'O')
    total = n_pd + n_o
    if total == 1:
        spec = 'Pd' if n_pd == 1 else 'O'
        site = f'atop_{spec}'
    elif total == 2:
        if n_pd == 2: site = 'bridge_Pd-Pd'
        elif n_o == 2: site = 'bridge_O-O'
        else:          site = 'bridge_Pd-O'
    elif total == 3:
        if n_pd == 3: site = 'hollow_3Pd'
        elif n_o == 3: site = 'hollow_3O'
        else:          site = f'hollow_3(Pd{n_pd}O{n_o})'
    elif total == 4:
        site = f'4f_(Pd{n_pd}O{n_o})'
    else:
        site = f'over_{total}(Pd{n_pd}O{n_o})'

    # z above substrate top
    sub_z = [pos[i, 2] for i in sub_indices]
    z_top = max(sub_z)

    return {'site': site, 'n': total, 'n_Pd': n_pd, 'n_O': n_o,
            'nearest_d': float(nbrs[0][1]),
            'nearest_species': syms[nbrs[0][0]],
            'z_above': float(pos[anchor_idx, 2] - z_top)}


def analyze_one(vasp_path, kind, sid):
    atoms = read(vasp_path)
    sub_indices, ads = find_substrate_and_adsorbate(atoms, kind)
    if not ads:
        return {'verdict': 'unrecognized', 'kind': kind, 'surface': sid,
                'path': str(vasp_path.relative_to(ROOT))}
    out = {'kind': kind, 'surface': sid, 'path': str(vasp_path.relative_to(ROOT))}

    if kind == 'single_CO':
        out['C_anchor'] = classify_anchor(atoms, ads['C_CO'], sub_indices)
    elif kind == 'single_CH3O':
        out['O_anchor'] = classify_anchor(atoms, ads['O_methoxy'], sub_indices)
    elif kind == 'coads_SetA':
        out['C_CO_anchor'] = classify_anchor(atoms, ads['C_CO'], sub_indices)
        out['O_methoxy_anchor'] = classify_anchor(atoms, ads['O_methoxy'], sub_indices)
    return out


# Process all 47
records = []
for sid in SURFACES:
    for kind in ['single_CO', 'single_CH3O', 'coads_SetA']:
        d = SHORTLIST / sid / kind
        if not d.exists(): continue
        for vasp in sorted(d.glob('*.vasp')):
            if '.broken_bak' in vasp.name: continue
            records.append(analyze_one(vasp, kind, sid))

# Print table
print(f"{'Sur':<5} {'kind':<14} {'rank':<5} {'site (anchor)':<28} {'n_nbr':<6} {'d_nrst':<7} {'spec':<5} {'z':<7}")
print('=' * 90)
for r in records:
    rank = Path(r['path']).stem.split('_')[0]
    if r['kind'] == 'single_CO':
        s = r['C_anchor']
        site = s['site']
        n = f"{s['n']}(Pd{s['n_Pd']}/O{s['n_O']})"
        d = f"{s['nearest_d']:.2f}" if s['nearest_d'] else '—'
        spec = s['nearest_species'] or '—'
        z = f"{s['z_above']:+.2f}"
    elif r['kind'] == 'single_CH3O':
        s = r['O_anchor']
        site = s['site']
        n = f"{s['n']}(Pd{s['n_Pd']}/O{s['n_O']})"
        d = f"{s['nearest_d']:.2f}" if s['nearest_d'] else '—'
        spec = s['nearest_species'] or '—'
        z = f"{s['z_above']:+.2f}"
    else:
        sc = r['C_CO_anchor']; so = r['O_methoxy_anchor']
        site = f"CO:{sc['site'][:12]}/OMe:{so['site'][:12]}"
        n = f"{sc['n']}/{so['n']}"
        d = f"{sc['nearest_d']:.2f}/{so['nearest_d']:.2f}" if sc['nearest_d'] and so['nearest_d'] else '—'
        spec = f"{sc['nearest_species']}/{so['nearest_species']}"
        z = f"{sc['z_above']:+.2f}"
    print(f"{r['surface']:<5} {r['kind']:<14} {rank:<5} {site:<28} {n:<6} {d:<7} {spec:<5} {z:<7}")

json.dump(records, open(OUT / 'site_types.json', 'w'), indent=2)

# Per-surface distribution
print()
print('='*70)
print('Site type 분포 per surface')
print('='*70)
for sid in SURFACES:
    print(f'\n{sid}:')
    for kind in ['single_CO', 'single_CH3O', 'coads_SetA']:
        rs = [r for r in records if r['surface'] == sid and r['kind'] == kind]
        if not rs: continue
        if kind == 'single_CO':
            sites = [r['C_anchor']['site'] for r in rs]
        elif kind == 'single_CH3O':
            sites = [r['O_anchor']['site'] for r in rs]
        else:
            sites = []
            for r in rs:
                sites.append(f"CO:{r['C_CO_anchor']['site']}")
                sites.append(f"OMe:{r['O_methoxy_anchor']['site']}")
        c = Counter(sites)
        print(f"  {kind:<14}  {', '.join(f'{k}×{v}' for k, v in c.most_common())}")

# Build figure
print()
fig, axes = plt.subplots(1, 3, figsize=(17, 6))
titles = ['(a) CO* anchor (C atom)', '(b) CH₃O* anchor (methoxy O)', '(c) co-ads (both anchors counted)']
kinds = ['single_CO', 'single_CH3O', 'coads_SetA']

# Collect all site labels
all_sites = set()
for r in records:
    if r['kind'] == 'single_CO':
        all_sites.add(r['C_anchor']['site'])
    elif r['kind'] == 'single_CH3O':
        all_sites.add(r['O_anchor']['site'])
    elif r['kind'] == 'coads_SetA':
        all_sites.add(r['C_CO_anchor']['site'])
        all_sites.add(r['O_methoxy_anchor']['site'])
all_sites = sorted(all_sites)
cmap = plt.cm.tab20
site_color = {s: cmap(i / max(len(all_sites), 1)) for i, s in enumerate(all_sites)}

for ax, kind, title in zip(axes, kinds, titles):
    surf_data = {}
    for sid in SURFACES:
        c = Counter()
        for r in records:
            if r['surface'] != sid or r['kind'] != kind: continue
            if kind == 'single_CO':
                c[r['C_anchor']['site']] += 1
            elif kind == 'single_CH3O':
                c[r['O_anchor']['site']] += 1
            else:
                c[r['C_CO_anchor']['site']] += 1
                c[r['O_methoxy_anchor']['site']] += 1
        surf_data[sid] = c

    x = np.arange(len(SURFACES))
    bottoms = np.zeros(len(SURFACES))
    for site in all_sites:
        vals = [surf_data[sid].get(site, 0) for sid in SURFACES]
        if sum(vals) == 0: continue
        ax.bar(x, vals, bottom=bottoms, color=site_color[site], edgecolor='black', label=site, linewidth=0.5)
        bottoms += vals
    ax.set_xticks(x); ax.set_xticklabels(SURFACES)
    ax.set_ylabel('Count')
    ax.set_title(title)
    ax.legend(fontsize=7, loc='upper right', ncol=2)
    ax.grid(True, alpha=0.3, axis='y')
plt.suptitle('DFT Shortlist Adsorption Site Distribution (47 candidates)', fontsize=14, weight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUT / 'A_site_distribution.png', dpi=140, bbox_inches='tight')
plt.close()
print(f'Figure saved: {OUT / "A_site_distribution.png"}')
