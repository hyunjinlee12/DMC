"""Render representative initial heuristic structures (CO*, CH3O*, co-ads)
to verify adsorbate floating height before DFT/MLIP relaxation.

NOTE: SetA.traj stores substrate trimmed (only bottom subsurface layer).
Adsorbate xyz is absolute, so we MUST reattach to the full G2 slab to get
a faithful render. CO/CH3O candidates.traj also have the same trimming,
so do the same reattach for them.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read
from ase import Atoms
from ase.visualize.plot import plot_atoms

plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12})

ROOT_G3 = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G3_adsorption')
ROOT_G2 = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab')
OUT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/floating_check')
OUT.mkdir(parents=True, exist_ok=True)

surfaces = {
    'S1':  ('S1_Pd100',            'Pd(100)'),
    'S2':  ('S2_PdO101_Pd100',     'PdO(101)/Pd(100)'),
    'S3':  ('S3_PdO100',           'PdO(100) O-term'),
    'S3b': ('S3b_PdO100_PdOterm',  'PdO(100) PdO-term'),
    'S4':  ('S4_PdO2_110',         'PdO2(110)'),
}

def reattach(cand_atoms, slab, n_ads):
    """Take last n_ads atoms of cand_atoms (adsorbate, absolute xyz)
    and concatenate to the full G2 slab."""
    ads = cand_atoms[-n_ads:]
    full = slab + ads
    full.set_pbc(slab.get_pbc())
    full.set_cell(slab.cell)
    return full

def ads_min_dist(atoms, n_ads):
    n = len(atoms)
    sub = np.arange(n - n_ads)
    ads = np.arange(n - n_ads, n)
    return atoms.get_all_distances(mic=True)[np.ix_(ads, sub)].min()

def ads_dz(atoms, n_ads):
    return atoms.positions[-n_ads:, 2].min() - atoms.positions[:-n_ads, 2].max()

def render_one(atoms, title, outpath, n_ads):
    dmin = ads_min_dist(atoms, n_ads)
    dz = ads_dz(atoms, n_ads)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    plot_atoms(atoms, axes[0], rotation='0x,0y,0z', radii=0.85, show_unit_cell=2)
    axes[0].set_title(f'{title}  —  top')
    axes[0].set_xticks([]); axes[0].set_yticks([])
    plot_atoms(atoms, axes[1], rotation='-90x', radii=0.85, show_unit_cell=2)
    axes[1].set_title(f'side  —  d(min ads–sub) = {dmin:.2f} Å,  Δz_top = {dz:.2f} Å')
    axes[1].set_xticks([]); axes[1].set_yticks([])
    plt.tight_layout()
    plt.savefig(outpath, dpi=130, bbox_inches='tight')
    plt.close()
    return dmin, dz

print(f"{'Surface':<5} {'Ads':<8} {'d_min':>7} {'Δz_top':>7}")
print('-' * 35)
rows = []

for sid, (sdir, label) in surfaces.items():
    base = ROOT_G3 / sdir
    slab = read(ROOT_G2 / sdir / 'CONTCAR')

    # CO* (2 atoms)
    co = list(read(base / 'CO/candidates.traj', index=':'))
    a = reattach(co[len(co)//2], slab, 2)
    d, dz = render_one(a, f'{sid} {label} — CO*', OUT / f'{sid}_CO.png', 2)
    print(f"{sid:<5} {'CO*':<8} {d:7.2f} {dz:7.2f}")
    rows.append((sid, 'CO*', d, dz))

    # CH3O* (5 atoms)
    ch = list(read(base / 'CH3O/candidates.traj', index=':'))
    a = reattach(ch[len(ch)//2], slab, 5)
    d, dz = render_one(a, f'{sid} {label} — CH$_3$O*', OUT / f'{sid}_CH3O.png', 5)
    print(f"{sid:<5} {'CH3O*':<8} {d:7.2f} {dz:7.2f}")
    rows.append((sid, 'CH3O*', d, dz))

    # co-ads SetA (7 atoms: CO + CH3O)
    coa = list(read(base / 'coads_guide/SetA.traj', index=':'))
    a = reattach(coa[len(coa)//2], slab, 7)
    d, dz = render_one(a, f'{sid} {label} — CO* + CH$_3$O* (Set A reactive)', OUT / f'{sid}_coads.png', 7)
    print(f"{sid:<5} {'co-ads':<8} {d:7.2f} {dz:7.2f}")
    rows.append((sid, 'co-ads', d, dz))

print()
print('Saved 15 figures →', OUT)
print()
print(f"{'Sur':<5} {'Ads':<8} {'d_min (Å)':>10} {'Δz_top (Å)':>12}")
for sid, ad, d, dz in rows:
    print(f"{sid:<5} {ad:<8} {d:10.2f} {dz:12.2f}")
