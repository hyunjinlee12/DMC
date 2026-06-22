"""Manually construct S4 CH3O* replacements for the 2 BROKEN candidates.

Strategy: place intact CH3O at known high-symmetry sites on S4 PdO2(110):
  R1: atop top-layer Pd
  R2: atop bridging O (lattice O site, not surface bonded — methoxy O coordinates to Pd nearby)

CH3O geometry (intact reference):
  O-C bond: 1.42 Å
  C-H bonds: 1.10 Å (3 H)
  H-C-H angle: ~110° (tetrahedral)
  O-C tilt: O down toward surface

Replacements numbered 90/91 to keep separate from rank 0/1/2 of MLIP.
"""
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
from pathlib import Path

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab'
OUT = ROOT / 'calculations/G3_adsorption/DFT_shortlist/S4/single_CH3O'

def build_ch3o(o_pos, axis_up=np.array([0, 0, 1])):
    """Build CH3O fragment with O at o_pos, C above, 3 H rotated symmetrically.
    Returns positions for [O, C, H, H, H] (5 atoms).
    """
    # O-C bond 1.42 along axis (up direction)
    c_pos = o_pos + 1.42 * axis_up
    # 3 H at tetrahedral positions around C (away from O)
    # H-C-O angle ~110°, H-C-H 109°
    # Place 3 H at 120° rotations around the O-C axis, tilted by ~70° from up
    tilt_deg = 70
    tilt_rad = np.deg2rad(tilt_deg)
    # base unit vectors perp to up
    perp1 = np.array([1, 0, 0])
    perp2 = np.array([0, 1, 0])
    h_positions = []
    for theta_deg in [0, 120, 240]:
        theta = np.deg2rad(theta_deg)
        # H direction: tilt by tilt_deg from +axis, rotated by theta around axis
        h_dir = np.cos(tilt_rad) * axis_up + np.sin(tilt_rad) * (np.cos(theta) * perp1 + np.sin(theta) * perp2)
        h_pos = c_pos + 1.10 * h_dir
        h_positions.append(h_pos)
    return np.array([o_pos, c_pos] + h_positions)


def fix_bottom_50(atoms):
    syms = atoms.get_chemical_symbols()
    n_sub = len(atoms) - 5   # CH3O = 5
    z = atoms.positions[:n_sub, 2]
    z_med = np.median(z)
    fixed = [i for i in range(n_sub) if atoms.positions[i, 2] < z_med]
    atoms.set_constraint(FixAtoms(indices=fixed))


# Load S4 clean slab
slab = read(G2 / 'S4_PdO2_110/CONTCAR')
syms = np.array(slab.get_chemical_symbols())
z = slab.positions[:, 2]
z_top = z.max()
print(f'S4 slab: {len(slab)} atoms, z_top = {z_top:.2f}')

# Find top-layer Pd atoms (PdO2(110): top is O, Pd just below)
pd_idx_all = np.where(syms == 'Pd')[0]
# Widen window: top-layer Pd within 1.5 Å of topmost atom
pd_top = [i for i in pd_idx_all if z[i] > z_top - 1.5]
# Sort by z desc
pd_top = sorted(pd_top, key=lambda i: -z[i])
print(f'Top-layer Pd atoms (within 1.5 Å of top): {len(pd_top)}, z range {z[pd_top].min():.2f}-{z[pd_top].max():.2f}')

# Find top-layer O atoms
o_idx_all = np.where(syms == 'O')[0]
o_top = [i for i in o_idx_all if z[i] > z_top - 0.5]
print(f'Top-layer O atoms (top 0.5 Å): {len(o_top)}')

# ---------- R1: CH3O atop a top-layer Pd ----------
pd_anchor = pd_top[0]
pd_pos = slab.positions[pd_anchor]
# Place methoxy O at Pd_top.z + 2.05 Å (Pd-O chemisorbed)
o_pos_R1 = np.array([pd_pos[0], pd_pos[1], pd_pos[2] + 2.05])
ch3o_R1 = build_ch3o(o_pos_R1)
atoms_R1 = slab.copy()
atoms_R1 += Atoms('OCHHH', positions=ch3o_R1)
fix_bottom_50(atoms_R1)
print(f'R1 (atop Pd): O placed at {o_pos_R1}')

# ---------- R2: CH3O atop bridge between two top Pd ----------
pd_a, pd_b = pd_top[0], pd_top[1] if len(pd_top) > 1 else pd_top[0]
mid = (slab.positions[pd_a] + slab.positions[pd_b]) / 2.0
o_pos_R2 = np.array([mid[0], mid[1], max(slab.positions[pd_a, 2], slab.positions[pd_b, 2]) + 2.00])
ch3o_R2 = build_ch3o(o_pos_R2)
atoms_R2 = slab.copy()
atoms_R2 += Atoms('OCHHH', positions=ch3o_R2)
fix_bottom_50(atoms_R2)
print(f'R2 (bridge Pd-Pd): O placed at {o_pos_R2}')

# Save POSCARs
out_R1 = OUT / '90_single_CH3O_manual_atop_Pd.vasp'
out_R2 = OUT / '91_single_CH3O_manual_bridge.vasp'
write(str(out_R1), atoms_R1, format='vasp', direct=True, sort=True, vasp5=True)
write(str(out_R2), atoms_R2, format='vasp', direct=True, sort=True, vasp5=True)
print(f'\nSaved:')
print(f'  {out_R1}')
print(f'  {out_R2}')

# Quick sanity check
for fpath in [out_R1, out_R2]:
    a = read(fpath)
    n_C = a.get_chemical_symbols().count('C')
    n_H = a.get_chemical_symbols().count('H')
    n_O = a.get_chemical_symbols().count('O')
    n_Pd = a.get_chemical_symbols().count('Pd')
    # Find methoxy
    syms_a = np.array(a.get_chemical_symbols())
    c_idx = np.where(syms_a == 'C')[0][0]
    h_indices = np.where(syms_a == 'H')[0]
    d_co = min(a.get_distances(c_idx, np.where(syms_a == 'O')[0], mic=True))
    d_ch = [a.get_distance(c_idx, h, mic=True) for h in h_indices]
    print(f'  {fpath.name}: n={len(a)} (Pd:{n_Pd} O:{n_O} C:{n_C} H:{n_H})  d(O-C)={d_co:.2f}  d(C-H)={[f"{d:.2f}" for d in d_ch]}')

# Mark broken ones for archive (rename to .bak)
broken_files = ['00_single_CH3O_rank0_idx00065.vasp', '02_single_CH3O_rank2_idx00112.vasp']
for bf in broken_files:
    src = OUT / bf
    if src.exists():
        dst = OUT / (bf + '.broken_bak')
        src.rename(dst)
        print(f'Archived: {src.name} -> {dst.name}')
