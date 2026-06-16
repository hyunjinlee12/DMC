#!/usr/bin/env python
"""Quick test of CO* population with fixes"""

import sys
sys.stdout = sys.stderr = open('test_co.log', 'w', buffering=1)

from ase.io import read, write
from autoadsorbate import Surface, Fragment

print("Loading slab...")
slab = read("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab/S1_Pd100/CONTCAR")
print(f"Loaded: {len(slab)} Pd atoms")

print("Creating Surface...")
surface = Surface(slab, precision=0.25, mode='slab')
print(f"Surface created, {len(surface.site_df)} sites")

print("Creating CO* Fragment...")
co_smiles = 'ClC=O'
co_frag = Fragment(co_smiles, to_initialize=1, random_seed=2104)
print(f"Fragment created, {len(co_frag.conformers)} conformers")

# Apply fix
for conf in co_frag.conformers:
    if 'smiles' not in conf.info:
        conf.info['smiles'] = co_smiles
print("Applied smiles fix to conformers")

print("Populating sites (this will take several minutes)...")
co_structs = surface.get_populated_sites(
    co_frag,
    mode='all',
    sample_rotation=False,
    conformers_per_site_cap=1,
    overlap_thr=1.25
)

print(f"SUCCESS! Generated {len(co_structs)} CO* structures")
write('test_co_out.traj', co_structs)
print("Saved to test_co_out.traj")
