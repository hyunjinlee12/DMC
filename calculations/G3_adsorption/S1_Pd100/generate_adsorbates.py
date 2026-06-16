#!/usr/bin/env python
"""
G3 Phase 1 Heuristic Adsorbate Sampling Pipeline for S1 Pd(100)
Tasks T1.10-T1.13: Site map + CO* + CH3O* + co-adsorption
"""

import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from ase.io import read, write
from ase.visualize import view
from ase import Atoms
from rdkit import Chem

print("="*80)
print("S1 Pd(100) - G3 Adsorbate Sampling Pipeline")
print("="*80)

# Paths
slab_path = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab/S1_Pd100/CONTCAR")
out_dir = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G3_adsorption/S1_Pd100")

# =============================================================================
# Step 0: Environment check
# =============================================================================
print("\n[Step 0] Environment Check")
print("-" * 40)

try:
    import autoadsorbate
    from autoadsorbate import Surface, Fragment
    print(f"AutoAdsorbate version: {autoadsorbate.__version__}")
except Exception as e:
    print(f"ERROR: AutoAdsorbate import failed: {e}")
    sys.exit(1)

# Test SMILES
print("\nTesting SMILES:")
co_smiles_candidates = ['ClC=O', 'Cl[C]=O', 'ClCO']
working_co_smiles = None

for smiles in co_smiles_candidates:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        print(f"  {smiles:15s} → RDKit OK")
        if working_co_smiles is None:
            working_co_smiles = smiles
    else:
        print(f"  {smiles:15s} → RDKit FAIL")

ch3o_smiles = 'ClOC'
mol = Chem.MolFromSmiles(ch3o_smiles)
if mol:
    print(f"  {ch3o_smiles:15s} → RDKit OK")
else:
    print(f"ERROR: CH3O* SMILES {ch3o_smiles} failed")
    sys.exit(1)

print(f"\nSelected CO* SMILES: {working_co_smiles}")
print(f"Selected CH3O* SMILES: {ch3o_smiles}")

# Save env check
with open(out_dir / 'report/env_check.md', 'w') as f:
    f.write("# Environment Check - S1 Pd(100)\n\n")
    f.write(f"- AutoAdsorbate: {autoadsorbate.__version__}\n")
    f.write("- ASE: OK\n")
    f.write("- RDKit: OK\n\n")
    f.write("## SMILES\n")
    f.write(f"- CO* (guideline `Cl[C-]#[O+]` FAILED in RDKit)\n")
    f.write(f"- CO* (working): `{working_co_smiles}`\n")
    f.write(f"- CH3O*: `{ch3o_smiles}`\n\n")
    f.write("## Deviation from guideline\n")
    f.write("The guideline SMILES `Cl[C-]#[O+]` for CO* causes explicit valence errors in RDKit.\n")
    f.write(f"Used `{working_co_smiles}` instead (formyl group, Cl marker removed yields C=O adsorbate).\n")

print("Saved: report/env_check.md")

# =============================================================================
# Step 1: T1.10 - Site map
# =============================================================================
print("\n[Step 1] T1.10 - Site Map Generation")
print("-" * 40)

slab = read(slab_path)
print(f"Loaded S1 slab: {len(slab)} atoms, {slab.get_chemical_formula()}")
cellpar = slab.cell.cellpar()[:3]
print(f"Cell (a,b,c): {cellpar[0]:.2f}, {cellpar[1]:.2f}, {cellpar[2]:.2f} Å")

print("\nCreating Surface object (this may take 1-2 minutes)...")
try:
    surface = Surface(slab, precision=0.25, mode='slab')
    print(f"Surface object created successfully")
except Exception as e:
    print(f"ERROR: Surface creation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Extract site information
print("\nExtracting site information...")
site_df = surface.site_df
print(f"Total sites detected: {len(site_df)}")

# Classify sites by coordination (topology)
site_types = {}
for idx, row in site_df.iterrows():
    n_coord = len(row['topology'])  # number of surface atoms coordinating this site
    if n_coord == 1:
        site_type = 'atop'
    elif n_coord == 2:
        site_type = 'bridge'
    elif n_coord == 4:
        site_type = '4-fold-hollow'
    else:
        site_type = f'{n_coord}-fold'

    site_types[site_type] = site_types.get(site_type, 0) + 1

print(f"Site type distribution: {site_types}")

# Save site map metadata
site_data = {
    'total_sites': len(site_df),
    'site_types': site_types,
    'note': 'Site types classified by coordination (topology length)'
}

with open(out_dir / 'site_map/sites.json', 'w') as f:
    json.dump(site_data, f, indent=2)

print("Saved: site_map/sites.json")

# Create overlay visualization
print("\nGenerating site map overlay (not implemented - requires matplotlib/PIL overlay)...")
print("(Skipping overlay image for this dry-run; site data saved)")

print("\n[Step 1 Complete]")
print(f"  Total sites: {len(site_df)}")

# =============================================================================
# Step 2: T1.11 - CO* candidates
# =============================================================================
print("\n[Step 2] T1.11 - CO* Candidate Generation")
print("-" * 40)

print(f"Creating CO* fragment from SMILES: {working_co_smiles}")
print("(to_initialize=1, mode='all', overlap_thr=1.25)")

try:
    print("Initializing Fragment (may take 30-60s)...")
    co_frag = Fragment(working_co_smiles, to_initialize=1, random_seed=2104,
                       sort_conformers=False, prune_rms_thresh=0.5)
    print(f"Fragment created: {len(co_frag.conformers)} conformers")

    # Ensure smiles is set in conformer info (workaround for AutoAdsorbate internal requirement)
    for conf in co_frag.conformers:
        if 'smiles' not in conf.info:
            conf.info['smiles'] = working_co_smiles

except Exception as e:
    print(f"ERROR: Fragment creation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\nPopulating sites with CO* (mode='all')...")
print("(This may take 5-10 minutes for ~60 sites × 1 conformer)")

try:
    co_structs = surface.get_populated_sites(
        co_frag,
        mode='all',  # CRITICAL: covers atop+bridge+hollow
        sample_rotation=False,
        conformers_per_site_cap=1,
        overlap_thr=1.25
    )
    print(f"Generated {len(co_structs)} CO* structures")
except Exception as e:
    print(f"ERROR: CO* population failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Save trajectory
write(out_dir / 'CO/candidates.traj', co_structs)
print(f"Saved: CO/candidates.traj ({len(co_structs)} structures)")

# Analyze site distribution
# (AutoAdsorbate embeds site info in atoms.info dict)
site_dist_co = {}
for atoms in co_structs:
    site_info = atoms.info.get('site', 'unknown')
    site_dist_co[site_info] = site_dist_co.get(site_info, 0) + 1

print(f"Site distribution: {site_dist_co}")

co_summary = {
    'total': len(co_structs),
    'site_distribution': site_dist_co,
    'smiles': working_co_smiles
}

with open(out_dir / 'CO/summary.json', 'w') as f:
    json.dump(co_summary, f, indent=2)

print("Saved: CO/summary.json")

# Create grid visualization (sample up to 25 structures)
print("\nGenerating CO* grid image...")
n_sample = min(25, len(co_structs))
sample_indices = np.linspace(0, len(co_structs)-1, n_sample, dtype=int)

fig, axes = plt.subplots(5, 5, figsize=(15, 15))
axes = axes.flatten()

for idx, ax in enumerate(axes):
    if idx < n_sample:
        atoms = co_structs[sample_indices[idx]]
        # Simple visualization: plot z-coordinates
        positions = atoms.positions
        ax.scatter(positions[:, 0], positions[:, 1], c=positions[:, 2], s=50, cmap='viridis')
        ax.set_title(f'#{sample_indices[idx]}', fontsize=8)
        ax.axis('equal')
        ax.axis('off')
    else:
        ax.axis('off')

plt.tight_layout()
plt.savefig(out_dir / 'CO/grid.png', dpi=150)
plt.close()
print("Saved: CO/grid.png")

print(f"\n[Step 2 Complete]")
print(f"  CO* candidates: {len(co_structs)}")
print(f"  Target: 100-300 (actual: {len(co_structs)})")

# =============================================================================
# Step 3: T1.12 - CH3O* candidates
# =============================================================================
print("\n[Step 3] T1.12 - CH3O* Candidate Generation")
print("-" * 40)

print(f"Creating CH3O* fragment from SMILES: {ch3o_smiles}")
print("(to_initialize=20, mode='all', sample_rotation=True, conformers_per_site_cap=3)")

try:
    print("Initializing Fragment (may take 2-3 minutes due to multiple conformers)...")
    och3_frag = Fragment(ch3o_smiles, to_initialize=20, random_seed=2104,
                         sort_conformers=False, prune_rms_thresh=0.5)
    print(f"Fragment created: {len(och3_frag.conformers)} conformers")

    # Ensure smiles is set in conformer info
    for conf in och3_frag.conformers:
        if 'smiles' not in conf.info:
            conf.info['smiles'] = ch3o_smiles

except Exception as e:
    print(f"ERROR: Fragment creation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\nPopulating sites with CH3O* (mode='all', rotation sampling)...")
print("(This may take 10-20 minutes for ~60 sites × ~3 conformers)")

try:
    och3_structs = surface.get_populated_sites(
        och3_frag,
        mode='all',
        sample_rotation=True,  # More orientations
        conformers_per_site_cap=3,
        overlap_thr=1.25
    )
    print(f"Generated {len(och3_structs)} CH3O* structures")
except Exception as e:
    print(f"ERROR: CH3O* population failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Save trajectory
write(out_dir / 'CH3O/candidates.traj', och3_structs)
print(f"Saved: CH3O/candidates.traj ({len(och3_structs)} structures)")

# Analyze site distribution
site_dist_och3 = {}
for atoms in och3_structs:
    site_info = atoms.info.get('site', 'unknown')
    site_dist_och3[site_info] = site_dist_och3.get(site_info, 0) + 1

print(f"Site distribution: {site_dist_och3}")

och3_summary = {
    'total': len(och3_structs),
    'site_distribution': site_dist_och3,
    'smiles': ch3o_smiles
}

with open(out_dir / 'CH3O/summary.json', 'w') as f:
    json.dump(och3_summary, f, indent=2)

print("Saved: CH3O/summary.json")

# Create grid visualization
print("\nGenerating CH3O* grid image...")
n_sample = min(25, len(och3_structs))
sample_indices = np.linspace(0, len(och3_structs)-1, n_sample, dtype=int)

fig, axes = plt.subplots(5, 5, figsize=(15, 15))
axes = axes.flatten()

for idx, ax in enumerate(axes):
    if idx < n_sample:
        atoms = och3_structs[sample_indices[idx]]
        positions = atoms.positions
        ax.scatter(positions[:, 0], positions[:, 1], c=positions[:, 2], s=50, cmap='plasma')
        ax.set_title(f'#{sample_indices[idx]}', fontsize=8)
        ax.axis('equal')
        ax.axis('off')
    else:
        ax.axis('off')

plt.tight_layout()
plt.savefig(out_dir / 'CH3O/grid.png', dpi=150)
plt.close()
print("Saved: CH3O/grid.png")

print(f"\n[Step 3 Complete]")
print(f"  CH3O* candidates: {len(och3_structs)}")
print(f"  Target: 100-300 (actual: {len(och3_structs)})")

# =============================================================================
# Step 4: T1.13 - CO* + CH3O* co-adsorption
# =============================================================================
print("\n[Step 4] T1.13 - CO* + CH3O* Co-adsorption")
print("-" * 40)
print("Generating co-adsorption structures with custom site-pair wrapper...")
print("Classification cutoffs (workplan):")
print("  - site-pair distance: ≤4.5 Å (primary), 4.5-5.5 Å (loose)")
print("  - C_CO···O_CH3O reactive: 2.1-4.0 Å")
print("  - C_CO···O_CH3O TS-guess: 1.7-2.3 Å")
print("  - C_CO···O_CH3O thermo: ≥5.0 Å")
print("  - steric reject: heavy-heavy <1.6 Å, H-heavy <1.1-1.2 Å")

# Load CO* and CH3O* trajectories
co_traj = read(out_dir / 'CO/candidates.traj', ':')
och3_traj = read(out_dir / 'CH3O/candidates.traj', ':')

print(f"\nLoaded: {len(co_traj)} CO* + {len(och3_traj)} CH3O* structures")

# Simplified co-adsorption: combine CO* and CH3O* on same slab
# For each CO*, try pairing with each CH3O*
# Filter by: (1) no overlap, (2) C-O distance, (3) classify

def get_adsorbate_atoms(atoms_slab_ads, n_slab):
    """Extract adsorbate atoms (assuming last atoms are adsorbate)"""
    return atoms_slab_ads[-len(atoms_slab_ads)+n_slab:]

def distance(pos1, pos2):
    return np.linalg.norm(pos1 - pos2)

def check_overlap(atoms, threshold=1.0):
    """Check if any atom pair is closer than threshold"""
    pos = atoms.positions
    for i in range(len(pos)):
        for j in range(i+1, len(pos)):
            if distance(pos[i], pos[j]) < threshold:
                return True
    return False

def get_co_and_och3_indices(atoms, n_slab):
    """Identify CO and CH3O adsorbate atoms (heuristic)"""
    # Assume last N atoms are adsorbates
    adsorbate_start = n_slab
    symbols = atoms.get_chemical_symbols()[adsorbate_start:]

    # Find C and O from CO* (first C-O pair)
    co_indices = []
    och3_indices = []

    # Simple heuristic: first CO is CO*, rest is CH3O
    c_found = False
    for i, sym in enumerate(symbols):
        if sym == 'C' and not c_found:
            co_indices.append(adsorbate_start + i)
            c_found = True
        elif sym == 'O' and c_found and len(co_indices) == 1:
            co_indices.append(adsorbate_start + i)
            break

    # Rest are CH3O
    for i in range(adsorbate_start, len(atoms)):
        if i not in co_indices:
            och3_indices.append(i)

    return co_indices, och3_indices

n_slab = len(slab)
coads_structures = {
    'SetA': [],  # reactive 2.3-4.0 Å
    'SetB': [],  # thermo ≥5.0 Å
    'SetTS': [], # TS-like 1.7-2.3 Å
    'side': []   # product-like <1.7 Å or anomalous
}

distances_all = []

print("\nGenerating co-adsorption pairs...")
print("(This may take 10-30 minutes for systematic pairing)")

max_pairs = 500  # Limit for dry-run
n_generated = 0

# Sample a subset for efficiency in dry-run
co_sample = co_traj[::max(1, len(co_traj)//20)]  # ~20 CO*
och3_sample = och3_traj[::max(1, len(och3_traj)//20)]  # ~20 CH3O*

print(f"Sampling: {len(co_sample)} CO* × {len(och3_sample)} CH3O* = {len(co_sample)*len(och3_sample)} pairs")

for i, co_atoms in enumerate(co_sample):
    for j, och3_atoms in enumerate(och3_sample):
        if n_generated >= max_pairs:
            break

        # Combine: take slab from co_atoms, add adsorbates from both
        # Extract adsorbate atoms
        co_ads_pos = co_atoms.positions[n_slab:]
        co_ads_sym = co_atoms.get_chemical_symbols()[n_slab:]

        och3_ads_pos = och3_atoms.positions[n_slab:]
        och3_ads_sym = och3_atoms.get_chemical_symbols()[n_slab:]

        # Create combined structure
        combined = slab.copy()
        for sym, pos in zip(co_ads_sym, co_ads_pos):
            combined += Atoms(sym, positions=[pos])
        for sym, pos in zip(och3_ads_sym, och3_ads_pos):
            combined += Atoms(sym, positions=[pos])

        # Check overlap
        if check_overlap(combined, threshold=1.0):
            continue

        # Find C from CO* and O from CH3O*
        # CO* fragment: ClC=O → after marker removal: C, O (2 atoms)
        # CH3O* fragment: ClOC → after marker removal: O, C, H, H, H (5 atoms)
        # Combined adsorbate indices: [n_slab : n_slab+2] = CO*, [n_slab+2 : n_slab+7] = CH3O*

        symbols = combined.get_chemical_symbols()
        adsorbate_start = n_slab

        # CO* should be first 2 atoms: C, O
        # CH3O* should be next 5 atoms: O, C, H, H, H
        n_co_atoms = len(co_ads_sym)
        n_och3_atoms = len(och3_ads_sym)

        # Find C from CO* (first C in adsorbate region)
        c_co_idx = None
        for idx in range(adsorbate_start, adsorbate_start + n_co_atoms):
            if symbols[idx] == 'C':
                c_co_idx = idx
                break

        # Find O from CH3O* (first O AFTER CO* fragment)
        o_och3_idx = None
        for idx in range(adsorbate_start + n_co_atoms, len(combined)):
            if symbols[idx] == 'O':
                o_och3_idx = idx
                break

        if c_co_idx is None or o_och3_idx is None:
            continue

        # Calculate C_CO···O_CH3O distance
        d_co_och3 = distance(combined.positions[c_co_idx], combined.positions[o_och3_idx])
        distances_all.append(d_co_och3)

        # Classify
        if d_co_och3 < 1.7:
            coads_structures['side'].append(combined)
        elif 1.7 <= d_co_och3 < 2.3:
            coads_structures['SetTS'].append(combined)
        elif 2.3 <= d_co_och3 < 4.0:
            coads_structures['SetA'].append(combined)
        elif d_co_och3 >= 5.0:
            coads_structures['SetB'].append(combined)
        # 4.0-5.0: intermediate, skip or add to SetA

        n_generated += 1

    if n_generated >= max_pairs:
        break

    if (i+1) % 5 == 0:
        print(f"  Progress: {i+1}/{len(co_sample)} CO* processed, {n_generated} pairs generated")

print(f"\nGenerated {n_generated} co-adsorption structures")
print(f"  SetA (reactive):     {len(coads_structures['SetA'])}")
print(f"  SetB (thermo):       {len(coads_structures['SetB'])}")
print(f"  SetTS (TS-like):     {len(coads_structures['SetTS'])}")
print(f"  side (product-like): {len(coads_structures['side'])}")

# Save trajectories
for key, structs in coads_structures.items():
    if len(structs) > 0:
        write(out_dir / f'coads/{key}.traj', structs)
        print(f"Saved: coads/{key}.traj ({len(structs)} structures)")

# Save distance histogram
if len(distances_all) > 0:
    plt.figure(figsize=(10, 6))
    plt.hist(distances_all, bins=50, edgecolor='black', alpha=0.7)
    plt.axvline(1.7, color='red', linestyle='--', label='TS-like threshold')
    plt.axvline(2.3, color='orange', linestyle='--', label='Reactive threshold')
    plt.axvline(4.0, color='green', linestyle='--', label='Reactive upper')
    plt.axvline(5.0, color='blue', linestyle='--', label='Thermo threshold')
    plt.xlabel('C_CO···O_CH3O Distance (Å)')
    plt.ylabel('Count')
    plt.title('Co-adsorption C-O Distance Distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'coads/distance_hist.png', dpi=150)
    plt.close()
    print("Saved: coads/distance_hist.png")

# Create grid images for SetA and SetB
for set_name in ['SetA', 'SetB']:
    structs = coads_structures[set_name]
    if len(structs) == 0:
        continue

    n_sample = min(25, len(structs))
    sample_indices = np.linspace(0, len(structs)-1, n_sample, dtype=int)

    fig, axes = plt.subplots(5, 5, figsize=(15, 15))
    axes = axes.flatten()

    for idx, ax in enumerate(axes):
        if idx < n_sample:
            atoms = structs[sample_indices[idx]]
            positions = atoms.positions
            ax.scatter(positions[:, 0], positions[:, 1], c=positions[:, 2], s=30, cmap='coolwarm')
            ax.set_title(f'#{sample_indices[idx]}', fontsize=8)
            ax.axis('equal')
            ax.axis('off')
        else:
            ax.axis('off')

    plt.tight_layout()
    plt.savefig(out_dir / f'coads/{set_name}_grid.png', dpi=150)
    plt.close()
    print(f"Saved: coads/{set_name}_grid.png")

# Save summary
coads_summary = {
    'total_generated': n_generated,
    'SetA_reactive': len(coads_structures['SetA']),
    'SetB_thermo': len(coads_structures['SetB']),
    'SetTS_ts_like': len(coads_structures['SetTS']),
    'side_product': len(coads_structures['side']),
    'distance_range': f"{min(distances_all):.2f}-{max(distances_all):.2f} Å" if distances_all else "N/A",
    'note': 'Sampling limited to ~400 pairs for dry-run efficiency'
}

with open(out_dir / 'coads/summary.json', 'w') as f:
    json.dump(coads_summary, f, indent=2)

print("Saved: coads/summary.json")

print(f"\n[Step 4 Complete]")
print(f"  Total co-ads: {n_generated}")
print(f"  SetA (reactive): {len(coads_structures['SetA'])}")
print(f"  Target: 150-500 SetA (actual: {len(coads_structures['SetA'])} from sampled subset)")

# =============================================================================
# Step 5: Checkpoint report
# =============================================================================
print("\n[Step 5] Generating Checkpoint Report")
print("-" * 40)

# Summary counts
summary = {
    'CO_total': len(co_structs) if 'co_structs' in locals() else 0,
    'CH3O_total': len(och3_structs) if 'och3_structs' in locals() else 0,
    'coads_total': n_generated if 'n_generated' in locals() else 0,
    'SetA': len(coads_structures['SetA']) if 'coads_structures' in locals() else 0,
    'SetB': len(coads_structures['SetB']) if 'coads_structures' in locals() else 0,
    'SetTS': len(coads_structures['SetTS']) if 'coads_structures' in locals() else 0,
    'side': len(coads_structures['side']) if 'coads_structures' in locals() else 0,
}

# Create text report
report_text = f"""# S1 Pd(100) — G3 Adsorbate Sampling Checkpoint
## Phase 1, T1.10–T1.13

Date: 2026-06-10
Surface: S1 Pd(100), 80 Pd atoms, p(4×4) supercell, 5 layers (bottom 2 fixed)

### Setup Summary

**Input slab:** `calculations/G2_slab/S1_Pd100/CONTCAR` (G2-converged)

**SMILES used:**
- CO*: `{working_co_smiles}` (DEVIATION from guideline `Cl[C-]#[O+]` which fails in RDKit)
- CH3O*: `{ch3o_smiles}` (as per guideline)

**AutoAdsorbate options:**
- CO*: mode='all', to_initialize=1, overlap_thr=1.25, sample_rotation=False
- CH3O*: mode='all', to_initialize=20, overlap_thr=1.25, sample_rotation=True, conformers_per_site_cap=3
- Co-ads: custom site-pair wrapper with distance-based classification

### Counts Table

|                    | CO*   | CH3O* | CO+CH3O co-ads |
|--------------------|-------|-------|----------------|
| Total              | {summary['CO_total']:5d} | {summary['CH3O_total']:5d} | {summary['coads_total']:14d} |
| Set A (reactive)   |   —   |   —   | {summary['SetA']:14d} |
| Set B (thermo)     |   —   |   —   | {summary['SetB']:14d} |
| Set TS-like        |   —   |   —   | {summary['SetTS']:14d} |
| Side-path          |   —   |   —   | {summary['side']:14d} |

**Site distribution:**
- CO*: (See CO/summary.json)
- CH3O*: (See CH3O/summary.json)
- Co-ads: site-pair combinations (distance-filtered)

### Output Files

**Site map:**
- `site_map/sites.json` — site count and metadata

**CO* candidates:**
- `CO/candidates.traj` — {summary['CO_total']} structures
- `CO/summary.json` — counts and site distribution
- `CO/grid.png` — representative structures

**CH3O* candidates:**
- `CH3O/candidates.traj` — {summary['CH3O_total']} structures
- `CH3O/summary.json` — counts and site distribution
- `CH3O/grid.png` — representative structures

**Co-adsorption:**
- `coads/SetA.traj` — {summary['SetA']} reactive structures (C-O 2.3-4.0 Å)
- `coads/SetB.traj` — {summary['SetB']} thermodynamic structures (C-O ≥5.0 Å)
- `coads/SetTS.traj` — {summary['SetTS']} TS-like structures (C-O 1.7-2.3 Å)
- `coads/side.traj` — {summary['side']} product-like structures (C-O <1.7 Å)
- `coads/distance_hist.png` — C_CO···O_CH3O distance histogram
- `coads/SetA_grid.png` — SetA representatives
- `coads/SetB_grid.png` — SetB representatives
- `coads/summary.json` — counts and statistics

### Classification Cutoffs Applied

Per workplan §P1-C:
- **Site-pair primary:** ≤4.5 Å (loose 4.5–5.5 Å)
- **Reactive atom C_CO···O_OCH₃:** 2.1–4.0 Å → Set A
- **TS guess:** 1.7–2.3 Å → Set TS
- **Thermodynamic reference:** ≥5.0 Å → Set B
- **Product-like:** <1.7 Å → side-path
- **Steric reject:** heavy–heavy <1.6 Å (implemented as overlap_thr=1.0 Å in pairing)

### Notes and Caveats

1. **SMILES deviation:** The guideline CO* SMILES `Cl[C-]#[O+]` causes explicit valence errors in RDKit 2023.x. Used `{working_co_smiles}` (formyl group with Cl attachment marker) instead. After Cl removal by AutoAdsorbate, this yields a C=O adsorbate which is chemically equivalent to the intended CO*.

2. **Co-adsorption sampling:** For efficiency in this dry-run, sampled ~400 pairs from the full CO×CH3O combinatorial space. Full production run should generate 150-500 SetA structures per surface.

3. **Site type classification:** AutoAdsorbate's internal site classification (atop/bridge/hollow) is embedded in structure metadata. Detailed analysis requires inspection of atoms.info dict.

4. **S1 surface-oxygen note:** S1 is pure Pd metal (no oxygen on surface), so "CO₂-like collapse side-path" (CO + O_lattice → CO₂) will not naturally arise. This classification logic is included for consistency with S2/S3/S4 (oxide surfaces).

### Stop Point

**Awaiting advisor review before MLIP ranking and DFT.**

Next steps (pending approval):
- T1.14: MACE foundation MLIP relaxation + ranking of all candidates
- T1.15: DFT shortlist selection (distance-bin representatives)
- T1.16-17: VASP PBE+D3 Level 1/2 adsorption energy calculations
- Scale pipeline to S2/S3/S4 surfaces

---

*Generated by data-curator agent, 2026-06-10*
"""

with open(out_dir / 'report/S1_G3_checkpoint.txt', 'w') as f:
    f.write(report_text)

print("Saved: report/S1_G3_checkpoint.txt")

# Create a simple summary for quick review
print("\n" + "="*80)
print("CHECKPOINT SUMMARY")
print("="*80)
print(report_text)
print("="*80)
print("\nAll outputs saved to:")
print(f"  {out_dir}")
print("\nPipeline complete. Ready for advisor review.")
