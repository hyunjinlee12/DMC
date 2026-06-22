#!/usr/bin/env python3
"""
Co-adsorption candidate generation following guideline specification.
Replaces stride-based sampling with exact guideline cutoffs.

Input: existing CO/candidates.traj and CH3O/candidates.traj per surface
Output: coads_guide/{SetA,SetB,SetTS,side,rejected}.traj + summary + viz
"""

import numpy as np
from pathlib import Path
from ase.io import read, write
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

def get_adsorbate_indices(atoms):
    """
    Extract adsorbate atom indices from an adsorbed structure.
    Adsorbate = non-metal atoms (C, O, H) that are not part of slab.
    For Pd surfaces, adsorbate = C, O, H atoms.
    For PdO surfaces, adsorbate = C, H atoms (O could be slab).

    Simple heuristic: take the last N atoms where N is number of consecutive
    non-Pd atoms from the end (works because AutoAdsorbate adds adsorbate last).
    """
    syms = atoms.get_chemical_symbols()

    # Find adsorbate as trailing non-Pd atoms
    ads_indices = []
    for i in range(len(atoms) - 1, -1, -1):
        if syms[i] in ['C', 'H']:  # C and H are definitely adsorbate
            ads_indices.insert(0, i)
        elif syms[i] == 'O':
            # O could be adsorbate or slab - include if it's in the trailing group
            # Check if we already have C/H (adsorbate started)
            if ads_indices:
                ads_indices.insert(0, i)
            else:
                # First O from end - assume it's adsorbate if we're looking for CH3O or CO
                # Check the next atom
                if i < len(atoms) - 1 and syms[i+1] in ['C', 'H']:
                    ads_indices.insert(0, i)
                else:
                    # Isolated O at end might be adsorbate
                    ads_indices.insert(0, i)
        else:
            # Hit a Pd - stop
            break

    return ads_indices

def compute_reactive_distance(co_struct, ch3o_struct):
    """
    Compute C_CO - O_CH3O distance.
    CO: C is bonded to surface, O is terminal → take C (bottom of CO)
    CH3O: O is bonded to surface → take O (bottom of CH3O)
    """
    co_ads = get_adsorbate_indices(co_struct)
    ch3o_ads = get_adsorbate_indices(ch3o_struct)

    # CO: find C (should be closer to surface or heavier)
    co_atoms_only = co_struct[co_ads]
    c_idx_local = None
    for i, sym in enumerate(co_atoms_only.get_chemical_symbols()):
        if sym == 'C':
            c_idx_local = i
            break
    if c_idx_local is None:
        # Fallback: take lower Z position
        c_idx_local = np.argmin(co_atoms_only.positions[:, 2])
    c_pos = co_atoms_only[c_idx_local].position

    # CH3O: find O (should be lowest Z or first O)
    ch3o_atoms_only = ch3o_struct[ch3o_ads]
    o_idx_local = None
    for i, sym in enumerate(ch3o_atoms_only.get_chemical_symbols()):
        if sym == 'O':
            o_idx_local = i
            break
    if o_idx_local is None:
        # Fallback: lowest Z
        o_idx_local = np.argmin(ch3o_atoms_only.positions[:, 2])
    o_pos = ch3o_atoms_only[o_idx_local].position

    # Distance with PBC
    delta = c_pos - o_pos
    cell = co_struct.get_cell()
    # Minimum image convention
    for i in range(3):
        if cell[i, i] > 0:
            delta[i] -= cell[i, i] * np.round(delta[i] / cell[i, i])

    return np.linalg.norm(delta)

def check_steric_clash(co_struct, ch3o_struct, heavy_thr=1.6, h_heavy_thr=1.1):
    """
    Check steric overlap between CO and CH3O adsorbates.
    Returns True if clash detected (→ reject), False if OK.
    """
    co_ads = get_adsorbate_indices(co_struct)
    ch3o_ads = get_adsorbate_indices(ch3o_struct)

    co_pos = co_struct[co_ads].positions
    co_sym = co_struct[co_ads].get_chemical_symbols()
    ch3o_pos = ch3o_struct[ch3o_ads].positions
    ch3o_sym = ch3o_struct[ch3o_ads].get_chemical_symbols()

    cell = co_struct.get_cell()

    for i, (p1, s1) in enumerate(zip(co_pos, co_sym)):
        for j, (p2, s2) in enumerate(zip(ch3o_pos, ch3o_sym)):
            delta = p1 - p2
            # PBC
            for k in range(3):
                if cell[k, k] > 0:
                    delta[k] -= cell[k, k] * np.round(delta[k] / cell[k, k])
            dist = np.linalg.norm(delta)

            # Check thresholds
            if s1 == 'H' or s2 == 'H':
                if dist < h_heavy_thr:
                    return True
            else:
                if dist < heavy_thr:
                    return True
    return False

def merge_structures(co_struct, ch3o_struct):
    """
    Create combined structure with both CO and CH3O adsorbates.
    """
    from ase import Atoms

    # Start with clean slab (fixed atoms from co_struct)
    if hasattr(co_struct, 'constraints') and len(co_struct.constraints) > 0:
        from ase.constraints import FixAtoms
        fixed = set()
        for c in co_struct.constraints:
            if isinstance(c, FixAtoms):
                fixed.update(c.get_indices())
        slab_indices = sorted(fixed)
    else:
        # Assume slab is all atoms except adsorbate
        co_ads = get_adsorbate_indices(co_struct)
        slab_indices = [i for i in range(len(co_struct)) if i not in co_ads]

    slab_atoms = co_struct[slab_indices]
    co_ads_atoms = co_struct[get_adsorbate_indices(co_struct)]
    ch3o_ads_atoms = ch3o_struct[get_adsorbate_indices(ch3o_struct)]

    # Combine
    combined = slab_atoms + co_ads_atoms + ch3o_ads_atoms
    combined.set_cell(co_struct.get_cell())
    combined.set_pbc(co_struct.get_pbc())

    # Re-apply constraints
    from ase.constraints import FixAtoms
    combined.set_constraint(FixAtoms(indices=list(range(len(slab_atoms)))))

    return combined

def generate_coadsorption(surface_name, co_traj, ch3o_traj, out_dir):
    """
    Generate co-adsorption candidates following guideline cutoffs.

    Cutoffs (guideline §9):
    - site-pair primary: ≤ 4.5 Å
    - reactive C_CO···O_CH3O: 2.1 – 4.0 Å (Set A)
    - TS guess C···O: 1.7 – 2.3 Å (Set TS)
    - thermodynamic reference: ≥ 5.0 Å (Set B)
    - steric reject: heavy-heavy < 1.6 Å, H-heavy < 1.1 Å
    - side-path: < 1.7 Å

    Note: Set A (2.1-4.0) overlaps with Set TS (1.7-2.3) in range 2.1-2.3.
          Guideline allows this overlap.
    """
    print(f"\n{'='*60}")
    print(f"Processing {surface_name}")
    print(f"{'='*60}")

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Load candidates
    co_structs = read(co_traj, ':')
    ch3o_structs = read(ch3o_traj, ':')

    n_co = len(co_structs)
    n_ch3o = len(ch3o_structs)
    total_pairs = n_co * n_ch3o

    print(f"CO candidates: {n_co}")
    print(f"CH3O candidates: {n_ch3o}")
    print(f"Total pair combinations: {total_pairs}")

    # Step 1: Enumerate all pairs + compute reactive distance
    print("\nStep 1: Enumerating all pairs and computing C_CO-O_CH3O distances...")

    pairs = []
    for i, co_s in enumerate(co_structs):
        for j, ch3o_s in enumerate(ch3o_structs):
            dist = compute_reactive_distance(co_s, ch3o_s)
            pairs.append({
                'co_idx': i,
                'ch3o_idx': j,
                'distance': dist,
                'co_struct': co_s,
                'ch3o_struct': ch3o_s
            })

    print(f"Total pairs enumerated: {len(pairs)}")

    # Step 2: Primary filter (≤ 4.5 Å for reactive pool, ≥ 5.0 Å for thermo pool)
    print("\nStep 2: Applying site-pair primary filter...")

    primary_pool = [p for p in pairs if p['distance'] <= 4.5]
    thermo_pool = [p for p in pairs if p['distance'] >= 5.0]
    ambiguous = [p for p in pairs if 4.5 < p['distance'] < 5.0]

    print(f"Primary pool (≤ 4.5 Å): {len(primary_pool)}")
    print(f"Thermo pool (≥ 5.0 Å): {len(thermo_pool)}")
    print(f"Ambiguous (4.5-5.0 Å, excluded): {len(ambiguous)}")

    # Step 3: max_per_pair=2 (limit same site-pair combinations)
    # For simplicity, we group by (co_idx, ch3o_idx) which already represents unique pairs
    # Since we're doing all×all, each (i,j) is unique → this filter doesn't reduce count
    # But we keep it for consistency with guideline

    # In practice, max_per_pair would apply if we had multiple conformers per site
    # AutoAdsorbate's candidates.traj may have multiple structures per site
    # We'll assume co_idx and ch3o_idx already represent site indices
    # (CO likely 1 per site, CH3O has 3 conformers per site)

    # Extract site indices from candidate index
    # Heuristic: CO index = site index (1:1)
    #            CH3O index // 3 = site index (3 conformers per site)

    def get_site_pair(p):
        co_site = p['co_idx']
        ch3o_site = p['ch3o_idx'] // 3  # Assume 3 conformers per site
        return (co_site, ch3o_site)

    # Group primary pool by site-pair
    print("\nStep 3: Applying max_per_pair=2 filter...")

    primary_by_site = defaultdict(list)
    for p in primary_pool:
        site_pair = get_site_pair(p)
        primary_by_site[site_pair].append(p)

    # Keep max 2 per site-pair (closest + farthest for diversity)
    primary_filtered = []
    for site_pair, group in primary_by_site.items():
        if len(group) <= 2:
            primary_filtered.extend(group)
        else:
            # Keep closest and farthest
            sorted_group = sorted(group, key=lambda x: x['distance'])
            primary_filtered.extend([sorted_group[0], sorted_group[-1]])

    print(f"After max_per_pair=2: {len(primary_filtered)} (from {len(primary_pool)})")

    # Apply same to thermo pool
    thermo_by_site = defaultdict(list)
    for p in thermo_pool:
        site_pair = get_site_pair(p)
        thermo_by_site[site_pair].append(p)

    thermo_filtered = []
    for site_pair, group in thermo_by_site.items():
        if len(group) <= 2:
            thermo_filtered.extend(group)
        else:
            sorted_group = sorted(group, key=lambda x: x['distance'])
            thermo_filtered.extend([sorted_group[0], sorted_group[-1]])

    print(f"Thermo pool after max_per_pair=2: {len(thermo_filtered)} (from {len(thermo_pool)})")

    # Step 4: Steric clash filter
    print("\nStep 4: Checking steric clashes...")

    primary_ok = []
    primary_rejected = []
    for p in primary_filtered:
        if check_steric_clash(p['co_struct'], p['ch3o_struct']):
            primary_rejected.append(p)
        else:
            primary_ok.append(p)

    thermo_ok = []
    thermo_rejected = []
    for p in thermo_filtered:
        if check_steric_clash(p['co_struct'], p['ch3o_struct']):
            thermo_rejected.append(p)
        else:
            thermo_ok.append(p)

    print(f"Primary pool after steric filter: {len(primary_ok)} (rejected: {len(primary_rejected)})")
    print(f"Thermo pool after steric filter: {len(thermo_ok)} (rejected: {len(thermo_rejected)})")

    # Step 5: Classification (guideline cutoffs with overlap)
    print("\nStep 5: Classifying into Sets...")

    set_side = []      # < 1.7 Å
    set_ts = []        # 1.7 – 2.3 Å (TS guess)
    set_a = []         # 2.1 – 4.0 Å (reactive)
    set_b = []         # ≥ 5.0 Å (thermo)

    # Primary pool classification
    for p in primary_ok:
        d = p['distance']

        # Side-path (< 1.7 Å)
        if d < 1.7:
            set_side.append(p)

        # Set TS (1.7 – 2.3 Å)
        if 1.7 <= d <= 2.3:
            set_ts.append(p)

        # Set A (2.1 – 4.0 Å)
        if 2.1 <= d <= 4.0:
            set_a.append(p)

    # Thermo pool → Set B
    set_b = thermo_ok

    # All rejected
    rejected_all = primary_rejected + thermo_rejected

    print(f"\nFinal classification:")
    print(f"  Set A (reactive, 2.1-4.0 Å): {len(set_a)}")
    print(f"  Set TS (TS guess, 1.7-2.3 Å): {len(set_ts)}")
    print(f"  Set B (thermo, ≥5.0 Å): {len(set_b)}")
    print(f"  Side-path (<1.7 Å): {len(set_side)}")
    print(f"  Rejected (steric): {len(rejected_all)}")

    # Note overlap
    overlap_count = len([p for p in primary_ok if 2.1 <= p['distance'] <= 2.3])
    print(f"  (Note: {overlap_count} structures in both Set A and Set TS, 2.1-2.3 Å overlap)")

    # Build merged structures and save
    print("\nBuilding merged structures...")

    def build_set(pair_list, name):
        structs = []
        for p in pair_list:
            merged = merge_structures(p['co_struct'], p['ch3o_struct'])
            merged.info['co_idx'] = p['co_idx']
            merged.info['ch3o_idx'] = p['ch3o_idx']
            merged.info['reactive_distance'] = p['distance']
            structs.append(merged)
        if structs:
            write(out_path / f"{name}.traj", structs)
            print(f"  Saved {len(structs)} to {name}.traj")
        return structs

    set_a_structs = build_set(set_a, 'SetA')
    set_ts_structs = build_set(set_ts, 'SetTS')
    set_b_structs = build_set(set_b, 'SetB')
    set_side_structs = build_set(set_side, 'side')
    rejected_structs = build_set(rejected_all, 'rejected')

    # Summary statistics
    summary = {
        'surface': surface_name,
        'n_co_candidates': n_co,
        'n_ch3o_candidates': n_ch3o,
        'total_pairs_enumerated': total_pairs,
        'after_primary_filter': len(primary_pool),
        'after_max_per_pair': len(primary_filtered),
        'after_steric': len(primary_ok),
        'thermo_pool_initial': len(thermo_pool),
        'thermo_after_max_per_pair': len(thermo_filtered),
        'thermo_after_steric': len(thermo_ok),
        'final_counts': {
            'SetA_reactive': len(set_a),
            'SetTS_guess': len(set_ts),
            'SetB_thermo': len(set_b),
            'side_path': len(set_side),
            'rejected_steric': len(rejected_all)
        },
        'overlap_SetA_SetTS': overlap_count,
        'distances': {
            'SetA': {
                'min': float(min([p['distance'] for p in set_a])) if set_a else None,
                'max': float(max([p['distance'] for p in set_a])) if set_a else None,
                'mean': float(np.mean([p['distance'] for p in set_a])) if set_a else None
            },
            'SetTS': {
                'min': float(min([p['distance'] for p in set_ts])) if set_ts else None,
                'max': float(max([p['distance'] for p in set_ts])) if set_ts else None,
                'mean': float(np.mean([p['distance'] for p in set_ts])) if set_ts else None
            },
            'SetB': {
                'min': float(min([p['distance'] for p in set_b])) if set_b else None,
                'max': float(max([p['distance'] for p in set_b])) if set_b else None,
                'mean': float(np.mean([p['distance'] for p in set_b])) if set_b else None
            },
            'side': {
                'min': float(min([p['distance'] for p in set_side])) if set_side else None,
                'max': float(max([p['distance'] for p in set_side])) if set_side else None,
                'mean': float(np.mean([p['distance'] for p in set_side])) if set_side else None
            }
        }
    }

    with open(out_path / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nSaved summary to summary.json")

    # Generate visualizations
    print("\nGenerating visualizations...")

    # Distance histogram
    fig, ax = plt.subplots(figsize=(10, 6))

    all_dists = [p['distance'] for p in pairs]
    primary_dists = [p['distance'] for p in primary_ok]
    thermo_dists = [p['distance'] for p in thermo_ok]

    ax.hist(all_dists, bins=50, alpha=0.3, label='All pairs', color='gray')
    ax.hist(primary_dists, bins=50, alpha=0.6, label=f'Primary (≤4.5 Å, n={len(primary_ok)})', color='blue')
    ax.hist(thermo_dists, bins=50, alpha=0.6, label=f'Thermo (≥5.0 Å, n={len(thermo_ok)})', color='green')

    # Mark cutoff regions
    ax.axvline(1.7, color='red', linestyle='--', linewidth=1, label='side | TS')
    ax.axvline(2.1, color='orange', linestyle='--', linewidth=1, label='TS | Set A')
    ax.axvline(2.3, color='purple', linestyle='--', linewidth=1, label='Set TS | Set A')
    ax.axvline(4.0, color='brown', linestyle='--', linewidth=1, label='Set A | ambiguous')
    ax.axvline(4.5, color='black', linestyle='--', linewidth=2, label='primary | ambiguous')
    ax.axvline(5.0, color='darkgreen', linestyle='--', linewidth=2, label='ambiguous | thermo')

    ax.set_xlabel('C_CO - O_CH3O distance (Å)')
    ax.set_ylabel('Count')
    ax.set_title(f'{surface_name} Co-adsorption Distance Distribution')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path / 'distance_hist.png', dpi=150)
    plt.close()

    print(f"  Saved distance_hist.png")

    # Representative structure grids (25 each for Set A and Set B)
    def plot_grid(structs, title, filename):
        if not structs:
            return

        n_plot = min(25, len(structs))
        indices = np.linspace(0, len(structs)-1, n_plot, dtype=int)

        fig, axes = plt.subplots(5, 5, figsize=(15, 15))
        axes = axes.flatten()

        for idx, i in enumerate(indices):
            ax = axes[idx]
            atoms = structs[i]

            # Simple top-down view visualization
            pos = atoms.get_positions()
            syms = atoms.get_chemical_symbols()

            # Color by element
            colors = []
            for s in syms:
                if s == 'Pd':
                    colors.append('gray')
                elif s == 'O':
                    colors.append('red')
                elif s == 'C':
                    colors.append('black')
                elif s == 'H':
                    colors.append('white')
                else:
                    colors.append('blue')

            ax.scatter(pos[:, 0], pos[:, 1], c=colors, s=10, alpha=0.6)
            ax.set_aspect('equal')
            ax.set_title(f"#{i} d={atoms.info['reactive_distance']:.2f}Å", fontsize=8)
            ax.axis('off')

        # Hide unused subplots
        for idx in range(n_plot, 25):
            axes[idx].axis('off')

        plt.suptitle(title, fontsize=14)
        plt.tight_layout()
        plt.savefig(out_path / filename, dpi=100)
        plt.close()

        print(f"  Saved {filename}")

    plot_grid(set_a_structs, f'{surface_name} Set A Representatives', 'SetA_grid.png')
    plot_grid(set_b_structs, f'{surface_name} Set B Representatives', 'SetB_grid.png')

    # Generate report
    report = f"""# Co-adsorption Generation Report: {surface_name}

## Summary

Generated co-adsorption candidates following guideline specification (§P1-C, §9).

### Input
- CO candidates: {n_co}
- CH3O candidates: {n_ch3o}

### Processing Pipeline

**Step 1: Enumerate all pairs**
- Total pair combinations: {total_pairs}

**Step 2: Primary filter (site-pair distance)**
- Primary pool (≤ 4.5 Å): {len(primary_pool)}
- Thermo pool (≥ 5.0 Å): {len(thermo_pool)}
- Ambiguous (4.5-5.0 Å, excluded): {len(ambiguous)}

**Step 3: max_per_pair=2**
- Primary pool after cap: {len(primary_filtered)} (from {len(primary_pool)})
- Thermo pool after cap: {len(thermo_filtered)} (from {len(thermo_pool)})

**Step 4: Steric clash filter**
- Primary pool passed: {len(primary_ok)}
- Primary pool rejected: {len(primary_rejected)}
- Thermo pool passed: {len(thermo_ok)}
- Thermo pool rejected: {len(thermo_rejected)}

**Step 5: Classification (guideline cutoffs)**

| Set | Cutoff | Count | Distance (Å) |
|-----|--------|-------|--------------|
| Side-path | < 1.7 Å | {len(set_side)} | {'%.2f - %.2f (mean %.2f)' % (summary['distances']['side']['min'], summary['distances']['side']['max'], summary['distances']['side']['mean']) if set_side else 'N/A'} |
| Set TS | 1.7 - 2.3 Å | {len(set_ts)} | {'%.2f - %.2f (mean %.2f)' % (summary['distances']['SetTS']['min'], summary['distances']['SetTS']['max'], summary['distances']['SetTS']['mean']) if set_ts else 'N/A'} |
| Set A | 2.1 - 4.0 Å | {len(set_a)} | {'%.2f - %.2f (mean %.2f)' % (summary['distances']['SetA']['min'], summary['distances']['SetA']['max'], summary['distances']['SetA']['mean']) if set_a else 'N/A'} |
| Set B | ≥ 5.0 Å | {len(set_b)} | {'%.2f - %.2f (mean %.2f)' % (summary['distances']['SetB']['min'], summary['distances']['SetB']['max'], summary['distances']['SetB']['mean']) if set_b else 'N/A'} |
| Rejected | steric clash | {len(rejected_all)} | N/A |

**Note:** Set A and Set TS overlap in range 2.1-2.3 Å per guideline specification.
Number of structures in overlap region: {overlap_count}

### Target Compliance

Guideline target: 150-500 candidates per surface (Phase 1 co-adsorption)

- **Set A (reactive)**: {len(set_a)} ✓ {'within target' if 150 <= len(set_a) <= 500 else 'outside target'}
- **Set B (thermo)**: {len(set_b)} ✓ {'within target' if len(set_b) >= 10 else 'low count'}
- **Combined relevant (A+TS+B-overlap)**: {len(set_a) + len(set_ts) + len(set_b) - overlap_count}

### Outputs

- `SetA.traj`: {len(set_a)} reactive pair candidates (2.1-4.0 Å)
- `SetTS.traj`: {len(set_ts)} TS guess candidates (1.7-2.3 Å)
- `SetB.traj`: {len(set_b)} thermodynamic reference candidates (≥5.0 Å)
- `side.traj`: {len(set_side)} side-path candidates (<1.7 Å)
- `rejected.traj`: {len(rejected_all)} steric clash rejects
- `summary.json`: full statistics
- `distance_hist.png`: C_CO-O_CH3O distance distribution
- `SetA_grid.png`: 25 Set A representatives
- `SetB_grid.png`: 25 Set B representatives

### Cutoffs Applied (Guideline §9)

| Parameter | Cutoff |
|-----------|--------|
| site-pair primary | ≤ 4.5 Å |
| reactive C_CO···O_OCH3 | 2.1 - 4.0 Å |
| TS guess C···O | 1.7 - 2.3 Å |
| thermodynamic reference | ≥ 5.0 Å |
| steric heavy-heavy reject | < 1.6 Å |
| steric H-heavy reject | < 1.1 Å |

---
Generated: {surface_name} co-adsorption candidates (guideline specification)
"""

    with open(out_path / 'report.md', 'w') as f:
        f.write(report)

    print(f"  Saved report.md")

    print(f"\n{'='*60}")
    print(f"Completed {surface_name}")
    print(f"{'='*60}\n")

    return summary

if __name__ == '__main__':
    base_dir = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G3_adsorption')

    surfaces = [
        'S1_Pd100',
        'S2_PdO101_Pd100',
        'S3_PdO100',
        'S3b_PdO100_PdOterm',
        'S4_PdO2_110'
    ]

    results = {}

    for surf in surfaces:
        surf_dir = base_dir / surf
        co_traj = surf_dir / 'CO' / 'candidates.traj'
        ch3o_traj = surf_dir / 'CH3O' / 'candidates.traj'
        out_dir = surf_dir / 'coads_guide'

        if not co_traj.exists() or not ch3o_traj.exists():
            print(f"Skipping {surf}: missing candidates")
            continue

        try:
            summary = generate_coadsorption(surf, co_traj, ch3o_traj, out_dir)
            results[surf] = summary
        except Exception as e:
            print(f"ERROR processing {surf}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Overall summary
    print("\n" + "="*60)
    print("OVERALL SUMMARY")
    print("="*60)

    for surf, res in results.items():
        print(f"\n{surf}:")
        print(f"  Total pairs enumerated: {res['total_pairs_enumerated']}")
        print(f"  After primary filter: {res['after_primary_filter']}")
        print(f"  After max_per_pair: {res['after_max_per_pair']}")
        print(f"  After steric: {res['after_steric']}")
        print(f"  Final counts:")
        for k, v in res['final_counts'].items():
            print(f"    {k}: {v}")
        print(f"  Target (150-500): {'✓ PASS' if 150 <= res['final_counts']['SetA_reactive'] <= 500 else '✗ outside'}")

    # Save combined summary
    with open(base_dir / 'coads_guide_summary.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved combined summary to coads_guide_summary.json")
    print("\nAll surfaces processed!")
