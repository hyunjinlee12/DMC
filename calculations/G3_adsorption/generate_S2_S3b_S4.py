"""G3 adsorption sampling pipeline for S2, S3b, S4 surfaces.

Extends the S1/S3 workflow with updated guideline-compliant cutoffs:
- side-path: < 1.7 Å
- Set TS: 1.7 – 2.1 Å
- Set A reactive: 2.1 – 4.0 Å (guideline P1-C)
- ambiguous (reject): 4.0 – 5.0 Å
- Set B thermo: >= 5.0 Å

Steric reject: heavy-heavy < 1.6 Å OR H-heavy < 1.1 Å

NO MLIP, NO DFT — checkpoint stop after structure generation.
"""
import sys, os, json, numpy as np
from pathlib import Path
from ase import Atoms
from ase.io import read, write
from autoadsorbate import Surface, Fragment
import matplotlib.pyplot as plt
from collections import Counter
import random

random.seed(2104)
np.random.seed(2104)

ROOT = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc")
G3_BASE = ROOT / "calculations/G3_adsorption"

# Surface configurations: (name, slab_path, n_sub)
SURFACES = [
    ("S2_PdO101_Pd100", ROOT / "calculations/G2_slab/S2_PdO101_Pd100/CONTCAR", 112),
    ("S3b_PdO100_PdOterm", ROOT / "calculations/G2_slab/S3b_PdO100_PdOterm/CONTCAR", 104),
    ("S4_PdO2_110", ROOT / "calculations/G2_slab/S4_PdO2_110/CONTCAR", 144),
]

# Updated guideline-compliant cutoffs
CUTOFF_SIDE = 1.7
CUTOFF_TS_LOW = 1.7
CUTOFF_TS_HIGH = 2.1
CUTOFF_A_LOW = 2.1
CUTOFF_A_HIGH = 4.0
CUTOFF_AMB_HIGH = 5.0  # above this = Set B
STERIC_HEAVY = 1.6
STERIC_H = 1.1

def find_first_atom(atoms, element, n_sub):
    """Find first occurrence of element in adsorbate region (index >= n_sub)."""
    syms = atoms.get_chemical_symbols()
    for i in range(n_sub, len(atoms)):
        if syms[i] == element:
            return i
    return None

def check_steric_clash(combined, n_sub, n_co_total):
    """Check for steric overlap between CO* and CH3O* adsorbates.

    Args:
        combined: Atoms object with substrate + CO* + CH3O*
        n_sub: number of substrate atoms
        n_co_total: total atoms including substrate + CO* (= len(a_co) in generate_coads)

    Structure:
        atoms[0:n_sub] = substrate
        atoms[n_sub:n_co_total] = CO* fragment
        atoms[n_co_total:] = CH3O* fragment

    Returns True if reject (clash detected between CO* and CH3O*), False if OK.
    Heavy-heavy < 1.6 Å OR H-heavy < 1.1 Å → reject.
    """
    syms = combined.get_chemical_symbols()

    co_indices = list(range(n_sub, n_co_total))
    ch3o_indices = list(range(n_co_total, len(combined)))

    # Check distances between CO* atoms and CH3O* atoms
    for i in co_indices:
        for j in ch3o_indices:
            d = combined.get_distance(i, j, mic=True)

            # Heavy-heavy check
            if syms[i] != 'H' and syms[j] != 'H':
                if d < STERIC_HEAVY:
                    return True
            # H-heavy check
            elif (syms[i] == 'H' and syms[j] != 'H') or (syms[i] != 'H' and syms[j] == 'H'):
                if d < STERIC_H:
                    return True
    return False

def generate_site_map(surface_name, slab_path, n_sub):
    """Step 1: Generate site map."""
    print(f"\n{'='*60}")
    print(f"[{surface_name}] Step 1: Site map generation")
    print(f"{'='*60}")

    work_dir = G3_BASE / surface_name
    site_dir = work_dir / "site_map"
    site_dir.mkdir(parents=True, exist_ok=True)

    slab = read(slab_path)
    print(f"  Slab loaded: {len(slab)} atoms")
    print(f"  Formula: {slab.get_chemical_formula()}")

    surface = Surface(slab, mode='slab', precision=0.25)
    site_df = surface.site_df

    print(f"  Total sites detected: {len(site_df)}")

    # Classify sites by coordination (topology length)
    site_types = {}
    for idx, row in site_df.iterrows():
        n_coord = len(row['topology'])
        if n_coord == 1:
            site_type = 'atop'
        elif n_coord == 2:
            site_type = 'bridge'
        elif n_coord == 3:
            site_type = '3-fold-hollow'
        elif n_coord == 4:
            site_type = '4-fold-hollow'
        else:
            site_type = f'{n_coord}-fold'
        site_types[site_type] = site_types.get(site_type, 0) + 1

    print(f"  Site type distribution: {site_types}")

    # Save summary
    summary = {
        "total_sites": len(site_df),
        "site_types": site_types,
        "slab_atoms": len(slab),
        "slab_formula": slab.get_chemical_formula(),
    }

    with open(site_dir / "sites.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Saved: {site_dir / 'sites.json'}")

    # Simple visualization: top view (slab atoms only, site overlay requires AutoAdsorbate API check)
    try:
        fig, ax = plt.subplots(figsize=(10, 10))

        # Plot slab atoms (top view, xy projection)
        pos = slab.get_positions()
        syms = slab.get_chemical_symbols()

        colors = {'Pd': 'steelblue', 'O': 'red', 'C': 'gray', 'H': 'lightgray'}
        for i, (p, s) in enumerate(zip(pos, syms)):
            c = colors.get(s, 'black')
            ax.scatter(p[0], p[1], c=c, s=50, alpha=0.5, edgecolor='k', linewidth=0.5)

        ax.set_xlabel('x (Å)')
        ax.set_ylabel('y (Å)')
        ax.set_title(f'{surface_name} Site Map (top view, slab only)')
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig(site_dir / "overlay.png", dpi=150, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {site_dir / 'overlay.png'}")
    except Exception as e:
        print(f"  (overlay visualization skipped: {e})")

    return surface, summary

def generate_CO(surface_name, surface, n_sub):
    """Step 2: CO* candidates with Cl[C]=O."""
    print(f"\n{'='*60}")
    print(f"[{surface_name}] Step 2: CO* candidate generation")
    print(f"{'='*60}")

    work_dir = G3_BASE / surface_name
    co_dir = work_dir / "CO"
    co_dir.mkdir(parents=True, exist_ok=True)

    print(f"  Fragment: Cl[C]=O ...")
    co = Fragment('Cl[C]=O', to_initialize=1)

    # Workaround: explicitly set smiles in conformer.info
    for conf in co.conformers:
        if 'smiles' not in conf.info:
            conf.info['smiles'] = 'Cl[C]=O'

    co_structs = surface.get_populated_sites(
        co, mode='all', sample_rotation=False,
        conformers_per_site_cap=1, overlap_thr=1.25
    )

    print(f"  CO* candidates generated: {len(co_structs)}")

    # Validate: first structure should have only C+O in adsorbate
    if len(co_structs) > 0:
        syms = co_structs[0].get_chemical_symbols()
        nonsubstrate = syms[n_sub:]
        print(f"  First structure adsorbate atoms: {nonsubstrate}")
        if 'H' in nonsubstrate:
            print(f"  ❌ FAIL: contains H. Aborting.")
            sys.exit(1)
        print(f"  ✅ OK: pure C+O fragment, no H")

    # Save
    write(co_dir / "candidates.traj", co_structs)

    summary = {
        "total": len(co_structs),
        "smiles": "Cl[C]=O",
        "fix_note": "guideline Cl[C-]#[O+] failed RDKit valence; corrected to Cl[C]=O (brackets prevent implicit H)",
        "fragment_atoms": dict(Counter(co_structs[0].get_chemical_symbols()[n_sub:])) if len(co_structs) > 0 else {},
    }

    with open(co_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Saved: {co_dir / 'candidates.traj'}")
    print(f"  Saved: {co_dir / 'summary.json'}")

    # Grid visualization
    if len(co_structs) > 0:
        from ase.visualize.plot import plot_atoms
        n_show = min(25, len(co_structs))
        idx = np.linspace(0, len(co_structs)-1, n_show, dtype=int)

        fig, axes = plt.subplots(5, 5, figsize=(15, 15))
        for i, ax in enumerate(axes.flat):
            if i < len(idx):
                plot_atoms(co_structs[idx[i]], ax, rotation='0x,0y,0z', radii=0.85, show_unit_cell=0)
                ax.set_title(f'#{idx[i]}', fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])

        plt.suptitle(f"{surface_name} CO* (Cl[C]=O) — {n_show} of {len(co_structs)}", fontsize=14, y=1.005)
        plt.tight_layout()
        plt.savefig(co_dir / "grid.png", dpi=120, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {co_dir / 'grid.png'}")

    return co_structs, summary

def generate_CH3O(surface_name, surface, n_sub):
    """Step 3: CH3O* candidates with ClOC."""
    print(f"\n{'='*60}")
    print(f"[{surface_name}] Step 3: CH3O* candidate generation")
    print(f"{'='*60}")

    work_dir = G3_BASE / surface_name
    ch_dir = work_dir / "CH3O"
    ch_dir.mkdir(parents=True, exist_ok=True)

    print(f"  Fragment: ClOC ...")
    ch3o = Fragment('ClOC', to_initialize=20)

    # Workaround: explicitly set smiles
    for conf in ch3o.conformers:
        if 'smiles' not in conf.info:
            conf.info['smiles'] = 'ClOC'

    ch_structs = surface.get_populated_sites(
        ch3o, mode='all', sample_rotation=True,
        conformers_per_site_cap=3, overlap_thr=1.25
    )

    print(f"  CH3O* candidates generated: {len(ch_structs)}")

    # Validate first structure
    if len(ch_structs) > 0:
        syms = ch_structs[0].get_chemical_symbols()
        nonsubstrate = syms[n_sub:]
        print(f"  First structure adsorbate atoms: {Counter(nonsubstrate)}")

    # Save
    write(ch_dir / "candidates.traj", ch_structs)

    summary = {
        "total": len(ch_structs),
        "smiles": "ClOC",
        "fragment_atoms": dict(Counter(ch_structs[0].get_chemical_symbols()[n_sub:])) if len(ch_structs) > 0 else {},
        "generation_params": {
            "to_initialize": 20,
            "sample_rotation": True,
            "conformers_per_site_cap": 3,
        }
    }

    with open(ch_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Saved: {ch_dir / 'candidates.traj'}")
    print(f"  Saved: {ch_dir / 'summary.json'}")

    # Grid visualization
    if len(ch_structs) > 0:
        from ase.visualize.plot import plot_atoms
        n_show = min(25, len(ch_structs))
        idx = np.linspace(0, len(ch_structs)-1, n_show, dtype=int)

        fig, axes = plt.subplots(5, 5, figsize=(15, 15))
        for i, ax in enumerate(axes.flat):
            if i < len(idx):
                plot_atoms(ch_structs[idx[i]], ax, rotation='0x,0y,0z', radii=0.85, show_unit_cell=0)
                ax.set_title(f'#{idx[i]}', fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])

        plt.suptitle(f"{surface_name} CH3O* (ClOC) — {n_show} of {len(ch_structs)}", fontsize=14, y=1.005)
        plt.tight_layout()
        plt.savefig(ch_dir / "grid.png", dpi=120, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {ch_dir / 'grid.png'}")

    return ch_structs, summary

def generate_coads(surface_name, co_structs, ch_structs, n_sub):
    """Step 4: Co-adsorption with updated guideline cutoffs."""
    print(f"\n{'='*60}")
    print(f"[{surface_name}] Step 4: Co-adsorption generation")
    print(f"{'='*60}")

    work_dir = G3_BASE / surface_name
    coads_dir = work_dir / "coads"
    coads_dir.mkdir(parents=True, exist_ok=True)

    # Stride sampling: 20x20 = 400 pairs
    co_stride = max(1, len(co_structs) // 20)
    ch_stride = max(1, len(ch_structs) // 20)
    co_idx = list(range(0, len(co_structs), co_stride))[:20]
    ch_idx = list(range(0, len(ch_structs), ch_stride))[:20]

    print(f"  Pair sampling: {len(co_idx)} CO × {len(ch_idx)} CH3O = {len(co_idx)*len(ch_idx)} pairs")

    setA, setB, setTS, side_path, rejected = [], [], [], [], []
    dists = []

    for ci in co_idx:
        a_co = co_structs[ci]
        pos_co = a_co.get_positions()
        iC = find_first_atom(a_co, 'C', n_sub)
        if iC is None:
            continue

        for chi in ch_idx:
            a_ch = ch_structs[chi]
            pos_ch = a_ch.get_positions()
            iO_ch = find_first_atom(a_ch, 'O', n_sub)
            if iO_ch is None:
                continue

            # Combine: CO slab + CH3O adsorbate fragment
            combined = a_co.copy()
            ch_adsorbate = a_ch[n_sub:].copy()
            combined += ch_adsorbate

            # Distance: C_CO ⋯ O_CH3O (with mic)
            d = combined.get_distance(iC, len(a_co) + (iO_ch - n_sub), mic=True)
            dists.append(d)

            # Steric check (only between CO* and CH3O* fragments)
            if check_steric_clash(combined, n_sub, len(a_co)):
                rejected.append(combined)
                continue

            # Updated guideline-compliant classification
            if d < CUTOFF_SIDE:
                side_path.append(combined)
            elif d < CUTOFF_TS_LOW:
                # Should not happen (CUTOFF_SIDE == CUTOFF_TS_LOW)
                pass
            elif d < CUTOFF_TS_HIGH:
                setTS.append(combined)
            elif d < CUTOFF_A_LOW:
                # Should not happen
                pass
            elif d < CUTOFF_A_HIGH:
                setA.append(combined)
            elif d < CUTOFF_AMB_HIGH:
                # Ambiguous range: skip
                pass
            else:  # d >= CUTOFF_AMB_HIGH (5.0 Å)
                setB.append(combined)

    print(f"  ✅ Classification complete:")
    print(f"    Set A (2.1–4.0 Å, reactive): {len(setA)}")
    print(f"    Set B (≥5.0 Å, thermo):      {len(setB)}")
    print(f"    Set TS (1.7–2.1 Å):          {len(setTS)}")
    print(f"    side-path (<1.7 Å):          {len(side_path)}")
    print(f"    rejected (steric):           {len(rejected)}")
    print(f"    distance range: {min(dists):.2f}–{max(dists):.2f} Å, mean {np.mean(dists):.2f}")

    # Save trajectories
    if setA:
        write(coads_dir / "SetA.traj", setA)
    if setB:
        write(coads_dir / "SetB.traj", setB)
    if setTS:
        write(coads_dir / "SetTS.traj", setTS)
    if side_path:
        write(coads_dir / "side.traj", side_path)
    if rejected:
        write(coads_dir / "rejected.traj", rejected)

    summary = {
        "total_pairs_generated": len(dists),
        "SetA_reactive_2.1-4.0": len(setA),
        "SetB_thermo_>=5.0": len(setB),
        "SetTS_1.7-2.1": len(setTS),
        "side_path_<1.7": len(side_path),
        "rejected_steric": len(rejected),
        "distance_min": float(min(dists)) if dists else 0.0,
        "distance_max": float(max(dists)) if dists else 0.0,
        "distance_mean": float(np.mean(dists)) if dists else 0.0,
        "cutoffs": {
            "side": CUTOFF_SIDE,
            "TS_low": CUTOFF_TS_LOW,
            "TS_high": CUTOFF_TS_HIGH,
            "A_low": CUTOFF_A_LOW,
            "A_high": CUTOFF_A_HIGH,
            "amb_high": CUTOFF_AMB_HIGH,
            "steric_heavy": STERIC_HEAVY,
            "steric_H": STERIC_H,
        },
        "co_smiles": "Cl[C]=O",
        "ch3o_smiles": "ClOC",
    }

    with open(coads_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Saved: {coads_dir / 'summary.json'}")

    # Distance histogram
    if dists:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(dists, bins=40, color='steelblue', edgecolor='black', alpha=0.7)

        # Overlay classification regions
        ax.axvspan(0, CUTOFF_SIDE, alpha=0.2, color='red', label=f'side-path (<{CUTOFF_SIDE})')
        ax.axvspan(CUTOFF_TS_LOW, CUTOFF_TS_HIGH, alpha=0.2, color='orange',
                   label=f'Set TS ({CUTOFF_TS_LOW}–{CUTOFF_TS_HIGH})')
        ax.axvspan(CUTOFF_A_LOW, CUTOFF_A_HIGH, alpha=0.3, color='green',
                   label=f'Set A reactive ({CUTOFF_A_LOW}–{CUTOFF_A_HIGH})')
        ax.axvspan(CUTOFF_AMB_HIGH, max(dists)*1.05, alpha=0.2, color='gray',
                   label=f'Set B thermo (≥{CUTOFF_AMB_HIGH})')

        ax.set_xlabel('d(C$_{CO}$ ⋯ O$_{CH_3O}$)  /  Å')
        ax.set_ylabel('count')
        ax.set_title(f'{surface_name} co-adsorption distance distribution')
        ax.legend(loc='upper right', fontsize=9)
        plt.tight_layout()
        plt.savefig(coads_dir / "distance_hist.png", dpi=140, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {coads_dir / 'distance_hist.png'}")

    # Grid visualizations for Set A and Set B
    from ase.visualize.plot import plot_atoms

    if setA:
        n_show = min(25, len(setA))
        idx_grid = np.linspace(0, len(setA)-1, n_show, dtype=int)
        fig, axes = plt.subplots(5, 5, figsize=(15, 15))
        for i, ax in enumerate(axes.flat):
            if i < len(idx_grid):
                plot_atoms(setA[idx_grid[i]], ax, rotation='0x,0y,0z', radii=0.85, show_unit_cell=0)
                ax.set_title(f'#{idx_grid[i]}', fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
        plt.suptitle(f"{surface_name} Set A reactive ({len(setA)} structures, C–O 2.1–4.0 Å)", fontsize=14)
        plt.tight_layout()
        plt.savefig(coads_dir / "SetA_grid.png", dpi=120, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {coads_dir / 'SetA_grid.png'}")

    if setB:
        n_show = min(25, len(setB))
        idx_grid = np.linspace(0, len(setB)-1, n_show, dtype=int)
        fig, axes = plt.subplots(5, 5, figsize=(15, 15))
        for i, ax in enumerate(axes.flat):
            if i < len(idx_grid):
                plot_atoms(setB[idx_grid[i]], ax, rotation='0x,0y,0z', radii=0.85, show_unit_cell=0)
                ax.set_title(f'#{idx_grid[i]}', fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
        plt.suptitle(f"{surface_name} Set B thermo ({len(setB)} structures, ≥5.0 Å)", fontsize=14)
        plt.tight_layout()
        plt.savefig(coads_dir / "SetB_grid.png", dpi=120, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {coads_dir / 'SetB_grid.png'}")

    return summary

def generate_report(surface_name, site_summary, co_summary, ch_summary, coads_summary):
    """Step 5: Generate markdown report."""
    print(f"\n{'='*60}")
    print(f"[{surface_name}] Step 5: Report generation")
    print(f"{'='*60}")

    work_dir = G3_BASE / surface_name
    report_dir = work_dir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    report = f"""# {surface_name} G3 Adsorption Sampling Report

Generated: {np.datetime64('now')}

## 1. Site Map

- Total sites detected: {site_summary['total_sites']}
- Slab atoms: {site_summary['slab_atoms']}
- Formula: {site_summary['slab_formula']}

### By site type:
```
{json.dumps(site_summary['site_types'], indent=2)}
```

## 2. CO* Candidates

- Total: {co_summary['total']}
- SMILES: `{co_summary['smiles']}`
- Fragment atoms: {co_summary['fragment_atoms']}
- Note: {co_summary['fix_note']}

## 3. CH3O* Candidates

- Total: {ch_summary['total']}
- SMILES: `{ch_summary['smiles']}`
- Fragment atoms: {ch_summary['fragment_atoms']}
- Generation params: {ch_summary['generation_params']}

## 4. Co-adsorption

### Classification (guideline-compliant cutoffs):

- **Set A (reactive, 2.1–4.0 Å)**: {coads_summary['SetA_reactive_2.1-4.0']}
- **Set B (thermo, ≥5.0 Å)**: {coads_summary['SetB_thermo_>=5.0']}
- **Set TS (1.7–2.1 Å)**: {coads_summary['SetTS_1.7-2.1']}
- **side-path (<1.7 Å)**: {coads_summary['side_path_<1.7']}
- **rejected (steric)**: {coads_summary['rejected_steric']}

### Distance statistics:

- Min: {coads_summary['distance_min']:.2f} Å
- Max: {coads_summary['distance_max']:.2f} Å
- Mean: {coads_summary['distance_mean']:.2f} Å

### Cutoffs applied:
```
{json.dumps(coads_summary['cutoffs'], indent=2)}
```

## 5. Output Files

- Site map: `site_map/sites.json`, `site_map/overlay.png`
- CO*: `CO/candidates.traj`, `CO/summary.json`, `CO/grid.png`
- CH3O*: `CH3O/candidates.traj`, `CH3O/summary.json`, `CH3O/grid.png`
- Co-ads: `coads/SetA.traj`, `coads/SetB.traj`, `coads/SetTS.traj`, `coads/side.traj`, `coads/rejected.traj`
- Co-ads plots: `coads/distance_hist.png`, `coads/SetA_grid.png`, `coads/SetB_grid.png`
- Report: `report/summary.md`

## 6. Notes

- Pipeline: AutoAdsorbate mode='all', overlap_thr=1.25
- CO* SMILES corrected from guideline to avoid RDKit valence error
- Steric rejection: heavy-heavy < {coads_summary['cutoffs']['steric_heavy']} Å OR H-heavy < {coads_summary['cutoffs']['steric_H']} Å
- Ambiguous range (4.0–5.0 Å) excluded from all sets

**Checkpoint stop: NO MLIP, NO DFT performed.**
"""

    with open(report_dir / "summary.md", "w") as f:
        f.write(report)

    print(f"  Saved: {report_dir / 'summary.md'}")

def main():
    """Main pipeline execution."""
    print("="*80)
    print("G3 Adsorption Sampling Pipeline: S2, S3b, S4")
    print("="*80)
    print(f"\nSurfaces to process: {len(SURFACES)}")
    for name, path, n in SURFACES:
        print(f"  - {name}: {n} substrate atoms")

    print(f"\nCutoffs (guideline-compliant):")
    print(f"  side-path: < {CUTOFF_SIDE} Å")
    print(f"  Set TS: {CUTOFF_TS_LOW} – {CUTOFF_TS_HIGH} Å")
    print(f"  Set A (reactive): {CUTOFF_A_LOW} – {CUTOFF_A_HIGH} Å")
    print(f"  ambiguous (reject): {CUTOFF_A_HIGH} – {CUTOFF_AMB_HIGH} Å")
    print(f"  Set B (thermo): ≥ {CUTOFF_AMB_HIGH} Å")
    print(f"  steric reject: heavy-heavy < {STERIC_HEAVY} Å OR H-heavy < {STERIC_H} Å")

    all_results = {}

    for surface_name, slab_path, n_sub in SURFACES:
        try:
            # Step 1: Site map
            surface, site_summary = generate_site_map(surface_name, slab_path, n_sub)

            # Step 2: CO*
            co_structs, co_summary = generate_CO(surface_name, surface, n_sub)

            # Step 3: CH3O*
            ch_structs, ch_summary = generate_CH3O(surface_name, surface, n_sub)

            # Step 4: Co-adsorption
            coads_summary = generate_coads(surface_name, co_structs, ch_structs, n_sub)

            # Step 5: Report
            generate_report(surface_name, site_summary, co_summary, ch_summary, coads_summary)

            all_results[surface_name] = {
                "site": site_summary,
                "CO": co_summary,
                "CH3O": ch_summary,
                "coads": coads_summary,
            }

        except Exception as e:
            print(f"\n❌ ERROR in {surface_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Final summary
    print("\n" + "="*80)
    print("PIPELINE COMPLETE — Summary")
    print("="*80)

    for surf_name, results in all_results.items():
        print(f"\n{surf_name}:")
        print(f"  Sites: {results['site']['total_sites']}")
        print(f"  CO*: {results['CO']['total']}")
        print(f"  CH3O*: {results['CH3O']['total']}")
        print(f"  Co-ads SetA: {results['coads']['SetA_reactive_2.1-4.0']}")
        print(f"  Co-ads SetB: {results['coads']['SetB_thermo_>=5.0']}")
        print(f"  Co-ads SetTS: {results['coads']['SetTS_1.7-2.1']}")
        print(f"  Co-ads side: {results['coads']['side_path_<1.7']}")
        print(f"  Co-ads rejected: {results['coads']['rejected_steric']}")
        print(f"  Distance range: {results['coads']['distance_min']:.2f}–{results['coads']['distance_max']:.2f} Å")

    print("\n✅ All surfaces processed. Checkpoint stop (NO MLIP, NO DFT).")
    print(f"\nOutput base: {G3_BASE}")

if __name__ == "__main__":
    main()
