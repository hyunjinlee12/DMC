"""Diagnostic script to understand steric rejection."""
import numpy as np
from ase.io import read
from pathlib import Path

# Load one rejected structure from each surface
surfaces = [
    ("S2_PdO101_Pd100", 112),
    ("S3b_PdO100_PdOterm", 104),
    ("S4_PdO2_110", 144),
]

for surf_name, n_sub in surfaces:
    print(f"\n{'='*60}")
    print(f"{surf_name}")
    print(f"{'='*60}")

    traj_path = Path(f"/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G3_adsorption/{surf_name}/coads/rejected.traj")

    if not traj_path.exists():
        print("  No rejected.traj found")
        continue

    structs = list(read(traj_path, index=':'))
    print(f"  Total rejected: {len(structs)}")

    # Analyze first structure
    atoms = structs[0]
    syms = atoms.get_chemical_symbols()
    pos = atoms.get_positions()

    # Adsorbate region
    ads_indices = list(range(n_sub, len(atoms)))
    ads_syms = [syms[i] for i in ads_indices]

    print(f"  Adsorbate atoms: {len(ads_indices)}, symbols: {set(ads_syms)}")

    # Find all adsorbate-adsorbate distances
    min_dists = []
    for i in ads_indices:
        for j in ads_indices:
            if i >= j:
                continue
            d = atoms.get_distance(i, j, mic=True)
            min_dists.append((d, syms[i], syms[j]))

    # Sort by distance
    min_dists.sort()

    print(f"  Adsorbate-adsorbate distances (first 10):")
    for d, si, sj in min_dists[:10]:
        flag = ""
        # Check against cutoffs
        if si != 'H' and sj != 'H':  # heavy-heavy
            if d < 1.6:
                flag = " ❌ REJECT (heavy-heavy < 1.6)"
        elif (si == 'H' and sj != 'H') or (si != 'H' and sj == 'H'):  # H-heavy
            if d < 1.1:
                flag = " ❌ REJECT (H-heavy < 1.1)"

        print(f"    {d:5.2f} Å  {si:2s}-{sj:2s}{flag}")

    # Count rejections by type
    heavy_heavy_reject = sum(1 for d, si, sj in min_dists
                             if si != 'H' and sj != 'H' and d < 1.6)
    h_heavy_reject = sum(1 for d, si, sj in min_dists
                         if ((si == 'H' and sj != 'H') or (si != 'H' and sj == 'H')) and d < 1.1)

    print(f"  Heavy-heavy violations (< 1.6 Å): {heavy_heavy_reject}")
    print(f"  H-heavy violations (< 1.1 Å): {h_heavy_reject}")

    # Also check substrate-adsorbate distances (might be issue)
    print(f"\n  Checking substrate-adsorbate distances:")
    sub_ads_min = []
    for i in range(n_sub):
        for j in ads_indices:
            d = atoms.get_distance(i, j, mic=True)
            sub_ads_min.append((d, syms[i], syms[j]))

    sub_ads_min.sort()
    print(f"  Substrate-adsorbate closest 10:")
    for d, si, sj in sub_ads_min[:10]:
        print(f"    {d:5.2f} Å  {si:2s}(sub)-{sj:2s}(ads)")

print("\n" + "="*60)
print("Diagnosis complete")
