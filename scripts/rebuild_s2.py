"""Rebuild S2: PdO(101)/Pd(100) surface oxide — v4 (fixed wrapping bug).

Uses vacuum > 0 when building substrate to prevent layer wrapping.
Extracts PdO(101) monolayer from bulk PdO via rotation, preserving Pd-O = 2.039 Å.

Usage: conda run -n pddmc python /tmp/rebuild_s2_v4.py
"""
import numpy as np
from pathlib import Path
from ase.build import fcc100, make_supercell
from ase import Atoms
from ase.io import write as ase_write, read
from ase.constraints import FixAtoms
from ase.neighborlist import neighbor_list
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator

PROJECT = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc")
G1 = PROJECT / "calculations" / "G1_bulk"
OUTDIR = PROJECT / "calculations" / "G2_slab" / "S2_PdO101_Pd100"

VACUUM = 20.0
INTERFACE_GAP = 2.2

def main():
    pd_struct = Structure.from_file(str(G1 / "Pd" / "CONTCAR"))
    a_pd = pd_struct.lattice.a
    d_nn = a_pd / np.sqrt(2)
    d_100 = a_pd / 2  # (100) interlayer spacing

    pdo_struct = Structure.from_file(str(G1 / "PdO" / "CONTCAR"))

    print(f"Pd: a={a_pd:.4f}, d_nn={d_nn:.4f}, d_100={d_100:.4f}")
    print(f"PdO: a={pdo_struct.lattice.a:.4f}, c={pdo_struct.lattice.c:.4f}")

    # ── PdO(101) monolayer from bulk ──
    slabgen = SlabGenerator(pdo_struct, (1, 0, 1), min_slab_size=6,
                            min_vacuum_size=12, center_slab=False)
    pdo_slab = slabgen.get_slabs()[0]

    a1_s = pdo_slab.lattice.matrix[0]
    a2_s = pdo_slab.lattice.matrix[1]
    n_vec = np.cross(a1_s, a2_s)
    n_hat = n_vec / np.linalg.norm(n_vec)

    # Extract bottom 4 atoms (1 repeat: 2 Pd + 2 O)
    all_h = [np.dot(site.coords, n_hat) for site in pdo_slab]
    sorted_idx = np.argsort(all_h)[:4]

    mono_species = [str(pdo_slab[i].species_string) for i in sorted_idx]
    mono_coords = [pdo_slab[i].coords.copy() for i in sorted_idx]

    # Rotate monolayer so surface normal → z-axis
    cos_a, sin_a = n_hat[2], n_hat[1]
    R = np.array([[1, 0, 0], [0, cos_a, -sin_a], [0, sin_a, cos_a]])

    a1_flat = R @ a1_s  # (3.054, 0, 0)
    a2_flat = R @ a2_s  # (0, 6.209, 0)

    mono_flat = [R @ c for c in mono_coords]
    z_min_mono = min(c[2] for c in mono_flat)
    for c in mono_flat:
        c[2] -= z_min_mono

    # Convert to fractional coords of flat PdO(101) cell
    cell_2d = np.array([[a1_flat[0], a1_flat[1]], [a2_flat[0], a2_flat[1]]])
    cell_inv = np.linalg.inv(cell_2d)

    mono_frac = []
    for sp, c in zip(mono_species, mono_flat):
        fxy = (cell_inv @ c[:2]) % 1.0
        mono_frac.append((sp, fxy[0], fxy[1], c[2]))

    # 2×1 supercell of monolayer
    mono_2x1 = []
    for sp, fa, fb, h in mono_frac:
        mono_2x1.append((sp, fa / 2, fb, h))
        mono_2x1.append((sp, fa / 2 + 0.5, fb, h))

    print(f"\nMonolayer 2×1: {len(mono_2x1)} atoms")
    for sp, fa, fb, h in mono_2x1:
        print(f"  {sp:3s} frac=({fa:.4f}, {fb:.4f}) z_rel={h:.3f}")

    # ── Pd(100) √5×√5 substrate (WITH vacuum to avoid wrapping) ──
    sub_prim = fcc100("Pd", size=(1, 1, 4), a=a_pd, vacuum=15.0, periodic=True)
    P_r5 = np.array([[2, 1, 0], [-1, 2, 0], [0, 0, 1]])
    sub = make_supercell(sub_prim, P_r5)

    # Verify 4 distinct layers
    z_sub = sub.positions[:, 2]
    z_unique_sub = sorted(set(round(zi, 3) for zi in z_sub))
    print(f"\nSubstrate: {len(sub)} atoms")
    print(f"  Layers at z = {z_unique_sub}")
    assert len(z_unique_sub) == 4, f"Expected 4 layers, got {len(z_unique_sub)}"

    # Shift so bottom layer at z=1.5
    sub.positions[:, 2] += 1.5 - z_sub.min()
    z_top_sub = sub.positions[:, 2].max()

    z_layers = sorted(set(round(zi, 3) for zi in sub.positions[:, 2]))
    print(f"  After shift: layers at {z_layers}")
    print(f"  z_top = {z_top_sub:.3f}")

    # ── Place oxide on substrate ──
    a1_cell = sub.cell[0, :2]
    a2_cell = sub.cell[1, :2]

    all_sym = list(sub.get_chemical_symbols())
    all_pos = list(sub.positions)

    z_base = z_top_sub + INTERFACE_GAP

    for sp, fa, fb, h in mono_2x1:
        xy = fa * a1_cell + fb * a2_cell
        all_pos.append([xy[0], xy[1], z_base + h])
        all_sym.append(sp)

    z_max = max(p[2] for p in all_pos)
    cell_new = sub.cell.copy()
    cell_new[2] = [0, 0, z_max + VACUUM]

    combined = Atoms(symbols=all_sym, positions=all_pos, cell=cell_new, pbc=True)

    # ── 2×2 supercell ──
    P_2x2 = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]])
    final = make_supercell(combined, P_2x2)

    # ── Fix bottom 2 layers ──
    z = final.positions[:, 2]
    z_all_unique = sorted(set(round(zi, 3) for zi in z))
    layer_groups = [[z_all_unique[0]]]
    for zv in z_all_unique[1:]:
        if zv - np.mean(layer_groups[-1]) > 0.5:
            layer_groups.append([zv])
        else:
            layer_groups[-1].append(zv)

    centers = sorted(np.mean(L) for L in layer_groups)
    cutoff_z = centers[1] + 0.5
    fix_idx = [i for i, zi in enumerate(z) if zi <= cutoff_z]
    final.set_constraint(FixAtoms(indices=fix_idx))

    # ── Validation ──
    syms = final.get_chemical_symbols()
    n_pd = syms.count("Pd")
    n_o = syms.count("O")

    print(f"\n{'='*50}")
    print(f"FINAL S2 STRUCTURE")
    print(f"{'='*50}")
    print(f"Atoms: {len(final)} (Pd={n_pd}, O={n_o})")
    print(f"Cell: {np.linalg.norm(final.cell[0]):.3f} × "
          f"{np.linalg.norm(final.cell[1]):.3f} × "
          f"{np.linalg.norm(final.cell[2]):.3f}")
    print(f"Fixed: {len(fix_idx)} atoms (bottom 2 layers)")

    # Layer analysis
    print(f"\nLayers ({len(layer_groups)}):")
    for lg in layer_groups:
        z_center = np.mean(lg)
        mask = [i for i in range(len(final))
                if any(abs(z[i] - zv) < 0.15 for zv in lg)]
        pd_n = sum(1 for i in mask if syms[i] == "Pd")
        o_n = sum(1 for i in mask if syms[i] == "O")
        comp = []
        if pd_n: comp.append(f"Pd={pd_n}")
        if o_n: comp.append(f"O={o_n}")
        fixed = sum(1 for i in mask if i in fix_idx)
        fix_tag = " [FIXED]" if fixed else ""
        print(f"  z≈{z_center:6.2f}: {len(mask):3d} atoms  {' '.join(comp)}{fix_tag}")

    # Pd-O distances
    ii, jj, dd = neighbor_list('ijd', final, cutoff=3.0)
    pd_o_dists = []
    for a, b, d in zip(ii, jj, dd):
        if {syms[a], syms[b]} == {"Pd", "O"}:
            pd_o_dists.append(d)

    dists = np.array(pd_o_dists) if pd_o_dists else np.array([])
    if len(dists) > 0:
        print(f"\nPd-O distances: min={dists.min():.3f}, max={dists.max():.3f}, "
              f"mean={dists.mean():.3f}")
        short = dists[dists < 1.5]
        if len(short):
            print(f"  *** WARNING: {len(short)//2} pairs < 1.5 Å ***")
        else:
            print(f"  ✓ No short contacts")

        for lo, hi in [(1.5, 2.0), (2.0, 2.2), (2.2, 2.5), (2.5, 3.0)]:
            count = np.sum((dists >= lo) & (dists < hi))
            print(f"    [{lo:.1f}-{hi:.1f}): {count//2} pairs")

    # Substrate-oxide interface check
    sub_pd_z = [z[i] for i in range(len(final)) if syms[i] == "Pd" and z[i] < z_base]
    ox_o_z = [z[i] for i in range(len(final)) if syms[i] == "O"]
    if sub_pd_z and ox_o_z:
        z_gap = min(ox_o_z) - max(sub_pd_z)
        print(f"\nInterface gap (top Pd → bottom O): {z_gap:.3f} Å")

    # ── Save ──
    OUTDIR.mkdir(parents=True, exist_ok=True)
    poscar_path = OUTDIR / "POSCAR"
    ase_write(str(poscar_path), final, format="vasp", vasp5=True, sort=True)
    print(f"\n→ Saved: {poscar_path}")

    check = read(str(poscar_path))
    z_check = check.positions[:, 2]
    z_check_unique = sorted(set(round(zi, 3) for zi in z_check))
    print(f"  Verified: {len(check)} atoms, {len(z_check_unique)} z-levels")


if __name__ == "__main__":
    main()
