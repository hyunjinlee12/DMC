"""Build G2 clean slabs from G1 optimized bulk structures.

T1.5: S1 Pd(100) — p(4×4), 5 layers
T1.6: S2 1 ML PdO(101)/Pd(100) — (√5×√5)R27° surface oxide (2×2 supercell)
T1.7: S3 PdO(100) — O-rich termination (2×2)
T1.8: S4 PdO₂(110) — stoichiometric

Usage:
    conda run -n pddmc python scripts/prepare_g2_slabs.py
"""

import numpy as np
from pathlib import Path
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Poscar
from ase import Atoms
from ase.io import write as ase_write
from ase.build import fcc100, make_supercell
from ase.constraints import FixAtoms

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from incar_templates import write_incar
from vasp_utils import write_potcar

PROJECT = Path(__file__).resolve().parent.parent
G1_DIR = PROJECT / "calculations" / "G1_bulk"
G2_DIR = PROJECT / "calculations" / "G2_slab"

VACUUM = 20.0  # Å


# ── helpers ──────────────────────────────────────────────────────────────

def read_bulk(material):
    return Structure.from_file(str(G1_DIR / material / "CONTCAR"))


def fix_bottom(atoms, n_layers=None, z_fraction=None, tol=0.5):
    """Fix bottom atoms. n_layers for clean metals, z_fraction for oxides."""
    z = atoms.positions[:, 2]

    if z_fraction is not None:
        z_min, z_max = z.min(), z.max()
        cutoff = z_min + z_fraction * (z_max - z_min)
        indices = [i for i, zi in enumerate(z) if zi <= cutoff]

    elif n_layers is not None:
        z_sorted = np.sort(z)
        layers = [[z_sorted[0]]]
        for zz in z_sorted[1:]:
            if zz - np.mean(layers[-1]) > tol:
                layers.append([zz])
            else:
                layers[-1].append(zz)
        centers = sorted(np.mean(L) for L in layers)
        n = min(n_layers, len(centers))
        cutoff = centers[n - 1] + tol
        indices = [i for i, zi in enumerate(z) if zi <= cutoff]
    else:
        indices = []

    atoms.set_constraint(FixAtoms(indices=indices))
    return atoms, indices


def add_vacuum_asym(atoms, vacuum=VACUUM):
    """Place slab near bottom, vacuum on top (asymmetric)."""
    pos = atoms.positions.copy()
    z_min = pos[:, 2].min()
    pos[:, 2] -= z_min - 1.5          # 1.5 Å buffer from cell floor
    z_max_new = pos[:, 2].max()
    new_c = z_max_new + vacuum - 1.5   # vacuum above

    cell = atoms.cell.copy()
    cell[2] = [0.0, 0.0, new_c]
    atoms.set_cell(cell)
    atoms.set_positions(pos)
    atoms.pbc = True
    return atoms


def write_inputs(atoms, outdir, template):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ase_write(str(outdir / "POSCAR"), atoms, format="vasp",
              vasp5=True, sort=True)
    write_incar(outdir / "INCAR", template, overrides={"KSPACING": 0.25})

    poscar = Poscar.from_file(str(outdir / "POSCAR"))
    seen = []
    for sp in poscar.site_symbols:
        s = str(sp)
        if s not in seen:
            seen.append(s)
    write_potcar(seen, outdir)


def info(atoms, label):
    syms = atoms.get_chemical_symbols()
    n_fixed = 0
    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            n_fixed += len(c.index)
    a = np.linalg.norm(atoms.cell[0])
    b = np.linalg.norm(atoms.cell[1])
    c = np.linalg.norm(atoms.cell[2])
    z = atoms.positions[:, 2]
    thick = z.max() - z.min()
    print(f"  {label}")
    print(f"    Atoms: {len(atoms)} (Pd={syms.count('Pd')}, O={syms.count('O')})")
    print(f"    Cell: {a:.2f} × {b:.2f} × {c:.2f} Å, thickness: {thick:.2f} Å")
    print(f"    Fixed: {n_fixed}/{len(atoms)} ({n_fixed/len(atoms)*100:.0f}%)")


# ── S1: Pd(100) p(4×4) 5-layer ──────────────────────────────────────────

def build_s1():
    print("\n=== S1: Pd(100) p(4×4) 5-layer ===")
    a = read_bulk("Pd").lattice.a  # 3.8907

    # ASE fcc100: p(4×4) with 5 layers, proper square cell
    # surface cell = a/√2 × a/√2 per atom → 4×(a/√2) ≈ 11.0 Å per side
    atoms = fcc100("Pd", size=(4, 4, 5), a=a, vacuum=0, periodic=True)
    atoms = add_vacuum_asym(atoms, VACUUM)
    atoms, _ = fix_bottom(atoms, n_layers=2, tol=0.5)

    info(atoms, "S1 Pd(100)")
    return atoms, "slab_metal"


# ── S2: 1 ML PdO(101)/Pd(100) (√5×√5)R27° 2×2 ─────────────────────────

def build_s2():
    """Pd(100)-(√5×√5)R27°-O surface oxide model.

    Ref: Todorova et al., Surf. Sci. 541, 101 (2003) / cond-mat/0304107.

    1) Pd(100) 4-layer substrate in (√5×√5)R27° × 2×2 cell
       (20 Pd/layer → 80 substrate Pd)
    2) PdO(101) oxide monolayer: 16 Pd_oxide + 16 O on top

    APPROXIMATE positions — VESTA verification required!
    """
    print("\n=== S2: 1 ML PdO(101)/Pd(100) (√5×√5)R27° 2×2 ===")
    a = read_bulk("Pd").lattice.a

    # Pd(100) p(1×1) 4-layer substrate (1 Pd per layer)
    sub = fcc100("Pd", size=(1, 1, 4), a=a, vacuum=0, periodic=True)

    # (√5×√5)R27° × 2×2 combined transformation
    # √5×√5: [[2,1],[-1,2]] (det=5, 5 atoms/layer)
    # × 2×2: [[2,0],[0,2]] (det=4)
    # Combined: [[4,2,0],[-2,4,0],[0,0,1]] (det=20, 20 atoms/layer)
    P = np.array([[4, 2, 0], [-2, 4, 0], [0, 0, 1]])
    sub = make_supercell(sub, P)
    sub = add_vacuum_asym(sub, VACUUM)

    print(f"  Substrate: {len(sub)} Pd, 4 layers")

    z_top = sub.positions[:, 2].max()
    a1 = sub.cell[0, :2]   # in-plane vec 1
    a2 = sub.cell[1, :2]   # in-plane vec 2

    # PdO(101) oxide layer — 16 PdO units in the 2×2 supercell
    # 4 Pd_oxide + 4 O per √5×√5 unit, tiled across 4 quadrants
    base_pd = [[0.05, 0.15], [0.15, 0.45], [0.25, 0.25], [0.45, 0.05]]
    base_o_lo = [[0.00, 0.30], [0.30, 0.00]]   # bridging O (low)
    base_o_hi = [[0.10, 0.40], [0.40, 0.10]]   # on-top O (high)

    quadrants = [[0, 0], [0.5, 0], [0, 0.5], [0.5, 0.5]]

    new_pos = list(sub.positions)
    new_sym = list(sub.get_chemical_symbols())

    for sx, sy in quadrants:
        for fx, fy in base_pd:
            xy = (fx + sx) * a1 + (fy + sy) * a2
            new_pos.append([xy[0], xy[1], z_top + 2.0])
            new_sym.append("Pd")
        for fx, fy in base_o_lo:
            xy = (fx + sx) * a1 + (fy + sy) * a2
            new_pos.append([xy[0], xy[1], z_top + 1.0])
            new_sym.append("O")
        for fx, fy in base_o_hi:
            xy = (fx + sx) * a1 + (fy + sy) * a2
            new_pos.append([xy[0], xy[1], z_top + 3.0])
            new_sym.append("O")

    atoms = Atoms(symbols=new_sym, positions=new_pos,
                  cell=sub.cell, pbc=True)
    atoms = add_vacuum_asym(atoms, VACUUM)
    atoms, _ = fix_bottom(atoms, n_layers=2, tol=0.5)

    info(atoms, "S2 PdO(101)/Pd(100)")
    print("    ⚠ Oxide layer APPROXIMATE — VESTA 검증 필수!")
    return atoms, "slab_oxide"


# ── S3: PdO(100) O-rich termination ─────────────────────────────────────

def build_s3():
    print("\n=== S3: PdO(100) O-rich termination ===")
    bulk = read_bulk("PdO")

    slabgen = SlabGenerator(
        bulk, (1, 0, 0),
        min_slab_size=10.0,
        min_vacuum_size=VACUUM,
        center_slab=False,
        lll_reduce=True,
        reorient_lattice=True,
    )
    slabs = slabgen.get_slabs(symmetrize=False)
    print(f"  Terminations: {len(slabs)}")

    best, best_score = None, -1
    for i, s in enumerate(slabs):
        atmp = AseAtomsAdaptor.get_atoms(s)
        z = atmp.positions[:, 2]
        top_mask = z > (z.max() - 1.5)
        top_syms = [atmp[j].symbol for j in np.where(top_mask)[0]]
        o_frac = top_syms.count("O") / max(len(top_syms), 1)
        print(f"    Term {i}: {len(atmp)} atoms, top={top_syms}, O%={o_frac:.0%}")
        if o_frac > best_score:
            best_score, best = o_frac, s

    best.make_supercell([2, 2, 1])
    atoms = AseAtomsAdaptor.get_atoms(best)
    atoms = add_vacuum_asym(atoms, VACUUM)
    atoms, _ = fix_bottom(atoms, z_fraction=0.35)

    info(atoms, "S3 PdO(100)")
    return atoms, "slab_oxide"


# ── S4: PdO₂(110) stoichiometric ────────────────────────────────────────

def build_s4():
    print("\n=== S4: PdO₂(110) stoichiometric ===")
    bulk = read_bulk("PdO2")

    slabgen = SlabGenerator(
        bulk, (1, 1, 0),
        min_slab_size=10.0,
        min_vacuum_size=VACUUM,
        center_slab=False,
        lll_reduce=True,
        reorient_lattice=True,
    )
    slabs = slabgen.get_slabs(symmetrize=False)
    print(f"  Terminations: {len(slabs)}")

    best, best_dev = None, 999
    for i, s in enumerate(slabs):
        comp = s.composition.as_dict()
        n_pd, n_o = comp.get("Pd", 0), comp.get("O", 0)
        ratio = n_o / max(n_pd, 1)
        dev = abs(ratio - 2.0)
        tag = "← stoich" if dev < 0.15 else ""
        print(f"    Term {i}: {len(s)} atoms, Pd={n_pd:.0f} O={n_o:.0f}, "
              f"O/Pd={ratio:.2f} {tag}")
        if dev < best_dev:
            best_dev, best = dev, s

    # Supercell to get in-plane > ~8 Å
    atmp = AseAtomsAdaptor.get_atoms(best)
    a_len = np.linalg.norm(atmp.cell[0])
    b_len = np.linalg.norm(atmp.cell[1])
    nx = max(1, int(np.ceil(8.0 / a_len)))
    ny = max(1, int(np.ceil(8.0 / b_len)))
    if nx > 1 or ny > 1:
        best.make_supercell([nx, ny, 1])
        print(f"  Supercell: {nx}×{ny}×1")

    atoms = AseAtomsAdaptor.get_atoms(best)
    atoms = add_vacuum_asym(atoms, VACUUM)
    atoms, _ = fix_bottom(atoms, z_fraction=0.35)

    info(atoms, "S4 PdO₂(110)")
    return atoms, "slab_oxide"


# ── main ─────────────────────────────────────────────────────────────────

def main():
    G2_DIR.mkdir(parents=True, exist_ok=True)

    surfaces = [
        ("S1_Pd100",         build_s1),
        ("S2_PdO101_Pd100",  build_s2),
        ("S3_PdO100",        build_s3),
        ("S4_PdO2_110",      build_s4),
    ]

    for name, builder in surfaces:
        try:
            atoms, template = builder()
            write_inputs(atoms, G2_DIR / name, template)
            print(f"  → {name}/ INCAR({template}, KSPACING=0.25)\n")
        except Exception as e:
            print(f"  ERROR {name}: {e}")
            import traceback; traceback.print_exc()

    print("=" * 60)
    print("G2 Slab Summary")
    print("=" * 60)
    for d in sorted(G2_DIR.iterdir()):
        if d.is_dir() and (d / "POSCAR").exists():
            s = Structure.from_file(str(d / "POSCAR"))
            print(f"  {d.name:25s} {s.num_sites:4d} atoms  "
                  f"{s.lattice.a:7.2f} × {s.lattice.b:6.2f} × {s.lattice.c:6.2f} Å")

    print(f"\nOutput: {G2_DIR}")
    print("VESTA로 모든 POSCAR 육안 확인 후 DFT 제출!")


if __name__ == "__main__":
    main()
