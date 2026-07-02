"""S5 Pd(111) — built EXACTLY like S1 Pd(100), only the surface face changed (100 → 111).

Uses same:
  - bulk a from G1/Pd/CONTCAR (DFT-relaxed Pd)
  - size=(4, 4, 5) → 80 atoms
  - vacuum=0 in builder + add_vacuum_asym (20 Å asymmetric)
  - fix_bottom(n_layers=2)
  - hexagonal natural cell (orthogonal=False, like fcc100 default)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import write as ase_write
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

PROJECT = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc")
G1_DIR = PROJECT / "calculations/G1_bulk"
G2_DIR = PROJECT / "calculations/G2_slab"
VACUUM = 20.0


def fix_bottom(atoms, n_layers=2, tol=0.5):
    z = atoms.positions[:, 2]
    z_sorted = np.sort(z)
    layers = [[z_sorted[0]]]
    for zz in z_sorted[1:]:
        if zz - np.mean(layers[-1]) > tol:
            layers.append([zz])
        else:
            layers[-1].append(zz)
    centers = sorted(np.mean(L) for L in layers)
    cutoff = centers[n_layers - 1] + tol
    indices = [i for i, zi in enumerate(z) if zi <= cutoff]
    atoms.set_constraint(FixAtoms(indices=indices))
    return atoms


def add_vacuum_asym(atoms, vacuum=VACUUM):
    pos = atoms.positions.copy()
    z_min = pos[:, 2].min()
    pos[:, 2] -= z_min - 1.5
    z_max_new = pos[:, 2].max()
    new_c = z_max_new + vacuum - 1.5
    cell = atoms.cell.copy()
    cell[2] = [0.0, 0.0, new_c]
    atoms.set_cell(cell)
    atoms.set_positions(pos)
    atoms.pbc = True
    return atoms


# === Build S5 Pd(111) — same recipe as S1 Pd(100), face only changed ===
a = Structure.from_file(str(G1_DIR / "Pd/CONTCAR")).lattice.a
print(f"Pd bulk a (from G1/Pd/CONTCAR) = {a:.4f} Å")

atoms = fcc111("Pd", size=(4, 4, 5), a=a, vacuum=0, periodic=True)
atoms = add_vacuum_asym(atoms, VACUUM)
atoms = fix_bottom(atoms, n_layers=2)

print(f"S5 Pd(111): {len(atoms)} atoms, cell={atoms.cell.lengths()}, "
      f"angles={atoms.cell.angles()}")
print(f"  z range: {atoms.positions[:,2].min():.2f} - {atoms.positions[:,2].max():.2f}, "
      f"vacuum_gap={atoms.cell.lengths()[2]-(atoms.positions[:,2].max()-atoms.positions[:,2].min()):.2f} Å")

# Save (same outdir, overwrite POSCAR)
out = G2_DIR / "S5_Pd111"
out.mkdir(exist_ok=True)
ase_write(str(out / "POSCAR"), atoms, format="vasp", vasp5=True, sort=True)
print(f"Saved POSCAR → {out/'POSCAR'}")
