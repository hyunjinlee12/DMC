"""Build 2 additional surface models for advisor PPT:
   S5 = Pd(111)        — close-packed metal
   S6 = PdO2(110)-Ocus — S4 + extra O on cus-Pd sites
"""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ase.io import read, write
from ase.build import fcc111
from ase.constraints import FixAtoms
from render_structure import render

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab'
OUT_FIG = ROOT / 'reports/predft_advisor_figures/extra_surfaces'
OUT_FIG.mkdir(parents=True, exist_ok=True)

# ===== S5: Pd(111) =====
# Lattice constant from S1 (G2 relaxed Pd) — use NN distance
s1 = read(G2/'S1_Pd100/CONTCAR')
d = s1.get_all_distances(mic=True)
nn = sorted(d[(d>0.1)&(d<4)].flatten())[0]
a_bulk = float(nn * np.sqrt(2))
print(f'Pd NN = {nn:.3f} Å → bulk a = {a_bulk:.4f} Å')

pd111 = fcc111('Pd', size=(4, 4, 5), a=a_bulk, vacuum=10.0, orthogonal=False)  # vacuum=10 → 20 Å gap (matches S1-S4)
n_pd = len(pd111)
# fix bottom 2 layers
z = pd111.positions[:, 2]
fixed = [i for i in range(n_pd) if z[i] < np.median(z)]
pd111.set_constraint(FixAtoms(indices=fixed))
S5_dir = G2/'S5_Pd111'; S5_dir.mkdir(exist_ok=True)
write(str(S5_dir/'POSCAR'), pd111, format='vasp', sort=True, direct=True, vasp5=True)
print(f'S5 Pd(111): {len(pd111)} atoms, cell {pd111.cell.lengths()[:2]}')

# ===== S6: PdO2(110)-Ocus =====
# True cus-Pd = top-layer Pd with only 5 O neighbors (missing apical)
# Pd_full = already 6-coordinated, should NOT add more O
s4 = read(G2/'S4_PdO2_110/CONTCAR')
syms = s4.get_chemical_symbols()
pd_idx = [i for i,s in enumerate(syms) if s=='Pd']
o_idx = [i for i,s in enumerate(syms) if s=='O']
z_arr = s4.positions[:,2]
z_max = max(z_arr[i] for i in pd_idx)
top_pd = [i for i in pd_idx if z_arr[i] > z_max - 0.5]
# Identify cus-Pd: n_O < 6 within 2.5 Å (typical 5-coord)
cus_pd = []
for i in top_pd:
    n_o = sum(1 for j in o_idx if s4.get_distance(i, j, mic=True) < 2.5)
    if n_o == 5:
        cus_pd.append(i)
print(f'S4 top Pd: {len(top_pd)}, Pd_cus (5-coord): {len(cus_pd)}, Pd_full (6-coord): {len(top_pd)-len(cus_pd)}')
# Add O 1.85 Å above each cus-Pd only
from ase import Atoms
new_atoms = s4.copy()
for i in cus_pd:
    o_pos = s4.positions[i] + np.array([0, 0, 1.85])
    new_atoms += Atoms('O', positions=[o_pos])
S6_dir = G2/'S6_PdO2_110_Ocus'; S6_dir.mkdir(exist_ok=True)
write(str(S6_dir/'POSCAR'), new_atoms, format='vasp', sort=True, direct=True, vasp5=True)
print(f'S6 PdO2(110)-Ocus: {len(new_atoms)} atoms (+{len(cus_pd)} Ocus on Pd_cus only)')

# ===== Render side views =====
SURF = [('S5', S5_dir/'POSCAR', 'Pd(111), Pd⁰'),
        ('S6', S6_dir/'POSCAR', 'PdO₂(110)-Ocus, Pd⁴⁺')]
for tag, pos_f, lbl in SURF:
    atoms = read(pos_f) * (2,2,1)
    out = OUT_FIG / f'{tag}_side.png'
    render(atoms, out, rotation='-90x,-90y,0z', width=2400, show_cell=False)
    print(f'{tag} side view ✓')
