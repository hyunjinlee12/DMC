"""Render S4 CO idx=7 (CO + O_lat → CO2*) for comparison."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ase.io import read
from render_structure import render

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab/S4_PdO2_110'
G3 = ROOT / 'calculations/G3_adsorption/S4_PdO2_110'
OUT = ROOT / 'reports/predft_advisor_figures/top1_adsorbate_sideviews_v2'

traj = list(read(G3/'MLIP_phase1/relaxed_CO.traj', index=':'))
slab = read(G2/'CONTCAR')
atoms = traj[7]
if len(atoms) != len(slab)+2:
    ads = atoms[-2:]; atoms = slab.copy(); atoms += ads
atoms_22 = atoms * (2, 2, 1)
render(atoms_22, OUT/'S4_CO_idx7_CO2like.png',
       rotation='-90x,-90y,0z', width=2400, show_cell=False)
print('S4 CO idx=7 ✓')
