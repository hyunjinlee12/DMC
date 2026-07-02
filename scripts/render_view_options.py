"""Render 4 candidate side-view rotations for S2 (most complex cell) — pick which one matches the user's desired view (a out, b right, c up)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ase.io import read
from render_structure import render

G2 = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab')
OUT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/predft_advisor_figures/view_options')
OUT.mkdir(parents=True, exist_ok=True)

atoms = read(G2/'S2_PdO101_Pd100/CONTCAR')

OPTIONS = [
    ('A', '-90x,0y,0z'),          # current: a right, b depth, c up
    ('B', '-90x,0y,-90z'),        # rotate frame so b right
    ('C', '-90x,0y,90z'),
    ('D', '90y,-90x,0z'),
]
for tag, rot in OPTIONS:
    render(atoms, OUT/f'opt_{tag}_{rot.replace(",","_")}.png',
           rotation=rot, width=1200, show_cell=True)
    print(f'{tag} ({rot}) ✓')
