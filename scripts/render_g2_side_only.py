"""Use render_structure.render() with side rotation — same as original side.png files."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ase.io import read
from render_structure import render

G2 = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab")
OUT = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/predft_advisor_figures/slab_sideviews_final")
OUT.mkdir(parents=True, exist_ok=True)

for sid, sdir in [('S1','S1_Pd100'),('S2','S2_PdO101_Pd100'),
                  ('S3','S3_PdO100'),('S3b','S3b_PdO100_PdOterm'),
                  ('S4','S4_PdO2_110')]:
    atoms = read(G2 / sdir / 'CONTCAR') * (2, 2, 1)
    out = OUT / f'{sid}_side.png'
    r = render(atoms, out, rotation="-90x,-90y,0z", width=2400, show_cell=False)
    print(f'{sid} {"OK" if r else "FAIL"}')
