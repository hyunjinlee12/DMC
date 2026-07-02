"""POSCAR (initial) vs CONTCAR (converged) side-by-side, side view."""
import sys, subprocess
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ase.io import read
from render_structure import render
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

G2 = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab')
OUT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/predft_advisor_figures/slab_init_vs_relaxed')
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial','Liberation Sans','DejaVu Sans']

for sid, sdir in [('S1','S1_Pd100'),('S2','S2_PdO101_Pd100'),
                  ('S3','S3_PdO100'),('S3b','S3b_PdO100_PdOterm'),
                  ('S4','S4_PdO2_110')]:
    p_init = OUT / f'{sid}_init.png'
    p_conv = OUT / f'{sid}_conv.png'
    render(read(G2/sdir/'POSCAR'), p_init, rotation='-90x,0y,0z', width=1600, show_cell=False)
    render(read(G2/sdir/'CONTCAR'), p_conv, rotation='-90x,0y,0z', width=1600, show_cell=False)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, png, title in [(axes[0], p_init, 'POSCAR (initial)'),
                           (axes[1], p_conv, 'CONTCAR (converged)')]:
        ax.imshow(mpimg.imread(png)); ax.set_axis_off()
        ax.set_title(title, fontsize=13, fontweight='bold')
    fig.suptitle(sid, fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUT/f'{sid}_compare.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{sid} ✓')
