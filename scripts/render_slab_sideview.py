"""Pure side view (orthographic) of converged G2 slabs — no perspective."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read
from ase.visualize.plot import plot_atoms

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']

G2 = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab')
OUT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/predft_advisor_figures/slab_sideviews')
OUT.mkdir(parents=True, exist_ok=True)

SDIR = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
        'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

# Individual high-DPI side views
for sid, sdir in SDIR.items():
    atoms = read(G2 / sdir / 'CONTCAR')
    fig, ax = plt.subplots(figsize=(6, 5))
    plot_atoms(atoms, ax, rotation='-90x', radii=0.55, show_unit_cell=2)
    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(OUT/f'{sid}_side.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{sid} ✓')

# Combined 1x5 panel
fig, axes = plt.subplots(1, 5, figsize=(20, 5))
for i, (sid, sdir) in enumerate(SDIR.items()):
    atoms = read(G2 / sdir / 'CONTCAR')
    plot_atoms(atoms, axes[i], rotation='-90x', radii=0.55, show_unit_cell=2)
    axes[i].set_axis_off()
    axes[i].set_title(sid, fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(OUT/'ALL_side.png', dpi=300, bbox_inches='tight')
plt.close()
print('ALL ✓')
