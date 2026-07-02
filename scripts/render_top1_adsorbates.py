"""Top-1 CO* and CH3O* per surface, bc-plane side view, 1x5 layout."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ase.io import read
from render_structure import render
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial','Liberation Sans','DejaVu Sans']

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
G2 = ROOT / 'calculations/G2_slab'
OUT = ROOT / 'reports/predft_advisor_figures/top1_adsorbate_sideviews'
OUT.mkdir(parents=True, exist_ok=True)

SURF = [('S1','S1_Pd100','Pd(100)','Pd$^0$'),
        ('S2','S2_PdO101_Pd100','1ML PdO(101)/Pd(100)','Pd$^0$+Pd$^{2+}$'),
        ('S3b','S3b_PdO100_PdOterm','Pd-rich PdO(100)','Pd$^{2+}$ Pd-top'),
        ('S3','S3_PdO100','O-rich PdO(100)','Pd$^{2+}$ O-top'),
        ('S4','S4_PdO2_110','PdO$_2$(110)','Pd$^{4+}$')]


def top1_atoms(sid, sdir, ads):
    n_ads = 2 if ads == 'CO' else 5
    unique = json.load(open(G3 / sdir / f'MLIP_phase1/unique_{ads}.json'))
    traj = list(read(G3 / sdir / f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
    top = sorted(unique, key=lambda r: r['E'])[0]
    atoms = traj[top['idx']]
    slab = read(G2 / sdir / 'CONTCAR')
    if len(atoms) != len(slab) + n_ads:
        ads_atoms = atoms[-n_ads:]
        atoms = slab.copy(); atoms += ads_atoms
    return atoms * (2, 2, 1)


def make_row(ads, label):
    pngs = []
    for sid, sdir, name, ox in SURF:
        atoms = top1_atoms(sid, sdir, ads)
        png = OUT / f'{sid}_{ads}_top1.png'
        render(atoms, png, rotation='-90x,-90y,0z', width=1800, show_cell=False)
        pngs.append((sid, name, ox, png))

    fig, axes = plt.subplots(1, 5, figsize=(22, 5.5))
    for ax, (sid, name, ox, png) in zip(axes, pngs):
        ax.imshow(mpimg.imread(png)); ax.set_axis_off()
        ax.set_title(f'{sid}: {name}\n{ox}', fontsize=12, fontweight='bold')
    fig.suptitle(f'Top-1 {label} (MACE+D3, lowest E)', fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(OUT/f'ALL_{ads}_top1.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{ads} row ✓')

make_row('CO', r'CO$^*$')
make_row('CH3O', r'CH$_3$O$^*$')
