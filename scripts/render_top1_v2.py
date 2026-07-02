"""Top-1 VALID adsorbates + custom colors (C dark, H white) + larger H radii."""
import sys, json, re, subprocess
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ase.io import read
from ase.io.pov import write_pov
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial','Liberation Sans','DejaVu Sans']

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
G2 = ROOT / 'calculations/G2_slab'
OUT = ROOT / 'reports/predft_advisor_figures/top1_adsorbate_sideviews_v2'
OUT.mkdir(parents=True, exist_ok=True)
POVRAY = "/home/hyunjin/.conda/envs/pddmc/bin/povray"

COLORS = {
    'Pd': [0.165, 0.616, 0.561],
    'O':  [0.941, 0.443, 0.404],
    'C':  [0.20, 0.20, 0.20],    # dark gray
    'H':  [0.95, 0.95, 0.95],    # near-white
}
RADII = {'Pd': 0.85, 'O': 0.65, 'C': 0.55, 'H': 0.45}   # H bigger than typical

SURF = [('S1','S1_Pd100','Pd(100)','Pd$^0$'),
        ('S2','S2_PdO101_Pd100','1ML PdO(101)/Pd(100)','Pd$^0$+Pd$^{2+}$'),
        ('S3b','S3b_PdO100_PdOterm','Pd-rich PdO(100)','Pd$^{2+}$ Pd-top'),
        ('S3','S3_PdO100','O-rich PdO(100)','Pd$^{2+}$ O-top'),
        ('S4','S4_PdO2_110','PdO$_2$(110)','Pd$^{4+}$')]


def intramol_valid(atoms, ads):
    syms = atoms.get_chemical_symbols()
    c_idx = [i for i,s in enumerate(syms) if s=='C']
    h_idx = [i for i,s in enumerate(syms) if s=='H']
    o_idx = [i for i,s in enumerate(syms) if s=='O']
    if not c_idx: return False
    c = c_idx[-1]
    if ads == 'CO':
        d = min(atoms.get_distances(c, o_idx, mic=True))
        return 1.05 <= d <= 1.30
    else:
        if len(h_idx) != 3: return False
        d_oc = min(atoms.get_distances(c, o_idx, mic=True))
        d_ch = [atoms.get_distance(c,h,mic=True) for h in h_idx]
        return 1.30 <= d_oc <= 1.55 and all(0.90<=d<=1.25 for d in d_ch)


def slab_intact(atoms, slab, n_ads, max_disp=0.8):
    """Substrate Pd atoms must not displace > max_disp Å from clean slab."""
    n_sub = len(slab)
    if len(atoms) != n_sub + n_ads: return False
    sub = atoms[:n_sub]
    disp = np.linalg.norm(sub.positions - slab.positions, axis=1)
    return float(disp.max()) < max_disp


def top1_valid_atoms(sid, sdir, ads):
    n_ads = 2 if ads=='CO' else 5
    unique = json.load(open(G3 / sdir / f'MLIP_phase1/unique_{ads}.json'))
    traj = list(read(G3 / sdir / f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
    slab = read(G2 / sdir / 'CONTCAR')
    for r in sorted(unique, key=lambda r: r['E']):
        atoms = traj[r['idx']]
        if len(atoms) != len(slab)+n_ads:
            ads_atoms = atoms[-n_ads:]
            atoms = slab.copy(); atoms += ads_atoms
        if intramol_valid(atoms, ads) and slab_intact(atoms, slab, n_ads):
            print(f'  {sid} {ads}: idx={r["idx"]} E={r["E"]:.3f}')
            return atoms * (2,2,1)
    print(f'  {sid} {ads}: NO valid+intact candidate'); return None


def _unused_render_custom(atoms, out_png, rotation='-90x,-90y,0z', width=2400):
    syms = atoms.get_chemical_symbols()
    colors = [COLORS.get(s, [0.5,0.5,0.5]) for s in syms]
    radii = [RADII.get(s, 0.5) for s in syms]
    pov = out_png.with_suffix('.pov')
    ini = out_png.with_suffix('.ini')
    write_pov(str(pov), atoms, rotation=rotation, radii=radii, show_unit_cell=0,
              colors=colors,
              povray_settings={'textures':['ase3']*len(atoms),
                                'canvas_width': width, 'transparent': False})
    # tiny-angle perspective (orthographic-equivalent)
    pos = atoms.positions
    cell = atoms.cell.array
    # apply rotation to get rendered xy extents
    from ase.io.utils import rotate
    R = rotate(rotation)
    pp = pos @ R
    margin = 2.0
    xr = pp[:,0].max() - pp[:,0].min() + 2*margin
    yr = pp[:,1].max() - pp[:,1].min() + 2*margin
    cx = (pp[:,0].max()+pp[:,0].min())/2
    cy = (pp[:,1].max()+pp[:,1].min())/2
    half = max(xr, yr)/2
    angle = 1.0
    cam_dist = half / np.tan(np.radians(angle/2))
    txt = pov.read_text()
    cam_block = (f'camera {{perspective\n  angle {angle:.1f}\n'
                 f'  right -{xr/yr:.4f}*x up 1.0*y\n'
                 f'  location <{cx:.2f},{cy:.2f},{cam_dist:.1f}>'
                 f' look_at <{cx:.2f},{cy:.2f},0>}}')
    txt = re.sub(r'camera \{[^}]*\}', cam_block, txt)
    txt = re.sub(r'(light_source \{<)\s*[-\d.]+,\s*[-\d.]+,\s*[-\d.]+(>)',
                 lambda m: f'{m.group(1)}{cx+2:.2f},{cy+3:.2f},{cam_dist-10:.1f}{m.group(2)}', txt)
    pov.write_text(txt)
    # adjust ini width
    txt_ini = ini.read_text()
    h_px = int(width * yr / xr)
    txt_ini = re.sub(r'Width=\d+', f'Width={width}', txt_ini)
    txt_ini = re.sub(r'Height=\d+', f'Height={h_px}', txt_ini)
    ini.write_text(txt_ini)
    ret = subprocess.run([POVRAY, ini.name], cwd=str(out_png.parent),
                         capture_output=True, text=True)
    for f in [pov, ini]:
        if f.exists(): f.unlink()
    return out_png.exists()


def make_row(ads, label):
    pngs = []
    for sid, sdir, name, ox in SURF:
        atoms = top1_valid_atoms(sid, sdir, ads)
        if atoms is None:
            print(f'{sid} {ads}: no valid'); continue
        png = OUT / f'{sid}_{ads}_top1.png'
        from render_structure import render
        render(atoms, png, rotation='-90x,-90y,0z', width=2400, show_cell=False)
        pngs.append((sid, name, ox, png))
    fig, axes = plt.subplots(1, len(pngs), figsize=(4.4*len(pngs), 5.5))
    if len(pngs) == 1: axes = [axes]
    for ax, (sid, name, ox, png) in zip(axes, pngs):
        ax.imshow(mpimg.imread(png)); ax.set_axis_off()
        ax.set_title(f'{sid}: {name}\n{ox}', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUT/f'ALL_{ads}_top1.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{ads} row ✓')

make_row('CO', 'CO*')
make_row('CH3O', 'CH3O*')
