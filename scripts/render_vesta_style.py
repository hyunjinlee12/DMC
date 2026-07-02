"""VESTA-like render with bonds drawn between close atom pairs."""
import sys, re, subprocess, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
from pathlib import Path
from ase.io import read
from ase.io.pov import write_pov

OUT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/predft_advisor_figures/extra_surfaces_v2')
OUT.mkdir(parents=True, exist_ok=True)
POVRAY = "/home/hyunjin/.conda/envs/pddmc/bin/povray"

ATOM_COLORS = {
    "Pd": [0.165, 0.616, 0.561],   # teal
    "O":  [0.941, 0.443, 0.404],   # coral
    "C":  [0.18, 0.18, 0.18],
    "H":  [0.95, 0.95, 0.95],
}
ATOM_RADII = {"Pd": 0.95, "O": 0.55, "C": 0.5, "H": 0.4}
BOND_CUTOFF = {"Pd-O": 2.5, "Pd-Pd": 3.2, "C-O": 1.5, "C-H": 1.3, "O-H": 1.1}


def bond_list(atoms):
    """Find atom pairs within bonding distance."""
    syms = atoms.get_chemical_symbols()
    n = len(atoms)
    bonds = []
    for i in range(n):
        for j in range(i+1, n):
            pair = '-'.join(sorted([syms[i], syms[j]]))
            cutoff = BOND_CUTOFF.get(pair, 0)
            if cutoff == 0: continue
            d = atoms.get_distance(i, j, mic=True)
            if d < cutoff:
                bonds.append((i, j, (0,0,0), 1))   # (i, j, mic_offset, bond_order)
    return bonds


def render_vesta(atoms, out_png, rotation='-90x,-90y,0z', width=2400):
    syms = atoms.get_chemical_symbols()
    colors = [ATOM_COLORS.get(s, [0.5]*3) for s in syms]
    radii = [ATOM_RADII.get(s, 0.5) for s in syms]
    bonds = bond_list(atoms)
    pov = out_png.with_suffix('.pov')
    ini = out_png.with_suffix('.ini')
    write_pov(str(pov), atoms,
              rotation=rotation, radii=radii, show_unit_cell=0,
              colors=colors,
              povray_settings={'textures':['ase3']*len(atoms),
                                'canvas_width': width, 'transparent': False,
                                'bondatoms': bonds, 'bondlinewidth': 0.10})
    # fix camera (small-angle perspective)
    txt = pov.read_text()
    pos = atoms.positions
    from ase.io.utils import PlottingVariables
    # Use simple bounds from atoms
    R = np.array([[0,0,1],[1,0,0],[0,1,0]])   # just to estimate
    # Use ase.io.utils to get rotated coords -- simplified manual rotation
    # rotation -90x,-90y,0z: apply Rx(-90) then Ry(-90)
    import re as _re
    coords_match = _re.findall(r'atom\(<\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)>', txt)
    atom_coords = np.array([[float(x),float(y),float(z)] for x,y,z in coords_match])
    if len(atom_coords) == 0: return False
    xs, ys = atom_coords[:,0], atom_coords[:,1]
    margin = 2.0
    xr = xs.max() - xs.min() + 2*margin
    yr = ys.max() - ys.min() + 2*margin
    cx = (xs.max()+xs.min())/2
    cy = (ys.max()+ys.min())/2
    half = max(xr, yr) / 2
    angle = 1.0
    cam_dist = half / np.tan(np.radians(angle/2))
    cam = (f'camera {{perspective\n  angle {angle}\n'
           f'  right -{xr/yr:.4f}*x up 1.0*y\n'
           f'  location <{cx:.2f},{cy:.2f},{cam_dist:.1f}>'
           f' look_at <{cx:.2f},{cy:.2f},0>}}')
    txt = _re.sub(r'camera \{[^}]*\}', cam, txt)
    txt = _re.sub(r'(light_source \{<)\s*[-\d.]+,\s*[-\d.]+,\s*[-\d.]+(>)',
                  lambda m: f'{m.group(1)}{cx+2:.2f},{cy+3:.2f},{cam_dist-10:.1f}{m.group(2)}', txt)
    pov.write_text(txt)
    # aspect ini
    ini_txt = ini.read_text()
    h_px = int(width * yr / xr)
    ini_txt = re.sub(r'Width=\d+', f'Width={width}', ini_txt)
    ini_txt = re.sub(r'Height=\d+', f'Height={h_px}', ini_txt)
    ini.write_text(ini_txt)
    ret = subprocess.run([POVRAY, ini.name], cwd=str(out_png.parent),
                         capture_output=True, text=True)
    for f in [pov, ini]:
        if f.exists(): f.unlink()
    return out_png.exists()


for tag, p in [('S5_vesta', '/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab/S5_Pd111/POSCAR'),
                ('S6_vesta', '/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab/S6_PdO2_110_Ocus/POSCAR')]:
    a = read(p) * (2,2,1)
    ok = render_vesta(a, OUT/f'{tag}.png', width=2400)
    print(f'{tag}: {"OK" if ok else "FAIL"}')
