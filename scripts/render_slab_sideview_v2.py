"""Side view (orthographic, no perspective) — same ase3 style as render_all_povray.py."""
import re, subprocess
from pathlib import Path
from ase.io import read
from ase.io.pov import write_pov

OUT = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/predft_advisor_figures/slab_sideviews_v2")
OUT.mkdir(parents=True, exist_ok=True)
G2 = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/G2_slab")
POVRAY = "/home/hyunjin/.conda/envs/pddmc/bin/povray"

SLABS = [("S1_Pd100","S1"),("S2_PdO101_Pd100","S2"),("S3_PdO100","S3"),
         ("S3b_PdO100_PdOterm","S3b"),("S4_PdO2_110","S4")]

for dirname, tag in SLABS:
    atoms = read(G2 / dirname / "CONTCAR")
    pov = OUT / f"{tag}_side.pov"
    ini = OUT / f"{tag}_side.ini"
    png = OUT / f"{tag}_side.png"
    write_pov(str(pov), atoms,
              rotation="-90x,0y,0z",
              radii=0.85, show_unit_cell=0,
              povray_settings={"textures":["ase3"]*len(atoms),
                                "canvas_width":1200,"transparent":False})
    # Convert ASE's broken orthographic block → perspective with tiny angle from far away
    # (visually identical to orthographic, but POV-Ray accepts it).
    txt = pov.read_text()
    txt = re.sub(
        r'camera \{orthographic\s*\n\s*right (-?[\d.]+)\*x[^\n]*\n(\s*direction[^\n]*\n)?\s*location <([^>]+)>[^\}]*\}',
        lambda m: (f'camera {{perspective\n'
                   f'  angle 1\n'
                   f'  right {abs(float(m.group(1)))}*x\n'
                   f'  location <0,0,5000>\n  look_at <0,0,0>}}'),
        txt)
    pov.write_text(txt)
    ret = subprocess.run([POVRAY, ini.name], cwd=str(OUT),
                         capture_output=True, text=True)
    print(f'{tag} {"OK" if ret.returncode==0 else "FAIL"}')
