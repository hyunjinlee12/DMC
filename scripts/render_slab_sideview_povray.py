"""POV-Ray side views — orthographic, no perspective."""
import subprocess
from pathlib import Path
from ase.io import read
from ase.io.pov import write_pov

PROJECT = Path(__file__).resolve().parent.parent
G2 = PROJECT / "calculations" / "G2_slab"
OUT = PROJECT / "reports/predft_advisor_figures/slab_sideviews_povray"
OUT.mkdir(parents=True, exist_ok=True)

SLABS = [
    ("S1_Pd100",        "S1"),
    ("S2_PdO101_Pd100", "S2"),
    ("S3_PdO100",       "S3"),
    ("S3b_PdO100_PdOterm","S3b"),
    ("S4_PdO2_110",     "S4"),
]

for dirname, tag in SLABS:
    atoms = read(G2 / dirname / "CONTCAR")
    pov_path = OUT / f"{tag}_side.pov"
    ini_path = OUT / f"{tag}_side.ini"
    png_path = OUT / f"{tag}_side.png"
    write_pov(
        str(pov_path), atoms,
        rotation="-90x,0y,0z",
        radii=0.85, show_unit_cell=0,
        povray_settings={
            "textures": ["ase3"]*len(atoms),
            "canvas_width": 1200,
            "transparent": False,
        },
    )
    # Fake-orthographic via very small perspective angle (avoids ASE's broken orthographic block)
    txt = pov_path.read_text()
    import re
    txt = re.sub(
        r'camera \{orthographic\s*\n\s*right ([^\n]+)\n\s*direction[^\n]*\n\s*location ([^\}]+)\}',
        lambda m: f'camera {{perspective\n  angle 2\n  right {m.group(1)}\n  location {m.group(2)}}}',
        txt,
    )
    # Move camera back far for true-orthographic-look (angle 2, large distance)
    txt = re.sub(r'location <0,0,50.00>', 'location <0,0,2000>', txt)
    pov_path.write_text(txt)
    ret = subprocess.run(["/home/hyunjin/.conda/envs/pddmc/bin/povray", ini_path.name],
                         cwd=str(OUT), capture_output=True, text=True)
    if ret.returncode == 0 and png_path.exists():
        print(f"{tag} OK ({png_path.stat().st_size/1024:.0f} KB)")
    else:
        print(f"{tag} FAIL: {ret.stderr.strip().splitlines()[-1] if ret.stderr else '?'}")
