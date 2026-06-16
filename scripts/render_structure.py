"""POSCAR → POV-Ray rendered PNG (VESTA-quality, ase3 texture).

Accepts one or more POSCAR/CONTCAR paths and produces publication-quality
structure images using ASE + POV-Ray with the ase3 metallic texture.

Usage:
    # Single file
    conda run -n pddmc python scripts/render_structure.py POSCAR

    # Multiple files
    conda run -n pddmc python scripts/render_structure.py S1/POSCAR S2/POSCAR S3/POSCAR

    # Custom rotation & resolution
    conda run -n pddmc python scripts/render_structure.py POSCAR -r "-90x,0y,0z" -w 1200

    # Multi-view (perspective + top + side)
    conda run -n pddmc python scripts/render_structure.py POSCAR --multiview

    # Combined grid of multiple structures
    conda run -n pddmc python scripts/render_structure.py S1/POSCAR S2/POSCAR --grid
"""

import re
import sys
import argparse
import subprocess
import numpy as np
from pathlib import Path
from ase.io import read
from ase.io.pov import write_pov

TEXTURE = "ase3"
DEFAULT_ROTATION = "-75x,15y,3z"
DEFAULT_WIDTH = 800
DEFAULT_RADII = 0.85

ATOM_COLORS = {
    "Pd": [0.165, 0.616, 0.561],  # #2a9d8f
    "O":  [0.941, 0.443, 0.404],  # #f07167
}

VIEWS = {
    "perspective": "-75x,15y,3z",
    "top":         "0x,0y,0z",
    "side":        "-90x,0y,0z",
}


def _parse_pov_coords(pov_path):
    """Extract atom and cylinder coords from POV file (render coordinates)."""
    txt = pov_path.read_text()
    atoms = []
    for m in re.finditer(r'atom\(<\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)>', txt):
        atoms.append([float(m.group(1)), float(m.group(2)), float(m.group(3))])
    cyls = []
    for m in re.finditer(r'cylinder\s*\{<\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)>', txt):
        cyls.append([float(m.group(1)), float(m.group(2)), float(m.group(3))])
    return (np.array(atoms) if atoms else None,
            np.array(cyls) if cyls else None)


def _fix_pov_camera(pov_path, atoms, rotation):
    """Rewrite camera: near-orthographic perspective framed to atom positions.

    POV-Ray 3.7.0.10 has a bug: orthographic cameras always fail the
    "angle >= 180" check.  Workaround: use perspective with angle=1 and
    a far camera (< 2% distortion, lattice lines appear parallel).
    For tilted views use angle=20 for depth cue.
    """
    atom_coords, cyl_coords = _parse_pov_coords(pov_path)
    if atom_coords is None or len(atom_coords) == 0:
        return

    is_side_or_top = rotation in ("-90x,0y,0z", "0x,0y,0z",
                                  "-90x,0y,0z", "90x,0y,0z")
    angle = 1.0 if is_side_or_top else 20.0

    all_xy = atom_coords[:, :2]

    margin = 2.0
    x_range = all_xy[:, 0].max() - all_xy[:, 0].min() + 2 * margin
    y_range = all_xy[:, 1].max() - all_xy[:, 1].min() + 2 * margin
    cx = (all_xy[:, 0].max() + all_xy[:, 0].min()) / 2
    cy = (all_xy[:, 1].max() + all_xy[:, 1].min()) / 2

    half_span = max(x_range, y_range) / 2
    cam_dist = half_span / np.tan(np.radians(angle / 2))

    txt = pov_path.read_text()
    cam_block = (
        f"camera {{perspective\n"
        f"  angle {angle:.1f}\n"
        f"  right -{x_range / y_range:.4f}*x up 1.0*y\n"
        f"  location <{cx:.2f},{cy:.2f},{cam_dist:.1f}>"
        f" look_at <{cx:.2f},{cy:.2f},0>}}"
    )
    txt = re.sub(r'camera \{[^}]*\}', cam_block, txt)
    txt = re.sub(
        r'(light_source \{<)\s*[-\d.]+,\s*[-\d.]+,\s*[-\d.]+(>)',
        lambda m: f'{m.group(1)}{cx + 2:.2f},{cy + 3:.2f},{cam_dist - 10:.1f}{m.group(2)}',
        txt,
    )
    pov_path.write_text(txt)


def _fix_ini_aspect(ini_path, pov_path, rotation=""):
    """Set ini Width/Height to match the camera viewport aspect ratio."""
    atom_coords, cyl_coords = _parse_pov_coords(pov_path)
    if atom_coords is None or len(atom_coords) == 0:
        return

    all_xy = atom_coords[:, :2]

    margin = 2.0
    x_range = all_xy[:, 0].max() - all_xy[:, 0].min() + 2 * margin
    y_range = all_xy[:, 1].max() - all_xy[:, 1].min() + 2 * margin
    ratio = y_range / max(x_range, 0.01)

    txt = ini_path.read_text()
    lines = txt.split("\n")
    width = 800
    for line in lines:
        if line.startswith("Width="):
            width = int(float(line.split("=")[1]))
    height = max(100, int(width * ratio))

    new_lines = []
    for line in lines:
        if line.startswith("Width="):
            new_lines.append(f"Width={width}")
        elif line.startswith("Height="):
            new_lines.append(f"Height={height}")
        else:
            new_lines.append(line)
    ini_path.write_text("\n".join(new_lines))


def render(atoms, out_png, rotation=DEFAULT_ROTATION, width=DEFAULT_WIDTH,
           show_cell=True):
    """Render atoms to PNG via POV-Ray with ase3 texture.

    Returns Path to PNG on success, None on failure.
    """
    out_png = Path(out_png)
    work_dir = out_png.parent
    stem = out_png.stem

    pov_path = work_dir / f"{stem}.pov"
    ini_path = work_dir / f"{stem}.ini"

    syms = atoms.get_chemical_symbols()
    colors = [ATOM_COLORS.get(s, [0.5, 0.5, 0.5]) for s in syms]

    write_pov(
        str(pov_path), atoms,
        rotation=rotation,
        radii=DEFAULT_RADII,
        show_unit_cell=2 if show_cell else 0,
        colors=colors,
        povray_settings={
            "textures": [TEXTURE] * len(atoms),
            "canvas_width": width,
            "transparent": False,
        },
    )
    _fix_pov_camera(pov_path, atoms, rotation)
    _fix_ini_aspect(ini_path, pov_path, rotation)

    ret = subprocess.run(
        ["povray", ini_path.name],
        cwd=str(work_dir),
        capture_output=True, text=True,
    )

    # Cleanup temp files
    for ext in [".pov", ".ini"]:
        f = work_dir / f"{stem}{ext}"
        if f.exists():
            f.unlink()

    if ret.returncode != 0:
        err = [l for l in ret.stderr.strip().split("\n") if "Error" in l]
        print(f"  FAILED: {err[-1] if err else 'unknown error'}")
        return None

    if out_png.exists():
        return out_png
    return None


def render_multiview(atoms, out_dir, name, width=DEFAULT_WIDTH):
    """Render perspective / top / side views and combine into one PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    out_dir = Path(out_dir)
    panels = {}

    for view_name, rot in VIEWS.items():
        tmp_png = out_dir / f"_tmp_{name}_{view_name}.png"
        result = render(atoms, tmp_png, rotation=rot, width=width)
        if result:
            panels[view_name] = result

    if not panels:
        print(f"  No views rendered for {name}")
        return None

    fig, axes = plt.subplots(1, len(panels), figsize=(6 * len(panels), 6))
    if len(panels) == 1:
        axes = [axes]

    for ax, (view_name, png_path) in zip(axes, panels.items()):
        img = mpimg.imread(str(png_path))
        ax.imshow(img)
        ax.set_title(view_name.capitalize(), fontsize=13, fontweight="bold")
        ax.axis("off")
        png_path.unlink()

    plt.tight_layout()
    out = out_dir / f"{name}_views.png"
    plt.savefig(str(out), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    return out


def render_grid(file_list, out_path, width=DEFAULT_WIDTH):
    """Render multiple structures into a single grid image."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    out_path = Path(out_path)
    n = len(file_list)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]

    for idx, poscar_path in enumerate(file_list):
        row, col = idx // ncols, idx % ncols
        ax = axes[row, col]

        atoms = read(str(poscar_path))
        name = Path(poscar_path).parent.name or Path(poscar_path).stem
        tmp_png = out_path.parent / f"_tmp_grid_{idx}.png"

        result = render(atoms, tmp_png, width=width)
        if result:
            img = mpimg.imread(str(result))
            ax.imshow(img)
            result.unlink()

        syms = atoms.get_chemical_symbols()
        comp = {}
        for s in syms:
            comp[s] = comp.get(s, 0) + 1
        comp_str = ", ".join(f"{k}={v}" for k, v in comp.items())
        ax.set_title(f"{name}\n({len(atoms)} atoms: {comp_str})",
                     fontsize=10, fontweight="bold")
        ax.axis("off")

    for idx in range(n, nrows * ncols):
        axes[idx // ncols, idx % ncols].axis("off")

    plt.tight_layout()
    plt.savefig(str(out_path), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="POSCAR → POV-Ray PNG (ase3 texture, VESTA-quality)")
    parser.add_argument("files", nargs="+", help="POSCAR/CONTCAR file paths")
    parser.add_argument("-r", "--rotation", default=DEFAULT_ROTATION,
                        help=f"Rotation string (default: {DEFAULT_ROTATION})")
    parser.add_argument("-w", "--width", type=int, default=DEFAULT_WIDTH,
                        help=f"Image width in px (default: {DEFAULT_WIDTH})")
    parser.add_argument("-o", "--outdir", default=None,
                        help="Output directory (default: same as input)")
    parser.add_argument("--multiview", action="store_true",
                        help="Render 3 views (perspective/top/side)")
    parser.add_argument("--grid", action="store_true",
                        help="Combine all inputs into one grid image")
    args = parser.parse_args()

    if subprocess.run(["which", "povray"], capture_output=True).returncode != 0:
        print("ERROR: povray not found. Install: conda install -c conda-forge povray")
        sys.exit(1)

    files = [Path(f) for f in args.files]
    for f in files:
        if not f.exists():
            print(f"ERROR: {f} not found")
            sys.exit(1)

    if args.grid and len(files) > 1:
        out_dir = Path(args.outdir) if args.outdir else files[0].parent
        out_path = out_dir / "structure_grid.png"
        print(f"Rendering grid ({len(files)} structures)...")
        result = render_grid(files, out_path, width=args.width)
        if result:
            print(f"Saved: {result}")
        return

    for f in files:
        atoms = read(str(f))
        name = f.parent.name if f.parent.name != "." else f.stem
        out_dir = Path(args.outdir) if args.outdir else f.parent

        if args.multiview:
            print(f"Rendering {name} (3 views)...")
            result = render_multiview(atoms, out_dir, name, width=args.width)
        else:
            out_png = out_dir / f"{name}.png"
            print(f"Rendering {name}...", end=" ", flush=True)
            result = render(atoms, out_png, rotation=args.rotation,
                           width=args.width)

        if result:
            print(f"Saved: {result} ({result.stat().st_size/1024:.0f} KB)")
        else:
            print(f"Failed: {name}")


if __name__ == "__main__":
    main()
