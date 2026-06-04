"""Generate ENCUT and k-mesh convergence test inputs for G1 bulk (T1.4).

Single-point (NSW=0) calculations at the MP geometry.
Step 1: ENCUT sweep at fixed k-mesh
Step 2: k-mesh sweep at fixed ENCUT (520 eV)

Usage:
    conda run -n pddmc python scripts/prepare_convergence.py
"""

import json
import shutil
from pathlib import Path

import ase.io
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from incar_templates import get_incar
from vasp_utils import write_potcar, write_kpoints_gamma, get_unique_elements

PROJECT = Path(__file__).resolve().parent.parent
STRUCTURES = PROJECT / "structures"
CONV_DIR = PROJECT / "calculations" / "G1_bulk" / "convergence"

ENCUT_VALUES = [400, 450, 500, 520, 550, 600]

BULK_CONFIG = {
    "Pd": {
        "template": "bulk_metal",
        "base_kpts": (12, 12, 12),
        "kmesh_sweep": [
            (6, 6, 6), (8, 8, 8), (10, 10, 10),
            (12, 12, 12), (14, 14, 14), (16, 16, 16),
        ],
    },
    "PdO": {
        "template": "bulk_oxide",
        "base_kpts": (8, 8, 6),
        "kmesh_sweep": [
            (4, 4, 3), (6, 6, 4), (8, 8, 6),
            (10, 10, 8), (12, 12, 10),
        ],
    },
    "PdO2": {
        "template": "bulk_oxide",
        "base_kpts": (6, 6, 8),
        "kmesh_sweep": [
            (4, 4, 4), (4, 4, 6), (6, 6, 8),
            (8, 8, 10), (10, 10, 12),
        ],
    },
}


def prepare():
    summary_path = STRUCTURES / "bulk_fetch_summary.json"
    with open(summary_path) as f:
        summary = json.load(f)

    all_dirs = []

    for name, config in BULK_CONFIG.items():
        info = summary[name]
        struct = Structure.from_file(info["poscar"])
        atoms = AseAtomsAdaptor.get_atoms(struct)
        elements = get_unique_elements(atoms)

        # --- ENCUT sweep ---
        for encut in ENCUT_VALUES:
            tag = f"{name}/encut_{encut}"
            outdir = CONV_DIR / tag
            outdir.mkdir(parents=True, exist_ok=True)

            ase.io.write(str(outdir / "POSCAR"), atoms, format="vasp")
            incar = get_incar(config["template"], overrides={
                "NSW": 0,
                "IBRION": -1,
                "ISIF": 0,
                "ENCUT": encut,
                "LWAVE": ".FALSE.",
                "LCHARG": ".FALSE.",
            })
            (outdir / "INCAR").write_text(incar)
            write_kpoints_gamma(config["base_kpts"], outdir)
            write_potcar(elements, outdir)
            all_dirs.append(str(outdir))

        # --- k-mesh sweep ---
        for kpts in config["kmesh_sweep"]:
            ktag = f"{kpts[0]}x{kpts[1]}x{kpts[2]}"
            tag = f"{name}/kpts_{ktag}"
            outdir = CONV_DIR / tag
            outdir.mkdir(parents=True, exist_ok=True)

            ase.io.write(str(outdir / "POSCAR"), atoms, format="vasp")
            incar = get_incar(config["template"], overrides={
                "NSW": 0,
                "IBRION": -1,
                "ISIF": 0,
                "LWAVE": ".FALSE.",
                "LCHARG": ".FALSE.",
            })
            (outdir / "INCAR").write_text(incar)
            write_kpoints_gamma(kpts, outdir)
            write_potcar(elements, outdir)
            all_dirs.append(str(outdir))

        print(f"[{name}] ENCUT: {len(ENCUT_VALUES)} pts, k-mesh: {len(config['kmesh_sweep'])} pts")

    # Write batch submission script
    batch_script = CONV_DIR / "submit_all.sh"
    lines = ["#!/bin/bash", f"# G1 convergence tests: {len(all_dirs)} jobs total", ""]
    for d in all_dirs:
        jobname = d.split("convergence/")[-1].replace("/", "_")
        lines.append(f'sbatch -J "cv_{jobname}" --ntasks-per-node=4 --chdir="{d}" '
                      f'{PROJECT}/scripts/submit_vasp.sh')
    batch_script.write_text("\n".join(lines) + "\n")

    print(f"\nTotal: {len(all_dirs)} single-point calculations")
    print(f"Batch script: {batch_script}")
    print(f"Run: bash {batch_script}")


if __name__ == "__main__":
    prepare()
