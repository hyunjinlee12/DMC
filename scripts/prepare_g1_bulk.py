"""Prepare VASP input files for G1 bulk optimization (T1.1-T1.3).

Reads fetched structures from structures/, generates POSCAR/INCAR/KPOINTS/POTCAR
in calculations/G1_bulk/{Pd,PdO,PdO2}/.

Usage:
    conda run -n pddmc python scripts/prepare_g1_bulk.py

Does NOT submit jobs — review inputs first.
"""

import json
from pathlib import Path

import ase.io
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from incar_templates import write_incar
from vasp_utils import write_potcar, write_kpoints_gamma, get_unique_elements

PROJECT = Path(__file__).resolve().parent.parent
STRUCTURES = PROJECT / "structures"
CALC_DIR = PROJECT / "calculations" / "G1_bulk"

BULK_CONFIG = {
    "Pd": {
        "template": "bulk_metal",
        "kpoints": (12, 12, 12),
        "comment": "fcc Pd, ISMEAR=1 (metal)",
    },
    "PdO": {
        "template": "bulk_oxide",
        "kpoints": (8, 8, 6),
        "comment": "tetragonal PdO, ISMEAR=0 (oxide)",
    },
    "PdO2": {
        "template": "bulk_oxide",
        "kpoints": (6, 6, 8),
        "comment": "PdO2, ISMEAR=0 (oxide), k-mesh adjusted to cell shape",
    },
}


def prepare():
    summary_path = STRUCTURES / "bulk_fetch_summary.json"
    if not summary_path.exists():
        print("ERROR: Run fetch_bulk_structures.py first")
        return

    with open(summary_path) as f:
        summary = json.load(f)

    for name, config in BULK_CONFIG.items():
        if name not in summary:
            print(f"WARNING: {name} not in fetch summary, skipping")
            continue

        info = summary[name]
        outdir = CALC_DIR / name
        outdir.mkdir(parents=True, exist_ok=True)

        struct = Structure.from_file(info["poscar"])
        atoms = AseAtomsAdaptor.get_atoms(struct)
        elements = get_unique_elements(atoms)

        ase.io.write(str(outdir / "POSCAR"), atoms, format="vasp")
        write_incar(outdir / "INCAR", config["template"])
        write_kpoints_gamma(config["kpoints"], outdir)
        write_potcar(elements, outdir)

        print(f"[{name}] {config['comment']}")
        print(f"  Structure: {info['mp_id']} ({info['spacegroup']}), {info['n_atoms']} atoms")
        print(f"  KPOINTS: {config['kpoints']} (Gamma)")
        print(f"  Elements: {elements}")
        print(f"  Output: {outdir}")
        print()

    print("Review inputs, then submit with:")
    print("  sbatch -J Pd_bulk --chdir=calculations/G1_bulk/Pd scripts/submit_vasp.sh")
    print("  sbatch -J PdO_bulk --chdir=calculations/G1_bulk/PdO scripts/submit_vasp.sh")
    print("  sbatch -J PdO2_bulk --chdir=calculations/G1_bulk/PdO2 scripts/submit_vasp.sh")


if __name__ == "__main__":
    prepare()
