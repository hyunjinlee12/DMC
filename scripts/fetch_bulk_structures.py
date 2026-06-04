"""Fetch Pd, PdO, PdO2 bulk structures from Materials Project.

Usage:
    conda run -n pddmc python scripts/fetch_bulk_structures.py

Saves conventional standard structures as POSCAR and .json in structures/.
"""

import os
import json
from pathlib import Path
from mp_api.client import MPRester
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
import ase.io

STRUCTURES_DIR = Path(__file__).resolve().parent.parent / "structures"
STRUCTURES_DIR.mkdir(exist_ok=True)

TARGETS = {
    "Pd": {
        "formula": "Pd",
        "description": "fcc Pd bulk",
        "spacegroup": "Fm-3m",
    },
    "PdO": {
        "formula": "PdO",
        "description": "tetragonal PdO bulk (P4_2/mmc)",
        "spacegroup": "P4_2/mmc",
    },
    "PdO2": {
        "formula": "PdO2",
        "description": "PdO2 bulk (rutile or hydrophilite-like)",
        "spacegroup": None,
    },
}

EXPERIMENTAL_LATTICE = {
    "Pd": {"a": 3.89},
    "PdO": {"a": 3.04, "c": 5.33},
    "PdO2": {"note": "polymorph-dependent"},
}


def fetch_structures():
    api_key = os.environ.get("MP_API_KEY")
    if not api_key:
        raise EnvironmentError("MP_API_KEY not set")

    results = {}
    with MPRester(api_key) as mpr:
        for name, info in TARGETS.items():
            print(f"\n--- Fetching {name}: {info['description']} ---")
            docs = mpr.materials.summary.search(
                formula=info["formula"],
                fields=["material_id", "structure", "symmetry", "energy_above_hull",
                         "formation_energy_per_atom", "band_gap"],
            )

            if info["spacegroup"]:
                filtered = [d for d in docs
                            if d.symmetry.symbol == info["spacegroup"]]
            else:
                filtered = sorted(docs, key=lambda d: d.energy_above_hull)

            if not filtered:
                print(f"  WARNING: No structures found for {name}, using lowest E_hull")
                filtered = sorted(docs, key=lambda d: d.energy_above_hull)

            for i, doc in enumerate(filtered[:3]):
                print(f"  [{i}] {doc.material_id} | {doc.symmetry.symbol} | "
                      f"E_hull={doc.energy_above_hull:.4f} eV/atom | "
                      f"gap={doc.band_gap:.2f} eV")

            best = filtered[0]
            struct = best.structure.get_conventional_standard_structure()

            outname = f"{name}_bulk_{best.material_id}"
            poscar_path = STRUCTURES_DIR / f"{outname}.vasp"
            json_path = STRUCTURES_DIR / f"{outname}.json"

            struct.to(str(poscar_path), fmt="poscar")
            struct.to(str(json_path))

            atoms = AseAtomsAdaptor.get_atoms(struct)
            ase.io.write(STRUCTURES_DIR / f"{outname}_ase.traj", atoms)

            a, b, c = struct.lattice.abc
            alpha, beta, gamma = struct.lattice.angles
            print(f"  Selected: {best.material_id} ({best.symmetry.symbol})")
            print(f"  Lattice: a={a:.4f} b={b:.4f} c={c:.4f}")
            print(f"  Angles:  {alpha:.1f} {beta:.1f} {gamma:.1f}")
            print(f"  Atoms:   {struct.num_sites}")
            print(f"  Saved:   {poscar_path}")

            results[name] = {
                "mp_id": str(best.material_id),
                "spacegroup": best.symmetry.symbol,
                "e_hull": best.energy_above_hull,
                "lattice_abc": [a, b, c],
                "n_atoms": struct.num_sites,
                "poscar": str(poscar_path),
            }

    summary_path = STRUCTURES_DIR / "bulk_fetch_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSummary saved: {summary_path}")
    return results


if __name__ == "__main__":
    fetch_structures()
