"""Analyze ENCUT and k-mesh convergence results for G1 (T1.4).

Reads OSZICAR from each convergence directory, extracts final energy,
produces convergence table + plots.

Usage:
    conda run -n pddmc python scripts/analyze_convergence.py
"""

import re
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONV_DIR = Path(__file__).resolve().parent.parent / "calculations" / "G1_bulk" / "convergence"
OUT_DIR = CONV_DIR / "results"


def read_oszicar_energy(calcdir):
    oszicar = Path(calcdir) / "OSZICAR"
    if not oszicar.exists():
        return None
    lines = oszicar.read_text().strip().split("\n")
    for line in reversed(lines):
        match = re.search(r"E0=\s*([-\d.]+E?[+-]?\d*)", line)
        if match:
            return float(match.group(1))
    return None


def read_natoms(calcdir):
    poscar = Path(calcdir) / "POSCAR"
    if not poscar.exists():
        return 1
    lines = poscar.read_text().strip().split("\n")
    counts = lines[6].split() if len(lines) > 6 else lines[5].split()
    return sum(int(x) for x in counts)


def collect_results():
    OUT_DIR.mkdir(exist_ok=True)
    materials = ["Pd", "PdO", "PdO2"]
    all_data = {}

    for mat in materials:
        mat_dir = CONV_DIR / mat
        if not mat_dir.exists():
            print(f"WARNING: {mat_dir} not found, skipping")
            continue

        encut_data = []
        kpts_data = []

        for d in sorted(mat_dir.iterdir()):
            if not d.is_dir():
                continue
            energy = read_oszicar_energy(d)
            natoms = read_natoms(d)
            if energy is None:
                print(f"  {d.name}: no energy (not finished?)")
                continue

            e_per_atom = energy / natoms

            if d.name.startswith("encut_"):
                encut = int(d.name.split("_")[1])
                encut_data.append({"ENCUT": encut, "E_total": energy, "E_per_atom": e_per_atom})
            elif d.name.startswith("kpts_"):
                ktag = d.name.split("_")[1]
                parts = ktag.split("x")
                k_total = 1
                for p in parts:
                    k_total *= int(p)
                kpts_data.append({"k-mesh": ktag, "k_total": k_total,
                                  "E_total": energy, "E_per_atom": e_per_atom})

        if encut_data:
            df_e = pd.DataFrame(encut_data).sort_values("ENCUT")
            ref = df_e["E_per_atom"].iloc[-1]
            df_e["dE_meV"] = (df_e["E_per_atom"] - ref) * 1000
            all_data[f"{mat}_encut"] = df_e
            print(f"\n[{mat}] ENCUT convergence (ref = {ref:.6f} eV/atom at {df_e['ENCUT'].iloc[-1]} eV):")
            print(df_e.to_string(index=False))

        if kpts_data:
            df_k = pd.DataFrame(kpts_data).sort_values("k_total")
            ref = df_k["E_per_atom"].iloc[-1]
            df_k["dE_meV"] = (df_k["E_per_atom"] - ref) * 1000
            all_data[f"{mat}_kpts"] = df_k
            print(f"\n[{mat}] k-mesh convergence (ref = {ref:.6f} eV/atom at {df_k['k-mesh'].iloc[-1]}):")
            print(df_k[["k-mesh", "E_total", "E_per_atom", "dE_meV"]].to_string(index=False))

    # --- Plot ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle("G1 Convergence Tests: ENCUT (top) / k-mesh (bottom)", fontsize=14)

    for i, mat in enumerate(materials):
        key_e = f"{mat}_encut"
        if key_e in all_data:
            df = all_data[key_e]
            ax = axes[0, i]
            ax.plot(df["ENCUT"], df["dE_meV"], "o-", color="C0")
            ax.axhline(1.0, color="r", ls="--", alpha=0.5, label="1 meV")
            ax.axhline(-1.0, color="r", ls="--", alpha=0.5)
            ax.set_xlabel("ENCUT (eV)")
            ax.set_ylabel("ΔE (meV/atom)")
            ax.set_title(f"{mat} ENCUT")
            ax.legend()

        key_k = f"{mat}_kpts"
        if key_k in all_data:
            df = all_data[key_k]
            ax = axes[1, i]
            ax.plot(range(len(df)), df["dE_meV"], "s-", color="C1")
            ax.set_xticks(range(len(df)))
            ax.set_xticklabels(df["k-mesh"], rotation=45)
            ax.axhline(1.0, color="r", ls="--", alpha=0.5, label="1 meV")
            ax.axhline(-1.0, color="r", ls="--", alpha=0.5)
            ax.set_xlabel("k-mesh")
            ax.set_ylabel("ΔE (meV/atom)")
            ax.set_title(f"{mat} k-mesh")
            ax.legend()

    plt.tight_layout()
    plot_path = OUT_DIR / "convergence_plots.png"
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlot saved: {plot_path}")
    return all_data


if __name__ == "__main__":
    collect_results()
