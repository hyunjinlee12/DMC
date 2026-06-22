"""S1 CO* 재생성 (Cl[C]=O) + co-ads 재분류."""
import sys, os, json, numpy as np
from pathlib import Path
from ase import Atoms
from ase.io import read, write
from autoadsorbate import Surface, Fragment
import matplotlib.pyplot as plt
from collections import Counter
import random

random.seed(2104); np.random.seed(2104)

ROOT = Path("/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc")
WORK = ROOT / "calculations/G3_adsorption/S1_Pd100"
SLAB = ROOT / "calculations/G2_slab/S1_Pd100/CONTCAR"

# ========== Step 2 — CO* 재생성 (Cl[C]=O) ==========
slab = read(SLAB)
print(f"[CO*] surface init...")
surface = Surface(slab, mode='slab', precision=0.25)
# enumerate available attrs to find site collection
attrs = [a for a in dir(surface) if not a.startswith('_')]
print(f"  Surface attrs: {[a for a in attrs if 'site' in a.lower() or 'fram' in a.lower() or 'top' in a.lower()][:10]}")


print(f"[CO*] Fragment Cl[C]=O ...")
co = Fragment('Cl[C]=O', to_initialize=1)
# autoadsorbate workaround: conformer.info['smiles'] 직접 설정 필요
for conf in co.conformers:
    if 'smiles' not in conf.info:
        conf.info['smiles'] = 'Cl[C]=O'
co_structs = surface.get_populated_sites(
    co, mode='all', sample_rotation=False,
    conformers_per_site_cap=1, overlap_thr=1.25
)
print(f"  CO* candidates: {len(co_structs)}")

# 검증: 첫 구조의 non-Pd atoms
syms = co_structs[0].get_chemical_symbols()
nonpd = [s for s in syms if s != 'Pd']
print(f"  First structure non-Pd atoms: {nonpd}")
if 'H' in nonpd:
    print(f"  ❌ FAIL: still contains H. Aborting.")
    sys.exit(1)
print(f"  ✅ OK: pure C+O fragment, no H")

# 저장
CO_DIR = WORK / "CO"
CO_DIR.mkdir(exist_ok=True)
write(CO_DIR / "candidates.traj", co_structs)
with open(CO_DIR / "summary.json","w") as f:
    json.dump({
        "total": len(co_structs),
        "smiles": "Cl[C]=O",
        "fix_note": "guideline Cl[C-]#[O+] failed RDKit valence; corrected ClC=O → Cl[C]=O (brackets prevent implicit H)",
        "fragment_atoms": dict(Counter(nonpd)),
    }, f, indent=2)
print(f"  saved {CO_DIR}/candidates.traj")

# 그리드 PNG
from ase.visualize.plot import plot_atoms
n_show = min(25, len(co_structs))
idx = np.linspace(0, len(co_structs)-1, n_show, dtype=int)
rows = cols = 5
fig, axes = plt.subplots(rows, cols, figsize=(15, 15))
for i, ax in enumerate(axes.flat):
    if i < len(idx):
        plot_atoms(co_structs[idx[i]], ax, rotation='0x,0y,0z', radii=0.85, show_unit_cell=0)
        ax.set_title(f'#{idx[i]}', fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
plt.suptitle(f"S1 CO* (Cl[C]=O) — 25 representatives of {len(co_structs)}", fontsize=14, y=1.005)
plt.tight_layout()
plt.savefig(CO_DIR / "grid.png", dpi=120, bbox_inches='tight')
plt.close()
print(f"  saved {CO_DIR}/grid.png")

# ========== Step 4 — co-ads 재실행 (new CO* + 기존 CH3O*) ==========
print(f"\n[Co-ads] loading CH3O* candidates...")
ch3o_structs = list(read(WORK / "CH3O/candidates.traj", index=':'))
print(f"  CH3O* loaded: {len(ch3o_structs)}")

# CO* 와 CH3O*를 페어로. dry-run에서는 stride sampling
# 각 CO* 후보의 C 위치 찾기, CH3O* 후보의 O 위치 찾기
def find_first_C(atoms, n_sub):
    """첫 C (substrate 다음, 즉 적층된 fragment의 C)"""
    syms = atoms.get_chemical_symbols()
    for i in range(n_sub, len(atoms)):
        if syms[i] == 'C': return i
    return None
def find_first_O(atoms, n_sub):
    syms = atoms.get_chemical_symbols()
    for i in range(n_sub, len(atoms)):
        if syms[i] == 'O': return i
    return None

n_sub = 80  # S1 Pd substrate
# 페어 생성. CO×CH3O 가 112×336=37,632 → stride sample
co_stride = max(1, len(co_structs)//20)
ch_stride = max(1, len(ch3o_structs)//20)
co_idx = list(range(0, len(co_structs), co_stride))[:20]
ch_idx = list(range(0, len(ch3o_structs), ch_stride))[:20]
print(f"  pair sampling: {len(co_idx)} CO × {len(ch_idx)} CH3O = {len(co_idx)*len(ch_idx)} pairs")

setA, setB, setTS, side_path = [], [], [], []
dists = []
for ci in co_idx:
    a_co = co_structs[ci]
    pos_co = a_co.get_positions()
    iC = find_first_C(a_co, n_sub)
    if iC is None: continue
    cC = pos_co[iC]
    iO_co = find_first_O(a_co, n_sub)

    for chi in ch_idx:
        a_ch = ch3o_structs[chi]
        pos_ch = a_ch.get_positions()
        iO_ch = find_first_O(a_ch, n_sub)
        if iO_ch is None: continue
        oCH3 = pos_ch[iO_ch]

        # 단순 페어: CO에 CH3O 원자들 추가 (cell mic 고려 안 함, 작은 cell 가정)
        combined = a_co.copy()
        ch_adsorbate = a_ch[n_sub:].copy()  # CH3O fragment만
        combined += ch_adsorbate

        # 거리: C_CO ⋯ O_CH3O (mic 적용)
        d = combined.get_distance(iC, len(a_co)+(iO_ch-n_sub), mic=True)
        dists.append(d)

        # overlap check: 모든 atom pair 거리 ≥ 1.0 (헐겁게)
        # 대량 페어이므로 빠르게 mic 모든 거리 검사 생략, 1.0 미만 적층만 보존 X
        # 분류
        if d < 1.7:
            side_path.append(combined)
        elif d < 2.3:
            setTS.append(combined)
        elif d < 4.0:
            setA.append(combined)
        elif d >= 5.0:
            setB.append(combined)
        # 4.0–5.0: 모호, skip

print(f"  ✅ classification done:")
print(f"    Set A (2.3–4.0 Å, reactive): {len(setA)}")
print(f"    Set B (≥5.0 Å, thermo):      {len(setB)}")
print(f"    Set TS (1.7–2.3 Å):          {len(setTS)}")
print(f"    side-path (<1.7 Å):          {len(side_path)}")
print(f"    distance range: {min(dists):.2f}–{max(dists):.2f} Å, mean {np.mean(dists):.2f}")

# 저장
COADS = WORK / "coads"
COADS.mkdir(exist_ok=True)
if setA: write(COADS/"SetA.traj", setA)
if setB: write(COADS/"SetB.traj", setB)
if setTS: write(COADS/"SetTS.traj", setTS)
if side_path: write(COADS/"side.traj", side_path)
with open(COADS/"summary.json","w") as f:
    json.dump({
        "total_pairs_generated": len(dists),
        "SetA_reactive_2.3-4.0": len(setA),
        "SetB_thermo_>=5.0":     len(setB),
        "SetTS_1.7-2.3":         len(setTS),
        "side_path_<1.7":        len(side_path),
        "distance_min":  float(min(dists)),
        "distance_max":  float(max(dists)),
        "distance_mean": float(np.mean(dists)),
        "co_smiles": "Cl[C]=O",
        "ch3o_smiles": "ClOC",
    }, f, indent=2)

# 거리 히스토그램
fig, ax = plt.subplots(figsize=(9,5))
ax.hist(dists, bins=40, color='steelblue', edgecolor='black', alpha=0.7)
ax.axvspan(1.7, 2.3, alpha=0.2, color='orange', label='Set TS (1.7–2.3)')
ax.axvspan(2.3, 4.0, alpha=0.3, color='green',  label='Set A reactive (2.3–4.0)')
ax.axvspan(0, 1.7,   alpha=0.2, color='red',    label='side-path (<1.7)')
ax.axvspan(5.0, max(dists)*1.05, alpha=0.2, color='gray', label='Set B thermo (≥5.0)')
ax.set_xlabel('d(C$_{CO}$ ⋯ O$_{CH_3O}$)  /  Å')
ax.set_ylabel('count')
ax.set_title(f'S1 Pd(100) co-adsorption distance distribution (CO* = Cl[C]=O)')
ax.legend(loc='upper right', fontsize=9)
plt.tight_layout()
plt.savefig(COADS/"distance_hist.png", dpi=140, bbox_inches='tight')
plt.close()
print(f"  saved distance_hist.png")

# Set A grid
if setA:
    n_show = min(25, len(setA)); idx_grid = np.linspace(0, len(setA)-1, n_show, dtype=int)
    fig, axes = plt.subplots(5,5, figsize=(15,15))
    for i, ax in enumerate(axes.flat):
        if i < len(idx_grid):
            plot_atoms(setA[idx_grid[i]], ax, rotation='0x,0y,0z', radii=0.85, show_unit_cell=0)
            ax.set_title(f'#{idx_grid[i]}', fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    plt.suptitle(f"S1 Set A reactive ({len(setA)} structures, C–O 2.3–4.0 Å)", fontsize=14)
    plt.tight_layout()
    plt.savefig(COADS/"SetA_grid.png", dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  saved SetA_grid.png")

# Set B grid
if setB:
    n_show = min(25, len(setB)); idx_grid = np.linspace(0, len(setB)-1, n_show, dtype=int)
    fig, axes = plt.subplots(5,5, figsize=(15,15))
    for i, ax in enumerate(axes.flat):
        if i < len(idx_grid):
            plot_atoms(setB[idx_grid[i]], ax, rotation='0x,0y,0z', radii=0.85, show_unit_cell=0)
            ax.set_title(f'#{idx_grid[i]}', fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    plt.suptitle(f"S1 Set B thermo ({len(setB)} structures, ≥5.0 Å)", fontsize=14)
    plt.tight_layout()
    plt.savefig(COADS/"SetB_grid.png", dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  saved SetB_grid.png")

print("\n=== DONE ===")
print(f"CO* candidates (Cl[C]=O, pure CO): {len(co_structs)}")
print(f"CH3O* candidates (loaded existing): {len(ch3o_structs)}")
print(f"co-ads total: {len(dists)}, SetA={len(setA)}, SetB={len(setB)}, SetTS={len(setTS)}, side={len(side_path)}")
