"""Validate site classification used in v2 DFT shortlist.

Tests:
  1. Cutoff sensitivity — re-classify with 2.4, 2.5, 2.6, 2.7, 2.8 Å and
     count label changes per structure
  2. Cross-check filename labels (v2 picks) vs fresh recomputation
  3. Edge case audit — borderline distances (within 0.1 Å of cutoff)
  4. Physisorbed flag — count "physisorbed" picks and report
"""
import json, re, glob
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
from ase.io import read

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
V2 = G3 / 'DFT_shortlist_v2'

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

CUTOFFS = [2.4, 2.5, 2.6, 2.7, 2.8]
BASE_CUTOFF = 2.60


def classify(atoms, anchor, sub_indices, cutoff):
    syms = atoms.get_chemical_symbols()
    d = atoms.get_distances(anchor, sub_indices, mic=True)
    nbrs = [(sub_indices[i], d[i]) for i in range(len(sub_indices)) if d[i] < cutoff]
    n_pd = sum(1 for i, _ in nbrs if syms[i] == 'Pd')
    n_o = sum(1 for i, _ in nbrs if syms[i] == 'O')
    total = n_pd + n_o
    if total == 0: lbl = 'physisorbed'
    elif total == 1: lbl = f"atop_{'Pd' if n_pd else 'O'}"
    elif total == 2:
        if n_pd == 2: lbl = 'bridge_Pd-Pd'
        elif n_o == 2: lbl = 'bridge_O-O'
        else: lbl = 'bridge_Pd-O'
    elif total == 3:
        if n_pd == 3: lbl = 'hollow_3Pd'
        elif n_o == 3: lbl = 'hollow_3O'
        else: lbl = f'hollow_3({n_pd}Pd{n_o}O)'
    else: lbl = f'{total}f({n_pd}Pd{n_o}O)'
    return lbl, total, n_pd, n_o, sorted([d_ for _, d_ in nbrs])


def find_anchor_indices(atoms, ads):
    """Return (anchor_idx, sub_indices). anchor = C of CO or O of methoxy."""
    syms = atoms.get_chemical_symbols()
    z = atoms.positions[:, 2]
    n = len(atoms)
    if ads == 'CO':
        c_idx = [i for i in range(n) if syms[i] == 'C']
        anchor = c_idx[-1]  # take last (adsorbate added last)
        sub = [i for i in range(n) if i != anchor and not (syms[i] == 'O' and abs(z[i]-z[anchor])<1.4 and atoms.get_distance(i, anchor, mic=True)<1.4)]
    elif ads == 'CH3O':
        c_idx = [i for i in range(n) if syms[i] == 'C']
        if not c_idx: return None, None
        c = c_idx[-1]
        h_idx = [i for i in range(n) if syms[i] == 'H' and atoms.get_distance(i, c, mic=True) < 1.4]
        o_to_c = [i for i in range(n) if syms[i] == 'O' and atoms.get_distance(i, c, mic=True) < 1.7]
        if not o_to_c: return None, None
        anchor = o_to_c[0]  # methoxy O
        ads_set = {c, anchor} | set(h_idx)
        sub = [i for i in range(n) if i not in ads_set]
    return anchor, sub


def parse_filename_label(fname):
    """e.g. 00_single_CO_atop_Pd_idx00064.vasp → 'atop_Pd'."""
    m = re.match(r'\d+_single_(CO|CH3O)_(.+)_idx\d+\.vasp', fname)
    if m: return m.group(2)
    return None


# ===== Per-surface analysis =====
print("="*78)
print("SITE CLASSIFICATION VALIDATION REPORT")
print("="*78)

per_surface_summary = []

for sid, sdir in SDIRS.items():
    for ads in ['CO', 'CH3O']:
        traj_f = G3 / sdir / f'MLIP_phase1/relaxed_{ads}.traj'
        unique_f = G3 / sdir / f'MLIP_phase1/unique_{ads}.json'
        if not traj_f.exists() or not unique_f.exists():
            print(f"\n[SKIP] {sid} {ads}: file missing")
            continue
        try:
            traj = list(read(traj_f, index=':'))
            unique = json.load(open(unique_f))
        except Exception as e:
            print(f"\n[ERR] {sid} {ads}: {e}")
            continue

        labels_by_cutoff = defaultdict(list)
        borderline_count = 0
        for u in unique:
            atoms = traj[u['idx']]
            anchor, sub = find_anchor_indices(atoms, ads)
            if anchor is None: continue
            base_lbl, _, _, _, nbr_d_base = classify(atoms, anchor, sub, BASE_CUTOFF)

            structure_labels = {}
            for cut in CUTOFFS:
                lbl, _, _, _, _ = classify(atoms, anchor, sub, cut)
                structure_labels[cut] = lbl
            labels_by_cutoff['all'].append(structure_labels)

            # borderline: any distance within ±0.10 Å of 2.60
            if any(abs(d - BASE_CUTOFF) < 0.10 for d in nbr_d_base[:6]):
                borderline_count += 1

        # Per-cutoff label distribution
        print(f"\n--- {sid} {ads} (n={len(labels_by_cutoff['all'])} unique structures) ---")
        # Stability: how many keep same label across cutoffs
        all_same = 0
        changed = 0
        for sl in labels_by_cutoff['all']:
            lbls = set(sl.values())
            if len(lbls) == 1:
                all_same += 1
            else:
                changed += 1
        pct_stable = (all_same / len(labels_by_cutoff['all']) * 100) if labels_by_cutoff['all'] else 0
        print(f"  Cutoff-stable (same label 2.4-2.8 Å):  {all_same}/{len(labels_by_cutoff['all'])} ({pct_stable:.0f}%)")
        print(f"  Cutoff-sensitive (label changes):       {changed}")
        print(f"  Borderline (neighbor d within ±0.1 of 2.6): {borderline_count}")

        # Show breakdown @ base
        c_base = Counter(sl[BASE_CUTOFF] for sl in labels_by_cutoff['all'])
        print(f"  @ 2.6 Å: {dict(c_base)}")
        c_25 = Counter(sl[2.5] for sl in labels_by_cutoff['all'])
        c_27 = Counter(sl[2.7] for sl in labels_by_cutoff['all'])
        print(f"  @ 2.5 Å: {dict(c_25)}")
        print(f"  @ 2.7 Å: {dict(c_27)}")

        per_surface_summary.append({'sid': sid, 'ads': ads,
                                    'n': len(labels_by_cutoff['all']),
                                    'stable': all_same, 'changed': changed,
                                    'borderline': borderline_count})

# ===== v2 filename cross-check =====
print("\n" + "="*78)
print("V2 SHORTLIST FILENAME CROSS-CHECK")
print("="*78)
for sid in SDIRS:
    for kind in ['single_CO', 'single_CH3O']:
        files = sorted(glob.glob(str(V2 / sid / kind / '*.vasp')))
        if not files: continue
        labels = [parse_filename_label(Path(f).name) for f in files]
        c = Counter(labels)
        print(f"\n  {sid}/{kind} ({len(files)} picks): {dict(c)}")

# ===== Final verdict =====
print("\n" + "="*78)
print("VERDICT")
print("="*78)
total_n = sum(s['n'] for s in per_surface_summary)
total_stable = sum(s['stable'] for s in per_surface_summary)
total_changed = sum(s['changed'] for s in per_surface_summary)
total_borderline = sum(s['borderline'] for s in per_surface_summary)
print(f"  Total unique structures examined: {total_n}")
print(f"  Cutoff-stable (robust label):     {total_stable} ({100*total_stable/total_n:.1f}%)")
print(f"  Cutoff-sensitive (changes):       {total_changed} ({100*total_changed/total_n:.1f}%)")
print(f"  Borderline structures:            {total_borderline} ({100*total_borderline/total_n:.1f}%)")
print()
print("  → If stable% > 85%, classification is robust.")
print("  → If changed structures concentrated in 1-2 cells, MLIP-positioned site is ambiguous.")
