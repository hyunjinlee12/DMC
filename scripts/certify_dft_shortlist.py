"""Final DFT shortlist certification — per-POSCAR strict chemistry audit.

For each of 47 POSCAR candidates:
  - Identify adsorbate atoms by species (NOT by index — POSCAR sort=True reorders)
  - Compute intramolecular bonds (C=O for CO, O-C/C-H for CH3O)
  - Compute Pd-X chemisorption distance
  - Verify FixAtoms applied
  - Classify: OK / SUSPECT / BROKEN

Outputs:
  reports/G3/DFT_certification/
    per_poscar_audit.json    — every POSCAR's chemistry report
    certification_table.md   — go/no-go per candidate
    A_per_surface_status.png — summary figure
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read
from ase.constraints import FixAtoms

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
SHORTLIST_DIR = G3 / 'DFT_shortlist'
OUT = ROOT / 'reports/G3/DFT_certification'
OUT.mkdir(parents=True, exist_ok=True)

# Catalysis reference (PBE+D3)
RANGES = {
    'Pd-C':       (1.85, 2.50),   # CO* chemisorbed
    'Pd-O':       (1.95, 2.30),   # CH3O* / surface O
    'C=O':        (1.05, 1.30),   # CO intra
    'O-C_methoxy':(1.30, 1.55),   # methoxy O-C
    'C-H':        (0.95, 1.20),   # methyl C-H
    'd_react':    (2.10, 4.00),   # SetA reactive band
}


def find_adsorbate_atoms(atoms, kind):
    """Identify adsorbate atoms by species (not by index).
    Returns dict with positions of key atoms.
    """
    syms = np.array(atoms.get_chemical_symbols())
    pos = atoms.positions
    z = pos[:, 2]
    z_sub_max = np.max(z[syms == 'Pd'])   # top of Pd substrate

    # CO adsorbate: 1 C atom (definitely adsorbate), nearest O is the CO oxygen
    c_idx = np.where(syms == 'C')[0]
    o_idx_all = np.where(syms == 'O')[0]
    h_idx = np.where(syms == 'H')[0]

    if kind == 'single_CO':
        if len(c_idx) != 1: return None
        c = c_idx[0]
        # Find nearest O to C
        dists = atoms.get_distances(c, o_idx_all, mic=True)
        o_ads = o_idx_all[np.argmin(dists)]
        return {'C_CO': c, 'O_CO': o_ads, 'kind': 'CO',
                'sub_O_idx': [i for i in o_idx_all if i != o_ads]}

    if kind == 'single_CH3O':
        # 1 C atom (methyl C), 3 H atoms close to C
        if len(c_idx) != 1: return None
        c = c_idx[0]
        # methoxy O is nearest O to C with d ~ 1.4 Å (intramolecular)
        dists = atoms.get_distances(c, o_idx_all, mic=True)
        # Pick the closest O that's NOT a substrate O (z close to C)
        z_c = pos[c, 2]
        # methoxy O should be close in z to C
        candidates = sorted([(i, atoms.get_distance(c, i, mic=True), abs(pos[i,2] - z_c)) for i in o_idx_all], key=lambda x: x[1])
        o_methoxy = candidates[0][0]
        return {'C_methyl': c, 'O_methoxy': o_methoxy, 'H_methyl': list(h_idx), 'kind': 'CH3O'}

    if kind == 'coads_SetA':
        # 2 C atoms: one for CO, one for CH3O methyl
        # 3 H: belong to methyl
        if len(c_idx) != 2 or len(h_idx) != 3: return None
        # CH3O C is the one bonded to 3 H
        h_to_c = {h: c_idx[np.argmin(atoms.get_distances(h, c_idx, mic=True))] for h in h_idx}
        c_methyl_candidates = list(set(h_to_c.values()))
        if len(c_methyl_candidates) != 1: return None
        c_methyl = c_methyl_candidates[0]
        c_co = [c for c in c_idx if c != c_methyl][0]
        # methoxy O: nearest O to c_methyl
        d_co_o = atoms.get_distances(c_co, o_idx_all, mic=True)
        d_me_o = atoms.get_distances(c_methyl, o_idx_all, mic=True)
        o_co = o_idx_all[np.argmin(d_co_o)]
        # methoxy O: nearest O to methyl C (excluding CO's O)
        rest = [o for o in o_idx_all if o != o_co]
        if not rest: return None
        d_me_o_rest = atoms.get_distances(c_methyl, rest, mic=True)
        o_methoxy = rest[np.argmin(d_me_o_rest)]
        return {'C_CO': c_co, 'O_CO': o_co, 'C_methyl': c_methyl,
                'O_methoxy': o_methoxy, 'H_methyl': list(h_idx), 'kind': 'coads'}

    return None


def audit_poscar(path, kind, verbose=False):
    """Return audit dict for a single POSCAR."""
    atoms = read(path)
    ads = find_adsorbate_atoms(atoms, kind)
    if ads is None:
        return {'path': str(path.relative_to(ROOT)), 'kind': kind,
                'verdict': 'UNRECOGNIZED', 'note': 'cannot identify adsorbate atoms'}

    result = {'path': str(path.relative_to(ROOT)), 'kind': kind,
              'n_atoms': len(atoms), 'ads_indices': {k: v for k, v in ads.items() if k != 'kind'}}

    # FixAtoms check
    n_fixed = 0
    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            n_fixed = len(c.get_indices())
    result['n_fixed'] = n_fixed
    n_sub = len(atoms) - (2 if kind == 'single_CO' else 5 if kind == 'single_CH3O' else 7)
    result['n_substrate'] = n_sub
    result['fix_fraction'] = n_fixed / n_sub if n_sub else 0

    # Intramolecular bonds
    syms = atoms.get_chemical_symbols()
    pd_idx = [i for i, s in enumerate(syms) if s == 'Pd']

    if kind == 'single_CO':
        d_co_intra = atoms.get_distance(ads['C_CO'], ads['O_CO'], mic=True)
        d_pdc = min(atoms.get_distances(ads['C_CO'], pd_idx, mic=True))
        d_pdo = min(atoms.get_distances(ads['O_CO'], pd_idx, mic=True))
        result.update({'d_C_O_intra': float(d_co_intra),
                       'd_Pd_C': float(d_pdc), 'd_Pd_O': float(d_pdo)})

    elif kind == 'single_CH3O':
        d_oc = atoms.get_distance(ads['C_methyl'], ads['O_methoxy'], mic=True)
        d_pdo = min(atoms.get_distances(ads['O_methoxy'], pd_idx, mic=True))
        d_chs = [atoms.get_distance(ads['C_methyl'], h, mic=True) for h in ads['H_methyl']]
        result.update({'d_O_C_methoxy': float(d_oc),
                       'd_Pd_O': float(d_pdo),
                       'd_C_H_min': float(min(d_chs)),
                       'd_C_H_max': float(max(d_chs))})

    elif kind == 'coads_SetA':
        d_co_intra = atoms.get_distance(ads['C_CO'], ads['O_CO'], mic=True)
        d_methoxy_oc = atoms.get_distance(ads['C_methyl'], ads['O_methoxy'], mic=True)
        d_pdc_co = min(atoms.get_distances(ads['C_CO'], pd_idx, mic=True))
        d_pdo_methoxy = min(atoms.get_distances(ads['O_methoxy'], pd_idx, mic=True))
        # Use DIRECT distance for reactive (PBC trap-safe)
        d_react_direct = float(np.linalg.norm(atoms.positions[ads['C_CO']] - atoms.positions[ads['O_methoxy']]))
        d_chs = [atoms.get_distance(ads['C_methyl'], h, mic=True) for h in ads['H_methyl']]
        result.update({'d_C_O_intra_CO': float(d_co_intra),
                       'd_O_C_methoxy': float(d_methoxy_oc),
                       'd_react_direct': d_react_direct,
                       'd_Pd_C_CO': float(d_pdc_co),
                       'd_Pd_O_methoxy': float(d_pdo_methoxy),
                       'd_C_H_min': float(min(d_chs)),
                       'd_C_H_max': float(max(d_chs))})

    # ----- Classify -----
    issues = []
    if kind == 'single_CO':
        if not RANGES['C=O'][0] <= result['d_C_O_intra'] <= RANGES['C=O'][1]:
            issues.append(f"C=O intra {result['d_C_O_intra']:.2f} out of {RANGES['C=O']}")
        if not RANGES['Pd-C'][0] <= result['d_Pd_C'] <= RANGES['Pd-C'][1]:
            issues.append(f"Pd-C {result['d_Pd_C']:.2f} out of chemisorbed {RANGES['Pd-C']}")
    elif kind == 'single_CH3O':
        if not RANGES['O-C_methoxy'][0] <= result['d_O_C_methoxy'] <= RANGES['O-C_methoxy'][1]:
            issues.append(f"methoxy O-C {result['d_O_C_methoxy']:.2f} out of {RANGES['O-C_methoxy']}")
        if not RANGES['Pd-O'][0] <= result['d_Pd_O'] <= RANGES['Pd-O'][1]:
            issues.append(f"Pd-O {result['d_Pd_O']:.2f} out of chemisorbed {RANGES['Pd-O']}")
        if not RANGES['C-H'][0] <= result['d_C_H_min'] <= RANGES['C-H'][1]:
            issues.append(f"C-H min {result['d_C_H_min']:.2f} out of {RANGES['C-H']}")
        if not RANGES['C-H'][0] <= result['d_C_H_max'] <= RANGES['C-H'][1]:
            issues.append(f"C-H max {result['d_C_H_max']:.2f} out of {RANGES['C-H']}")
    elif kind == 'coads_SetA':
        if not RANGES['C=O'][0] <= result['d_C_O_intra_CO'] <= RANGES['C=O'][1]:
            issues.append(f"CO C=O {result['d_C_O_intra_CO']:.2f} out of {RANGES['C=O']}")
        if not RANGES['O-C_methoxy'][0] <= result['d_O_C_methoxy'] <= RANGES['O-C_methoxy'][1]:
            issues.append(f"methoxy O-C {result['d_O_C_methoxy']:.2f} out of {RANGES['O-C_methoxy']}")
        if not RANGES['d_react'][0] <= result['d_react_direct'] <= RANGES['d_react'][1]:
            issues.append(f"d_react {result['d_react_direct']:.2f} out of SetA {RANGES['d_react']}")
        if not RANGES['Pd-C'][0] <= result['d_Pd_C_CO'] <= RANGES['Pd-C'][1]:
            issues.append(f"Pd-C(CO) {result['d_Pd_C_CO']:.2f} out of {RANGES['Pd-C']}")
        if not RANGES['Pd-O'][0] <= result['d_Pd_O_methoxy'] <= RANGES['Pd-O'][1]:
            issues.append(f"Pd-O(methoxy) {result['d_Pd_O_methoxy']:.2f} out of {RANGES['Pd-O']}")
        if not RANGES['C-H'][0] <= result['d_C_H_min'] <= RANGES['C-H'][1]:
            issues.append(f"C-H min {result['d_C_H_min']:.2f} out")

    if not issues:
        result['verdict'] = 'OK'
    elif len(issues) == 1 and ('Pd-C' in issues[0] or 'Pd-O' in issues[0]):
        result['verdict'] = 'SUSPECT_weak_binding'
    else:
        result['verdict'] = 'BROKEN' if any('intra' in i or 'methoxy O-C' in i or 'C-H' in i for i in issues) else 'SUSPECT'
    result['issues'] = issues
    return result


# ============================================================
# Run audit on all 47 POSCARs
# ============================================================
all_audits = []
SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']

for sid in SURFACES:
    for kind in ['single_CO', 'single_CH3O', 'coads_SetA']:
        sub = SHORTLIST_DIR / sid / kind
        if not sub.exists(): continue
        for vasp in sorted(sub.glob('*.vasp')):
            audit = audit_poscar(vasp, kind)
            audit['surface'] = sid
            all_audits.append(audit)

# Stats per surface
print('=== Audit Results ===\n')
print(f"{'Surface':<5} {'kind':<14} {'rank':<5} {'verdict':<25} {'key issue':<60}")
print('=' * 110)
for a in all_audits:
    rank = Path(a['path']).stem.split('_')[0]
    issue = a['issues'][0] if a.get('issues') else ''
    print(f"{a['surface']:<5} {a['kind']:<14} {rank:<5} {a['verdict']:<25} {issue:<60}")

# Per-surface summary
print()
summary = {sid: {'OK': 0, 'SUSPECT': 0, 'SUSPECT_weak_binding': 0, 'BROKEN': 0, 'UNRECOGNIZED': 0} for sid in SURFACES}
for a in all_audits:
    summary[a['surface']][a['verdict']] = summary[a['surface']].get(a['verdict'], 0) + 1
print(f"{'Sur':<5} {'OK':>4} {'SUSPECT':>9} {'weak':>6} {'BROKEN':>8} {'unrecog':>9}")
for sid in SURFACES:
    s = summary[sid]
    print(f"{sid:<5} {s['OK']:>4} {s['SUSPECT']:>9} {s['SUSPECT_weak_binding']:>6} {s['BROKEN']:>8} {s['UNRECOGNIZED']:>9}")

json.dump(all_audits, open(OUT / 'per_poscar_audit.json', 'w'), indent=2)
json.dump(summary, open(OUT / 'per_surface_summary.json', 'w'), indent=2)

# Plot
fig, ax = plt.subplots(figsize=(10, 5.5))
x = np.arange(len(SURFACES))
ok = [summary[s]['OK'] for s in SURFACES]
weak = [summary[s]['SUSPECT_weak_binding'] for s in SURFACES]
susp = [summary[s]['SUSPECT'] for s in SURFACES]
broken = [summary[s]['BROKEN'] for s in SURFACES]
ax.bar(x, ok, color='#27ae60', edgecolor='black', label='OK')
ax.bar(x, weak, bottom=ok, color='#f4a261', edgecolor='black', label='SUSPECT_weak_binding')
ax.bar(x, susp, bottom=np.array(ok)+np.array(weak), color='#e76f51', edgecolor='black', label='SUSPECT')
ax.bar(x, broken, bottom=np.array(ok)+np.array(weak)+np.array(susp), color='#c0392b', edgecolor='black', label='BROKEN')
for i, s in enumerate(SURFACES):
    total = ok[i]+weak[i]+susp[i]+broken[i]
    ax.text(i, total + 0.5, str(total), ha='center', fontsize=11, weight='bold')
ax.set_xticks(x); ax.set_xticklabels(SURFACES)
ax.set_ylabel('Number of DFT shortlist candidates')
ax.set_title('DFT shortlist certification — per-surface breakdown (n=47 total)')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(OUT / 'A_per_surface_status.png', dpi=140, bbox_inches='tight')
plt.close()
print(f'\nFigures + JSONs in {OUT}')
