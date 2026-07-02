"""Verify v2 shortlist GEOMETRIC diversity (not just label diversity).

For each pick, compute anchor xy in fractional coords (mod 1).
Two picks at the same (a,b) frac → SAME site (regardless of label).
Two picks at different (a,b) frac → DIFFERENT site (regardless of label).

Distance threshold: 1.0 Å in projected xy (well below any nearest-neighbor in unit cell).
"""
import glob, re
from pathlib import Path
import numpy as np
from ase.io import read

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
V2 = ROOT / 'calculations/G3_adsorption/DFT_shortlist_v2'

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']


def find_anchor(atoms, kind):
    syms = atoms.get_chemical_symbols()
    if kind == 'CO':
        c_idx = [i for i, s in enumerate(syms) if s == 'C']
        return c_idx[-1]
    elif kind == 'CH3O':
        c_idx = [i for i, s in enumerate(syms) if s == 'C']
        c = c_idx[-1]
        o_to_c = [i for i, s in enumerate(syms) if s == 'O' and
                  atoms.get_distance(i, c, mic=True) < 1.7]
        return o_to_c[0]
    return None


def parse_label(fname):
    name = Path(fname).stem
    m = re.match(r'\d+_single_(CO|CH3O)_(.+)_idx\d+', name)
    if m: return m.group(2)
    m = re.match(r'\d+_coads_CO-(.+?)_OMe-(.+?)_d', name)
    if m: return f'CO[{m.group(1)}]+OMe[{m.group(2)}]'
    return '?'


def projected_xy(atoms, anchor_idx):
    """Project anchor xy onto unit cell, returning fractional (a, b)."""
    cell = atoms.cell.array
    pos = atoms.positions[anchor_idx]
    inv = np.linalg.inv(cell[:2, :2])
    frac = inv @ pos[:2]
    return frac % 1.0   # wrap into [0, 1)


def min_image_xy_dist(frac1, frac2, cell):
    """Minimum image distance in xy between two fractional positions."""
    dfrac = frac1 - frac2
    dfrac -= np.round(dfrac)   # MIC in fractional
    a, b = cell[0][:2], cell[1][:2]
    dvec = dfrac[0]*a + dfrac[1]*b
    return float(np.linalg.norm(dvec))


print('='*80)
print('SHORTLIST GEOMETRIC DIVERSITY AUDIT')
print('='*80)

total_hidden_dup = 0
total_checks = 0

for sid in SURFACES:
    for kind, subdir in [('CO', 'single_CO'), ('CH3O', 'single_CH3O')]:
        files = sorted(glob.glob(str(V2 / sid / subdir / '*.vasp')))
        if len(files) < 2: continue
        print(f'\n--- {sid} / {subdir} ({len(files)} picks) ---')
        recs = []
        for f in files:
            atoms = read(f)
            anchor = find_anchor(atoms, kind)
            frac = projected_xy(atoms, anchor)
            lbl = parse_label(f)
            recs.append({'file': Path(f).name, 'anchor_frac': frac,
                         'label': lbl, 'cell': atoms.cell.array})

        # Pairwise distance check (in projected xy, MIC)
        cell = recs[0]['cell']
        print(f"  {'pair':<25} {'label A':<22} {'label B':<22} {'xy_dist (Å)':<12} {'verdict':<15}")
        for i in range(len(recs)):
            for j in range(i+1, len(recs)):
                d = min_image_xy_dist(recs[i]['anchor_frac'],
                                       recs[j]['anchor_frac'], cell)
                # If d < 1.0 Å → effectively same site
                same_site = d < 1.0
                same_label = recs[i]['label'] == recs[j]['label']
                verdict = ''
                if same_site and not same_label:
                    verdict = 'FALSE-DIVERSITY'   # diff label but same site
                    total_hidden_dup += 1
                elif not same_site and same_label:
                    verdict = 'real-diversity-OK'   # same label but actually diff site (fine)
                elif same_site and same_label:
                    verdict = 'genuine-duplicate'   # both label + site same → bad
                else:
                    verdict = 'truly-diverse'      # diff label + diff site → good
                total_checks += 1
                print(f"  {i:02d}-{j:02d}                      {recs[i]['label'][:20]:<22} {recs[j]['label'][:20]:<22} {d:<12.2f} {verdict:<15}")

print('\n' + '='*80)
print(f'Total pairwise checks: {total_checks}')
print(f'False-diversity (different label, same site): {total_hidden_dup}')
print('='*80)
