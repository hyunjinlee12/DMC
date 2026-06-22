# Phase 1 + Phase 2 Ranked Structures Archive

Generated: see file mtime. Includes the **ranked relaxed structures** from the
MACE-MH + D3 ranking pipeline.

## Layout

```
phase1/
  S1_Pd100/
    summary.json                        # per-surface stats
    ranked_CO.json   unique_CO.json     # rank order, dedup'd
    ranked_CH3O.json unique_CH3O.json
    top10_CO.png     top10_CH3O.png     # distribution plot
    CO_structures/                       # extxyz of top-10, mid-5, bottom-5
      CO_rank0000_idxNNNNN.xyz
      ...
    CH3O_structures/
      ...
  S2_PdO101_Pd100/  ... (same layout)
  ...

phase2/
  S1_Pd100/
    summary.json
    ranked_SetA.json   unique_SetA.json
    shortlist_SetA.json shortlist_SetA.png
    SetA_structures/                     # top-10, mid-5, bottom-5 (extxyz)
    SetA_shortlist/                      # DFT shortlist (3-4 candidates)
  ...
  (S4 missing — Phase 2 still running)
```

## File metadata (in `info` dict of extxyz)

- `rank`: position in ranked_*.json (0 = lowest E)
- `E_MACE`: MACE+D3 energy (eV) — RANKING ONLY, not DFT ground truth
- `dE_rel_meV`: relative to rank-0 (in meV)
- `d_min`: nearest ads-substrate distance (Å)
- `d_reactive`: (Phase 2 only) C_CO ↔ O_CH3O distance (Å)
- `converged`: True/False (LBFGS fmax=0.05)

## Adsorbate atom order (last N atoms)

- CO*  (N=2): [C, O]
- CH3O* (N=5): [O, C, H, H, H]  (from SMILES `ClOC`)
- co-ads SetA (N=7): [C, O, O, C, H, H, H]  (CO followed by CH3O)

## Open with ASE

```python
from ase.io import read
atoms_list = read('phase1/S1_Pd100/CO_structures/CO_rank0000_idx00056.xyz')
# or all in one dir:
import glob; files = sorted(glob.glob('phase1/S1_Pd100/CO_structures/*.xyz'))
for f in files:
    a = read(f); print(f, a.info['rank'], a.info['dE_rel_meV'])
```

## Phase 2 status

S1, S2, S3, S3b — DONE.  S4 PdO2(110) — still running (~9 hr remaining as of archive).
