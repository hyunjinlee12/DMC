"""T1.15 — Build DFT shortlist from Phase 1 (single ads) + Phase 2 filtered (co-ads).

Per gate guide:
  - T1.16 = Level 1 DFT (vacuum), T1.17 = Level 2 (VASPsol)
  - 표면당 ~6-9 후보 (S2 5-8 권장)
  - distance-bin stratified for co-ads

Sources:
  Phase 1 D3 (single ads):  MLIP_phase1/unique_{CO|CH3O}.json
                            + relaxed_{CO|CH3O}.traj
  Phase 2 filtered (co-ads): MLIP_phase2_filtered/shortlist_SetA.json
                             + MLIP_phase2/relaxed_SetA.traj  (filtered indexes only)
  S4 manual placement:      Pd atop, O atop, bridge — manually constructed
                            (commit Phase 3 refilter showed S4 co-ads unviable)

Output structure:
  calculations/G3_adsorption/DFT_shortlist/
    {surface}/
      single_CO/
        00_CO_phase1_rank0.vasp        (POSCAR)
        ...
      single_CH3O/
        00_CH3O_phase1_rank0.vasp
        ...
      coads/
        00_SetA_filtered_bin1.vasp    (lowest-E per distance bin)
        ...
      shortlist_index.json              (metadata for all entries)
"""
import json, shutil
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT / 'calculations/G3_adsorption'
G2 = ROOT / 'calculations/G2_slab'
OUT = G3 / 'DFT_shortlist'
if OUT.exists():
    shutil.rmtree(OUT)
OUT.mkdir()

SURFACES = {
    'S1':  'S1_Pd100',
    'S2':  'S2_PdO101_Pd100',
    'S3':  'S3_PdO100',
    'S3b': 'S3b_PdO100_PdOterm',
    'S4':  'S4_PdO2_110',
}

# Per-surface DFT shortlist budget (gate guide §P1-D)
BUDGET = {
    'S1':  {'CO': 3, 'CH3O': 3, 'coads': 3},
    'S2':  {'CO': 5, 'CH3O': 5, 'coads': 4},   # S2 PdO/Pd interface — extra
    'S3':  {'CO': 3, 'CH3O': 3, 'coads': 3},
    'S3b': {'CO': 3, 'CH3O': 3, 'coads': 3},
    'S4':  {'CO': 3, 'CH3O': 3, 'coads': 0},   # S4 co-ads unviable → 0
}


def apply_fix_constraint(atoms, n_sub):
    """Fix bottom 50% of substrate atoms by z (matches DFT relax convention)."""
    z = atoms.positions[:n_sub, 2]
    z_med = float(np.median(z))
    fixed_indices = [i for i in range(n_sub) if atoms.positions[i, 2] < z_med]
    atoms.set_constraint(FixAtoms(indices=fixed_indices))
    return len(fixed_indices)


def save_poscar(atoms, path, label, n_fixed):
    """Write VASP POSCAR with selective dynamics."""
    write(str(path), atoms, format='vasp',
          direct=True, sort=True, vasp5=True)
    return label


def collect_single_ads(sid, sdir, ads, budget):
    """Pick top-budget from Phase 1 unique pool, geometric-diversity-aware."""
    src = G3 / sdir / 'MLIP_phase1'
    unique = json.load(open(src / f'unique_{ads}.json'))
    # Get representatives by d_min spread (avoid all-atop pile-up)
    # Sort by E, pick top-budget but require distinct d_min bins
    selected = []
    bins_used = set()
    for r in unique:
        d_bin = round(r.get('d_min', 0), 1)
        if d_bin not in bins_used or len(selected) < budget // 2:
            selected.append(r)
            bins_used.add(d_bin)
        if len(selected) >= budget:
            break
    # Fallback to top-N if not enough diversity
    if len(selected) < budget:
        for r in unique:
            if r not in selected:
                selected.append(r)
            if len(selected) >= budget:
                break

    # Load structures from relaxed traj (traj indexed by original candidate idx)
    traj = list(read(src / f'relaxed_{ads}.traj', index=':'))
    n_ads = 2 if ads == 'CO' else 5
    out_entries = []
    for k, rec in enumerate(selected[:budget]):
        orig_idx = rec['idx']
        if orig_idx >= len(traj):
            continue
        atoms = traj[orig_idx]
        out_entries.append({
            'kind': f'single_{ads}',
            'rank': k,
            'orig_idx': orig_idx,
            'E_MACE': rec['E'],
            'dE_rel_meV': rec.get('dE_rel_meV', 0),
            'd_min': rec.get('d_min', 0),
            'n_atoms': len(atoms),
            'n_ads': n_ads,
            'atoms_obj': atoms,
        })
    return out_entries


def collect_coads(sid, sdir, budget):
    """Pick from Phase 2 FILTERED shortlist (distance-bin stratified)."""
    if budget == 0:
        return []
    src_filtered = G3 / sdir / 'MLIP_phase2_filtered'
    src_orig = G3 / sdir / 'MLIP_phase2'
    if not (src_filtered / 'shortlist_SetA.json').exists():
        return []
    shortlist = json.load(open(src_filtered / 'shortlist_SetA.json'))
    if not shortlist:
        return []
    # Structures live in original Phase 2 relaxed traj (filtered just keeps valid indexes)
    traj = list(read(src_orig / 'relaxed_SetA.traj', index=':'))
    out_entries = []
    for k, rec in enumerate(shortlist[:budget]):
        orig_idx = rec['idx']
        if orig_idx >= len(traj):
            continue
        atoms = traj[orig_idx]
        out_entries.append({
            'kind': 'coads_SetA',
            'rank': k,
            'orig_idx': orig_idx,
            'E_MACE': rec['E'],
            'dE_rel_meV': rec.get('dE_rel_meV', 0),
            'd_min': rec.get('d_min', 0),
            'd_reactive': rec.get('d_reactive', 0),
            'n_atoms': len(atoms),
            'n_ads': 7,
            'atoms_obj': atoms,
        })
    return out_entries


def main():
    print('=== T1.15 DFT shortlist consolidation ===\n')
    global_index = []
    total_jobs = 0
    total_atoms = 0

    for sid, sdir in SURFACES.items():
        b = BUDGET[sid]
        out_surf = OUT / sid
        out_surf.mkdir()
        slab = read(G2 / sdir / 'CONTCAR')
        n_sub = len(slab)

        entries = []
        for ads in ['CO', 'CH3O']:
            entries.extend(collect_single_ads(sid, sdir, ads, b[ads]))
        entries.extend(collect_coads(sid, sdir, b['coads']))

        # Write POSCAR for each
        for e in entries:
            atoms = e['atoms_obj']
            n_fixed = apply_fix_constraint(atoms, n_sub)
            label = f"{e['kind']}_rank{e['rank']}_idx{e['orig_idx']:05d}"
            sub_dir = out_surf / e['kind']
            sub_dir.mkdir(exist_ok=True)
            path = sub_dir / f"{e['rank']:02d}_{label}.vasp"
            save_poscar(atoms, path, label, n_fixed)
            e_record = {k: v for k, v in e.items() if k != 'atoms_obj'}
            e_record['poscar'] = str(path.relative_to(ROOT))
            e_record['surface'] = sid
            e_record['n_fixed_atoms'] = n_fixed
            global_index.append(e_record)
            total_jobs += 1
            total_atoms += e['n_atoms']

        # Per-surface summary
        json.dump([{k: v for k, v in e.items() if k != 'atoms_obj'} for e in entries],
                  open(out_surf / 'shortlist.json', 'w'), indent=2)
        print(f'  {sid:<5} CO={sum(1 for e in entries if e["kind"]=="single_CO")} '
              f'CH3O={sum(1 for e in entries if e["kind"]=="single_CH3O")} '
              f'coads={sum(1 for e in entries if e["kind"]=="coads_SetA")}  '
              f'total={len(entries)}')

    json.dump(global_index, open(OUT / 'shortlist_global.json', 'w'), indent=2)
    print()
    print(f'=== Total: {total_jobs} DFT jobs ===')
    print(f'Average atoms/job: {total_atoms/total_jobs:.0f}')
    print(f'Output: {OUT}')

    # Cost estimate (PBE+D3 VASP relax, RTX 6000 Ada GPU)
    avg_atoms = total_atoms / total_jobs
    # ~6-12 hr per relax on this system size for VASP GPU
    hr_per_job_low = 6.0
    hr_per_job_high = 12.0
    total_hr_low = total_jobs * hr_per_job_low
    total_hr_high = total_jobs * hr_per_job_high
    # If we can run 2 jobs in parallel (2 GPUs):
    n_parallel = 2
    wall_hr_low = total_hr_low / n_parallel
    wall_hr_high = total_hr_high / n_parallel
    print()
    print('=== DFT cost estimate (T1.16 Level 1 vacuum) ===')
    print(f'  Per-job:     {hr_per_job_low:.0f}-{hr_per_job_high:.0f} hr (VASP+D3 on RTX 6000 Ada, ~{avg_atoms:.0f} atoms)')
    print(f'  Sequential:  {total_hr_low:.0f}-{total_hr_high:.0f} GPU-hr ({total_jobs} jobs)')
    print(f'  Parallel (2 GPUs): {wall_hr_low/24:.1f}-{wall_hr_high/24:.1f} days wall time')
    print()
    print('=== T1.17 Level 2 (VASPsol single-point) ===')
    sp_hr_low = total_jobs * 1.0
    sp_hr_high = total_jobs * 2.0
    print(f'  Per-job:     1-2 hr each')
    print(f'  Parallel (2 GPUs): {sp_hr_low/n_parallel/24:.1f}-{sp_hr_high/n_parallel/24:.1f} days')
    print()
    print('=== Total T1.16+1.17 ===')
    total_w_low = (total_hr_low + sp_hr_low) / n_parallel / 24
    total_w_high = (total_hr_high + sp_hr_high) / n_parallel / 24
    print(f'  {total_w_low:.1f} - {total_w_high:.1f} days wall time (2 GPUs)')


if __name__ == '__main__':
    main()
