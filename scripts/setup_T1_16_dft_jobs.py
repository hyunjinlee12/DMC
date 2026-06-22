"""Set up T1.16 Level 1 (vacuum) DFT jobs for 47 shortlist candidates.

For each candidate POSCAR in DFT_shortlist/{surface}/{kind}/:
  1. Create job directory: calculations/T1_16_DFT_L1/{surface}/{kind}/{candidate_name}/
  2. Copy POSCAR there
  3. Copy INCAR from G2 slab same surface (ISMEAR/SIGMA matches surface)
  4. Build POTCAR by cat'ing element POTCARs in same order as POSCAR
  5. Symlink submit_vasp_gpu.sh

DOES NOT submit. User must explicitly sbatch after review.

Output:
  calculations/T1_16_DFT_L1/                    # main calc dir
    submit_all.sh                                # batch submit (commented sbatch lines)
    job_index.json                               # full inventory
"""
import json
import shutil
from pathlib import Path
from ase.io import read

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT / 'calculations/G2_slab'
SHORTLIST = ROOT / 'calculations/G3_adsorption/DFT_shortlist'
T1_16 = ROOT / 'calculations/T1_16_DFT_L1'
POTCAR_LIB = Path('/home/hyunjin/POTENTIAL/potpaw_PBE')
SUBMIT_SCRIPT = ROOT / 'scripts/submit_vasp_gpu.sh'

SURFACES = ['S1', 'S2', 'S3', 'S3b', 'S4']
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

# Element -> POTCAR path
ELEM_POTCAR = {
    'Pd': POTCAR_LIB / 'Pd_pv/POTCAR',
    'O':  POTCAR_LIB / 'O/POTCAR',
    'C':  POTCAR_LIB / 'C/POTCAR',
    'H':  POTCAR_LIB / 'H/POTCAR',
}


def build_potcar(poscar_path, out_potcar):
    """Cat element POTCARs in the order POSCAR specifies."""
    atoms = read(poscar_path)
    # Get unique elements in order of first appearance
    syms = atoms.get_chemical_symbols()
    seen = []
    for s in syms:
        if s not in seen:
            seen.append(s)
    with open(out_potcar, 'wb') as fout:
        for elem in seen:
            potcar_src = ELEM_POTCAR[elem]
            if not potcar_src.exists():
                raise FileNotFoundError(f'POTCAR not found for {elem}: {potcar_src}')
            with open(potcar_src, 'rb') as fin:
                fout.write(fin.read())
    return seen


def main():
    if T1_16.exists():
        print(f'{T1_16} already exists — refusing to overwrite. Remove first if redo.')
        return
    T1_16.mkdir(parents=True)

    job_index = []
    submit_lines = ['#!/bin/bash', '# T1.16 DFT batch submit', '# Review then uncomment desired lines',
                    '# All paths relative to project root', '', 'cd /home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc', '']

    for sid in SURFACES:
        sdir = SDIRS[sid]
        src_incar = G2 / sdir / 'INCAR'
        if not src_incar.exists():
            print(f'[{sid}] missing INCAR — skip'); continue
        for kind in ['single_CO', 'single_CH3O', 'coads_SetA']:
            sub = SHORTLIST / sid / kind
            if not sub.exists(): continue
            for vasp in sorted(sub.glob('*.vasp')):
                if '.broken_bak' in vasp.name: continue
                job_name = vasp.stem   # e.g., "00_single_CO_rank0_idx00059"
                jdir = T1_16 / sid / kind / job_name
                jdir.mkdir(parents=True)
                # 1. POSCAR
                shutil.copy(vasp, jdir / 'POSCAR')
                # 2. INCAR
                shutil.copy(src_incar, jdir / 'INCAR')
                # 3. POTCAR
                elems = build_potcar(vasp, jdir / 'POTCAR')
                # 4. SLURM submit script (symlink)
                (jdir / 'submit_vasp_gpu.sh').symlink_to(SUBMIT_SCRIPT)
                jname = f'{sid}_{kind}_{job_name.split("_")[0]}'   # e.g., S1_single_CO_00
                submit_lines.append(f'# sbatch -J {jname} --chdir={jdir.relative_to(ROOT)} scripts/submit_vasp_gpu.sh')
                job_index.append({
                    'job_name': jname,
                    'surface': sid,
                    'kind': kind,
                    'jdir': str(jdir.relative_to(ROOT)),
                    'n_atoms': len(read(vasp)),
                    'elements': elems,
                })
    submit_lines.append('')
    submit_lines.append('echo "Done queueing $(grep -c sbatch ${0}) jobs (commented). Uncomment to submit."')

    (T1_16 / 'submit_all.sh').write_text('\n'.join(submit_lines))
    json.dump(job_index, open(T1_16 / 'job_index.json', 'w'), indent=2)
    print(f'\n✓ Set up {len(job_index)} DFT jobs in {T1_16}')

    # Quick verification: list per-surface job count
    from collections import Counter
    cnt = Counter((j['surface'], j['kind']) for j in job_index)
    print('\nPer-surface/kind:')
    for sid in SURFACES:
        for kind in ['single_CO', 'single_CH3O', 'coads_SetA']:
            c = cnt.get((sid, kind), 0)
            print(f'  {sid} {kind:<14}  {c}')
    print(f'  Total: {len(job_index)}')

    # Sanity: open one INCAR, one POTCAR
    sample = T1_16 / SURFACES[0] / 'single_CO'
    sample_jdir = next(sample.iterdir())
    print(f'\nSample job dir: {sample_jdir.relative_to(ROOT)}')
    print(f'  POSCAR exists: {(sample_jdir/"POSCAR").exists()}')
    print(f'  INCAR exists: {(sample_jdir/"INCAR").exists()}')
    print(f'  POTCAR exists: {(sample_jdir/"POTCAR").exists()}, size={(sample_jdir/"POTCAR").stat().st_size//1024} KB')
    print(f'  submit_vasp_gpu.sh symlink: {(sample_jdir/"submit_vasp_gpu.sh").is_symlink()}')


if __name__ == '__main__':
    main()
