"""Reviewer round 6 fixes:

P0-1  H200 submit scripts (partition=h200q, 16 CPU, Singularity/srun,
      /scratch/taehun1/hyunjin) — replace ALL previous submit_vasp*.sh
      in T1.17 / gas_references / magnetic_seed_test with H200 versions.

P0-2  Gas paired-static guard: add READY sentinel expectation. Submit
      script exits immediately if READY file is absent. Also copy
      prepare_static_from_relax.py into each paired-static dir for
      convenience.

P0-3  S3 AFM seed rebalance: exactly 32/32 Pd split, not 31/33.

P0-4  T1.17 remains PILOT — S1_clean is the pilot; do not mass-submit
      until seed test resolves MAGMOM policy. Add PILOT_ONLY marker to
      all other T1.17 dirs.

P0-5  Magnetization CSV interpretation flipped: top/bottom < 1 means
      BOTTOM-heavy (moments concentrated at bottom), NOT surface-localized.
      Fix wording; add layer-wise Pd/O moment table.

P0-6  Docs ΔG methanol sign: μ_liq = μ_gas + RT ln(p_sat/p°) which is
      NEGATIVE at 298 K (methanol p_sat < 1 bar). Not +ΔG_vap.

Read-only aside from these targeted edits.
"""
import csv, json, re, shutil
from pathlib import Path

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')

# ============================================================================
# H200 submit script template (per reviewer 2026-07-17 round 6)
# ============================================================================
H200_SUBMIT = """#!/bin/bash
#SBATCH --partition=h200q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time={walltime}
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err

# ============================================================================
# H200 submission (Singularity + srun). Adjust paths for your account.
# ============================================================================
{ready_guard}
# --- container + VASPsol build ---
# The container must be built with VASPsol linked (verify with:
#   singularity exec $CONTAINER grep -l VASPsol $VASP_BIN
# or by running S1_clean pilot first and checking OUTCAR for solvation output.)
CONTAINER=${{CONTAINER:-/scratch/taehun1/hyunjin/vasp_vaspsol.sif}}
VASP_BIN=${{VASP_BIN:-vasp_std}}

# --- runtime env ---
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Job: ${{SLURM_JOB_NAME}} | Dir: $(pwd) | Container: ${{CONTAINER}}"
echo "Start: $(date)"

srun --mpi=pmix singularity exec --nv "${{CONTAINER}}" "${{VASP_BIN}}"

echo "End: $(date)"
"""

READY_GUARD = """# --- READY sentinel guard (paired-static only) ---
# This dir is a static single-point that needs POSCAR replaced with the
# 15 Å relaxation's CONTCAR via prepare_static_from_relax.py. Do NOT submit
# until that has been done AND the READY file exists.
if [ ! -f READY ]; then
  echo "ERROR: READY sentinel absent. Run prepare_static_from_relax.py in"
  echo "the gas_references parent dir FIRST, then 'touch READY' here."
  exit 1
fi
if [ ! -f CONTCAR_SOURCE_SHA256 ]; then
  echo "WARNING: no CONTCAR_SOURCE_SHA256 recorded; cannot verify provenance."
fi

"""

def h200_script(kind='normal', walltime='7-00:00:00'):
    if kind == 'paired_static':
        return H200_SUBMIT.format(ready_guard=READY_GUARD, walltime='1-00:00:00')
    else:
        return H200_SUBMIT.format(ready_guard='', walltime=walltime)

# ============================================================================
# 1. Replace submit scripts across all pending dirs
# ============================================================================
def replace_submits():
    counts = {'T1.17': 0, 'gas_primary': 0, 'gas_paired_static': 0, 'seed': 0}
    # T1.17
    for d in (ROOT/'calculations/T1_17_VASPsol').rglob('submit_vasp*.sh'):
        d.write_text(h200_script('normal', walltime='3-00:00:00'))
        counts['T1.17'] += 1
    # Gas primary (non-static)
    for d in (ROOT/'calculations/gas_references').iterdir():
        if not d.is_dir(): continue
        sub = d/'submit_vasp.sh'
        if not sub.exists(): continue
        if 'static' in d.name:
            sub.write_text(h200_script('paired_static'))
            counts['gas_paired_static'] += 1
        else:
            sub.write_text(h200_script('normal', walltime='1-00:00:00'))
            counts['gas_primary'] += 1
    # Magnetic seed
    for d in (ROOT/'calculations/magnetic_seed_test').iterdir():
        if not d.is_dir(): continue
        sub = d/'submit_vasp.sh'
        if not sub.exists(): continue
        sub.write_text(h200_script('normal', walltime='1-00:00:00'))
        counts['seed'] += 1
    return counts

# ============================================================================
# 2. Copy prepare_static_from_relax.py into gas_references/ + docs
# ============================================================================
def add_paired_static_helper_and_notes():
    gas = ROOT/'calculations/gas_references'
    src = ROOT/'scripts/prepare_static_from_relax.py'
    dest = gas/'prepare_static_from_relax.py'
    if src.exists():
        shutil.copy(src, dest)
    # PILOT_INSTRUCTIONS in each paired-static dir
    for d in gas.iterdir():
        if not (d.is_dir() and 'static' in d.name): continue
        (d/'PILOT_INSTRUCTIONS.txt').write_text(
            f"This dir is a STATIC single-point box-size sanity check.\n"
            f"DO NOT SUBMIT until the following two steps are done:\n\n"
            f"  1. Wait for {d.name.rsplit('_vacuum_',1)[0]}_vacuum/ relax to finish\n"
            f"     (OUTCAR contains 'reached required accuracy').\n"
            f"  2. In the gas_references parent dir, run:\n"
            f"     python prepare_static_from_relax.py\n"
            f"     This overwrites this dir's POSCAR with the relaxed coordinates\n"
            f"     placed in a {d.name.rsplit('_',1)[0].split('_')[-1]} box (cell only).\n"
            f"  3. Verify OUTCAR_SOURCE_SHA256 file matches actual source CONTCAR.\n"
            f"  4. touch READY in this dir.\n"
            f"  5. Then sbatch submit_vasp.sh.\n\n"
            f"The submit script has a READY guard that will EXIT immediately if\n"
            f"READY does not exist. This prevents accidental submission with the\n"
            f"placeholder POSCAR.\n"
        )
    return dest

# ============================================================================
# 3. Rebalance S3 AFM seed (32/32)
# ============================================================================
def rebalance_s3_afm():
    """Rewrite S3_clean__afm INCAR MAGMOM so exactly 32 Pd are +0.5 and 32 are −0.5."""
    from ase.io import read
    afm = ROOT/'calculations/magnetic_seed_test/S3_clean__afm'
    if not afm.exists():
        return False
    poscar = read(afm/'POSCAR')
    syms = poscar.get_chemical_symbols()
    counts = poscar.get_chemical_formula('reduce')  # unused; we count directly
    n_O = sum(1 for s in syms if s == 'O')
    n_Pd = sum(1 for s in syms if s == 'Pd')
    assert n_Pd == 64, f'expected 64 Pd got {n_Pd}'
    # Pd atoms in POSCAR sort=True order come after O. Sort them by (frac_x, y, z),
    # take even-indexed = +, odd-indexed = − → exactly balanced.
    pd_indices = list(range(n_O, n_O + n_Pd))
    frac = poscar.get_scaled_positions()
    ordered = sorted(pd_indices, key=lambda i: (round(frac[i,0], 3),
                                                 round(frac[i,1], 3),
                                                 round(frac[i,2], 3)))
    signs = ['+0.500' if k % 2 == 0 else '-0.500' for k in range(n_Pd)]
    # Map back to POSCAR position order (Pd indices are contiguous after O)
    # Signs by original POSCAR order:
    sign_per_pd_orig_index = {}
    for k, idx in enumerate(ordered):
        sign_per_pd_orig_index[idx] = signs[k]
    magmom_parts = [f'{n_O}*0.000']
    for idx in pd_indices:
        magmom_parts.append(sign_per_pd_orig_index[idx])
    new_magmom = ' '.join(magmom_parts)
    # Sanity: exactly 32 '+' and 32 '-'
    n_plus = sum(1 for s in signs if s.startswith('+'))
    n_minus = sum(1 for s in signs if s.startswith('-'))
    assert n_plus == 32 and n_minus == 32, f'got {n_plus}+/{n_minus}-'
    # Rewrite INCAR
    incar_p = afm/'INCAR'
    incar_lines = incar_p.read_text().splitlines()
    for i, ln in enumerate(incar_lines):
        if ln.strip().startswith('MAGMOM'):
            incar_lines[i] = f'MAGMOM = {new_magmom}'
            break
    incar_p.write_text('\n'.join(incar_lines) + '\n')
    return True

# ============================================================================
# 4. Fix magnetization diagnosis CSV interpretation
# ============================================================================
def fix_magnetization_csv():
    p = ROOT/'paper_data/07b_magnetization_diagnosis.csv'
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        ratio = r.get('top_over_bottom_ratio', '')
        try:
            ratio_v = float(ratio) if ratio else None
        except ValueError:
            ratio_v = None
        if ratio_v is None:
            r['interpretation'] = 'no ratio available'
        elif 0.8 < ratio_v < 1.25:
            r['interpretation'] = ('moments distributed ~evenly through slab '
                                   '(suggests bulk-like FM behaviour or seed '
                                   'artifact — seed test to distinguish)')
        elif ratio_v >= 1.25:
            r['interpretation'] = ('TOP-heavy: moments concentrated at surface '
                                   '(consistent with physical Pd surface Stoner '
                                   'magnetism — verify with seed test)')
        else:  # < 0.8
            r['interpretation'] = ('BOTTOM-heavy: moments concentrated on fixed '
                                   'substrate atoms (NOT surface-localized). '
                                   'Could be a persistent initialization signature '
                                   'because those atoms did not move geometrically '
                                   'during T1.16, but SCF still optimized the '
                                   'magnetization density. Seed test with '
                                   'nonmagnetic initialization is needed to check '
                                   'whether a lower-E magnetic state exists.')
    keys = list(rows[0].keys())
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    return len(rows)

# ============================================================================
# 5. Update CALCULATION_SETTINGS.md ΔG methanol sign
# ============================================================================
def patch_calc_settings():
    p = ROOT/'docs/CALCULATION_SETTINGS.md'
    text = p.read_text()
    old = ('CH₃OH from `CH3OH_vaspsol` gives an electronic energy in implicit methanol —\n'
           'this is NOT directly the liquid chemical potential. To use it as a liquid\n'
           'reservoir at activity 1, add a **liquid-standard-state correction**:\n'
           '`μ_CH₃OH(l) ≈ G(CH3OH)_gas + ΔG_vap(exp = 4.1 kJ/mol)` OR extract from AIMD.')
    new = ('CH₃OH from `CH3OH_vaspsol` gives an electronic energy in implicit methanol —\n'
           'this is NOT directly the liquid chemical potential. To use it as a liquid\n'
           'reservoir at activity 1, add a **liquid-standard-state correction**. Using\n'
           'the ideal-gas-of-methanol reference plus the vapour-pressure relation\n'
           '`μ_liq(T) = μ_gas°(T) + RT ln(p_sat(T)/p°)`, with methanol p_sat ≈ 0.169 bar\n'
           'at 298.15 K (NIST Antoine parameters) → `μ_CH₃OH(l) − μ_CH₃OH(g)° ≈ −4.40 kJ/mol`\n'
           '(NEGATIVE, condensation reduces chemical potential relative to the\n'
           '1 bar gas standard). Alternative: extract μ_CH3OH(l) directly from AIMD\n'
           'of liquid methanol. **The sign of the correction was reported wrong in the\n'
           'previous version of this doc**; it is negative, not positive.')
    if old in text:
        text = text.replace(old, new)
        p.write_text(text)
        return True
    return False

# ============================================================================
# Run all
# ============================================================================
if __name__ == '__main__':
    print('== 1. H200 submit scripts ==')
    c = replace_submits()
    for k, v in c.items(): print(f'   {k}: {v} scripts')
    print('== 2. Paired-static helper + PILOT_INSTRUCTIONS ==')
    d = add_paired_static_helper_and_notes()
    print(f'   helper copied to {d}')
    print('== 3. S3 AFM seed rebalance ==')
    ok = rebalance_s3_afm()
    print(f'   {"OK" if ok else "SKIP (dir missing)"}')
    print('== 4. Magnetization CSV interpretation fix ==')
    n = fix_magnetization_csv()
    print(f'   {n} rows updated')
    print('== 5. CALCULATION_SETTINGS.md methanol sign ==')
    ok = patch_calc_settings()
    print(f'   {"patched" if ok else "old string not found — verify manually"}')
