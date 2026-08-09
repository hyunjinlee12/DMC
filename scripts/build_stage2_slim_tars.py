"""Build 3 slim tars for Stage 2 batches per reviewer 2026-07-17 round 8.

  stage2_magnetic_seed.tar.gz       — 11 static seed dirs
  stage2_gas_vacuum_primary.tar.gz  — 4 gas vacuum dirs
  stage2_vaspsol_pilot.tar.gz       — 1 STATIC S1_clean_pilot dir

Each tar:
  - manifest.csv (per-dir tabulated)
  - README.md
  - sha256sums.txt (all included files)

Every candidate dir:
  - metadata.json extended with { purpose, stage, analysis_use }

Not included:
  - mask_sensitivity_test  (wait for magnetic seed result)
  - gas *_vaspsol         (wait for VASPsol pilot)
  - gas *_static           (wait for gas vacuum relax)
  - T1.17 non-pilot dirs   (wait for Stage 3 decisions)
  - setup_t1_16_final_vacuum.py / setup_t1_17_vaspsol_final.py
    (write only after Stage 2 decisions land)

Repo dirs are NOT deleted — only the 3 slim tars are distinct artifacts.
"""
import hashlib, json, subprocess, shutil, csv, os, re
from pathlib import Path

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
STAGE = Path('/tmp/stage2_slim_tars')
if STAGE.exists(): shutil.rmtree(STAGE)
STAGE.mkdir()

# --------- helpers ---------
def sha256_of(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for chunk in iter(lambda: fh.read(65536), b''): h.update(chunk)
    return h.hexdigest()

def write_sha256_manifest(root_dir):
    lines = []
    for p in sorted(root_dir.rglob('*')):
        if not p.is_file(): continue
        if p.name == 'sha256sums.txt': continue
        lines.append(f'{sha256_of(p)}  {p.relative_to(root_dir)}')
    (root_dir/'sha256sums.txt').write_text('\n'.join(lines) + '\n')

def update_metadata(dest_dir, extra):
    m = dest_dir/'metadata.json'
    if m.exists():
        meta = json.loads(m.read_text())
    else:
        meta = {}
    meta.update(extra)
    m.write_text(json.dumps(meta, indent=2))

def make_tar(work_dir, tar_path):
    """work_dir contains one top-level named dir; tar it up."""
    subprocess.run(['tar', '-czf', str(tar_path), '-C', str(work_dir.parent), work_dir.name],
                   check=True)

# =====================================================================
# TAR 1: stage2_magnetic_seed
# =====================================================================
print('\n== stage2_magnetic_seed ==')
seed_dir = STAGE/'stage2_magnetic_seed'
seed_dir.mkdir()
src_seed = ROOT/'calculations/magnetic_seed_test'
shutil.copytree(src_seed, seed_dir/'magnetic_seed_test')

# Update each candidate dir's metadata with purpose/stage/analysis_use
seed_metadata_common = {
    'stage': 'Stage 2 — magnetic seed test (electronic state only, no ionic relaxation)',
    'analysis_use': ('Pilot only. Do NOT use these energies as final. Purpose: '
                     'compare final total/local magnetization + E across 3 initial '
                     'MAGMOM seeds per system to identify true magnetic ground '
                     'state. If a seed gives lower E by > 0.02 eV, that state is '
                     'adopted; if within 0.02 eV, run a NUPDOWN scan. Results '
                     'drive the MAGMOM policy for downstream T1_16_final_vacuum '
                     'and T1_17_VASPsol.'),
}
for d in sorted((seed_dir/'magnetic_seed_test').iterdir()):
    if not d.is_dir(): continue
    name = d.name
    if 'S1_clean' in name and 'ispin1' in name:
        purpose = 'Reference nonmagnetic (ISPIN=1). Not same as ISPIN=2, MAGMOM=0.'
    elif 'S1_clean' in name and 'lowmag' in name:
        purpose = 'Low-moment seed (Pd 0.2 μB); tests weakly-magnetic solution.'
    elif 'S1_clean' in name and 'highmag' in name:
        purpose = 'High-moment seed (Pd 1.0 μB); reproduces current T1.16 default.'
    elif 'S1_CO_idx3' in name and 'ispin1' in name:
        purpose = 'S1/CO nonmagnetic seed (ISPIN=1). Reference for even-electron system.'
    elif 'S1_CO_idx3' in name and 'lowmag' in name:
        purpose = 'S1/CO low-moment seed (Pd 0.2 μB, C/O = 0).'
    elif 'S1_CO_idx3' in name and 'highmag' in name:
        purpose = 'S1/CO high-moment seed (Pd 1.0 μB, C/O = 0); T1.16 default reproduction.'
    elif 'S3_clean' in name and 'nonmag' in name:
        purpose = 'S3 clean nonmagnetic (ISPIN=1); PdO low-spin baseline.'
    elif 'S3_clean' in name and 'fm' in name:
        purpose = 'S3 clean FM seed (Pd 0.5 μB).'
    elif 'S3_clean' in name and 'afm' in name:
        purpose = 'S3 clean AFM seed (Pd ±0.5, exactly 32/32 by frac-x).'
    elif 'CH3O' in name and 'doublet' in name:
        purpose = 'S3/CH3O radical doublet (NUPDOWN=1, MAGMOM on ads O).'
    elif 'CH3O' in name and 'unconstrained' in name:
        purpose = 'S3/CH3O radical unconstrained (no NUPDOWN, Pd 0.3 + O 1).'
    else:
        purpose = 'magnetic seed test'
    update_metadata(d, {**seed_metadata_common, 'purpose': purpose})

# Manifest
manifest_rows = []
for d in sorted((seed_dir/'magnetic_seed_test').iterdir()):
    if not d.is_dir(): continue
    meta = json.loads((d/'metadata.json').read_text()) if (d/'metadata.json').exists() else {}
    incar = (d/'INCAR').read_text()
    ispin = next((int(l.split('=')[1].strip()) for l in incar.splitlines()
                  if l.strip().startswith('ISPIN')), None)
    nupdown = next((int(l.split('=')[1].strip()) for l in incar.splitlines()
                    if l.strip().startswith('NUPDOWN')), None)
    magmom = next((l.split('=',1)[1].strip() for l in incar.splitlines()
                   if l.strip().startswith('MAGMOM')), '')
    manifest_rows.append({'dir':d.name, 'ISPIN':ispin, 'NUPDOWN':nupdown,
                          'MAGMOM': (magmom[:80] + '...') if len(magmom) > 80 else magmom,
                          'source': meta.get('source',''), 'purpose': meta.get('purpose','')})

# Verify S3 AFM balance in manifest
afm_incar = (seed_dir/'magnetic_seed_test/S3_clean__afm/INCAR').read_text()
afm_line = next((l for l in afm_incar.splitlines() if l.strip().startswith('MAGMOM')), '')
p = m = 0
for tok in afm_line.split('=',1)[1].strip().split():
    if '*' in tok:
        n, v = tok.split('*'); v = float(v); n = int(n)
    else:
        n, v = 1, float(tok)
    if v > 0.01: p += n
    elif v < -0.01: m += n
assert p == 32 and m == 32, f'S3 AFM not balanced: +={p}, -={m}'
print(f'  S3 AFM balance verified: +0.5 × {p}, -0.5 × {m}')

with open(seed_dir/'manifest.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['dir','ISPIN','NUPDOWN','MAGMOM','source','purpose'])
    w.writeheader()
    for r in manifest_rows: w.writerow(r)

# Analysis helper script
(seed_dir/'analyze_seed_results.py').write_text('''"""After all 11 seed dirs finish, run this to summarize:
- final energy(sigma->0)
- final total magnetization
- per-orbital tot magnetization (s, p, d)
- convergence status

Usage:
    python analyze_seed_results.py > summary.txt
"""
import re
from pathlib import Path

for d in sorted(Path('magnetic_seed_test').iterdir()):
    if not d.is_dir(): continue
    outcar = d/'OUTCAR'
    if not outcar.exists():
        print(f'{d.name:<38} STATUS = not yet run')
        continue
    text = outcar.read_text()
    E = re.findall(r'energy\\(sigma->0\\)\\s*=\\s*(-?\\d+\\.\\d+)', text)
    mag_line = re.findall(r'number of electron.*?magnetization\\s+(-?\\d+\\.\\d+)', text)
    conv = 'reached required accuracy' in text
    E_final = float(E[-1]) if E else None
    m_final = float(mag_line[-1]) if mag_line else None
    print(f'{d.name:<38} E={E_final} μ_tot={m_final} conv={conv}')
''')

# README
(seed_dir/'README.md').write_text('''# Stage 2 — magnetic_seed_test

11 static single-point calcs (IBRION=-1, NSW=0) resolving the MAGMOM
initialization ambiguity for the T1.16/T1.17 pipeline.

## Systems

Even-electron (3 seeds each):
- S1_clean          (Pd metal, 80 atoms)
- S1_CO_idx3        (S1 top-1 CO adsorbate)
- S3_clean          (PdO 128 atoms; 3 seeds = nonmag / FM / AFM)

Odd-electron (2 seeds each):
- S3_CH3O_idx315    (radical; doublet NUPDOWN=1 / unconstrained)

Total: 3 + 3 + 3 + 2 = 11 dirs.

## Common INCAR

    IBRION = -1
    NSW = 0
    ISTART = 0
    ICHARG = 2
    ENCUT = 520
    PREC = Accurate
    IVDW = 12
    LDIPOL = .TRUE.  (slab dirs)
    KSPACING = 0.25
    LREAL = Auto

INCAR is IDENTICAL to T1.17 modulo the spin block (ISPIN/NUPDOWN/MAGMOM per seed)
and LSOL (off — vacuum). No WAVECAR/CHGCAR shared → seed independence.

## S3 AFM balance (verified)

MAGMOM alternates ± 0.5 μB across Pd atoms by fractional-x ordering.
Exactly 32 Pd with +0.5 and 32 Pd with -0.5. Net moment = 0. Substrate O atoms = 0.

## Post-run analysis

    python analyze_seed_results.py > summary.txt

Then read `summary.txt` and apply the decision rule:

- Seed ΔE > 0.02 eV, one clearly lower → adopt that state. Run one ionic
  relaxation on the representative structure to verify the magnetic
  solution survives geometry optimization.
- Seed ΔE < 0.02 eV, different moments → NUPDOWN fixed-spin scan.
- All seeds converge to same E + moment → magnetic solution stable.

S2/S3b/S4 magnetic ground states are NOT covered by these 11 dirs. After
T1.16 completes, inspect each surface's final total magnetization from
OUTCAR; add spot-check seeds for surfaces where the value looks anomalous.

## H200 submit

    sbatch submit_vasp.sh   (in each dir)

Each dir has partition=h200q, 16 CPU, Singularity+srun. VASPsol NOT
required (LSOL not set).
''')

write_sha256_manifest(seed_dir)
make_tar(seed_dir, STAGE/'stage2_magnetic_seed.tar.gz')

# =====================================================================
# TAR 2: stage2_gas_vacuum_primary
# =====================================================================
print('\n== stage2_gas_vacuum_primary ==')
gas_dir = STAGE/'stage2_gas_vacuum_primary'
gas_dir.mkdir()

for mol in ['CO', 'CH3O', 'CH3OH', 'H2']:
    src = ROOT/f'calculations/gas_references/{mol}_vacuum'
    dst = gas_dir/f'{mol}_vacuum'
    shutil.copytree(src, dst)
    purpose_map = {
        'CO':    ('Isolated CO in cubic box (edge 15 Å). Provides μ_CO(g) for '
                  'ΔG_ads formulae, applicable regardless of slab environment '
                  '(vacuum or VASPsol) as long as CO is supplied as gas feed.'),
        'CH3O':  ('Isolated CH3O radical (open-shell doublet, NUPDOWN=1). '
                  'Provides G(CH3O_radical, gas) reference for the radical '
                  'convention of ΔG_CH3O*. Not needed if only the CHE MeOH(U) '
                  'convention is used.'),
        'CH3OH': ('Isolated CH3OH (closed-shell). Gas-phase electronic energy '
                  'reference. NOTE: not equal to μ(CH3OH liquid, activity=1); '
                  'liquid standard state requires RT ln(p_sat/p°) correction '
                  '(≈ -4.40 kJ/mol at 298 K) or AIMD.'),
        'H2':    ('Isolated H2 for the CHE reference ½ G(H2). Standard SHE '
                  'convention: μ(H+) + e- = ½ μ(H2, g, 1 bar) at U=0 SHE.'),
    }
    update_metadata(dst, {
        'stage': 'Stage 2 — gas vacuum primary reference',
        'purpose': purpose_map[mol],
        'analysis_use': ('CORE reference for T1.18 ΔG computation. Provides '
                         'E_DFT; ZPE + vibrational free energy + trans/rot '
                         'entropy + pressure/activity corrections must be added '
                         'via ASE Thermochemistry to obtain full G.'),
    })

manifest = []
for d in sorted(gas_dir.iterdir()):
    if not d.is_dir(): continue
    meta = json.loads((d/'metadata.json').read_text())
    manifest.append({'mol':d.name, 'box_A':meta.get('box_A'),
                     'ISPIN':meta.get('ISPIN'), 'NUPDOWN':meta.get('NUPDOWN'),
                     'MAGMOM':meta.get('MAGMOM'), 'purpose':meta.get('purpose','')[:80]+'...'})
with open(gas_dir/'manifest.csv','w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['mol','box_A','ISPIN','NUPDOWN','MAGMOM','purpose'])
    w.writeheader()
    for r in manifest: w.writerow(r)

(gas_dir/'README.md').write_text('''# Stage 2 — gas_vacuum_primary

4 isolated-molecule vacuum relaxations. These are the **core reference pool**
for T1.18 adsorption-energy computations. Not all four appear in every ΔG
formula — the ones used depend on the reservoir convention chosen for each
species.

## Included

| dir | ISPIN | NUPDOWN | MAGMOM | role |
|---|---|---|---|---|
| CO_vacuum    | 1 | – | –             | μ_CO(g) — for any convention where CO is gas feed |
| CH3O_vacuum  | 2 | 1 | 0 0 0 0 1 (POSCAR order C H H H O) | radical reference (only if radical convention chosen) |
| CH3OH_vacuum | 1 | – | –             | CHE MeOH(U) numerator |
| H2_vacuum    | 1 | – | –             | ½ G(H2) for CHE |

## reservoir vs environment (important)

The molecular reservoir chemical potential μ_species is determined by the
**physical supply** — gas feed, dissolved, or liquid — and is INDEPENDENT of
whether the slab pair is calculated in vacuum or VASPsol. Do NOT force the
molecular reference to match the slab environment.

Example:
    slab pair in VASPsol + CO supplied as gas feed
    → μ_CO from CO_vacuum, not from a solvated CO calculation

VASPsol calculations of the same molecules (`CO_vaspsol`, `CH3O_vaspsol`, ...)
are needed only if a **dissolved-species** convention is chosen (rare in this
project), or if a `Δμ_solvation` correction is being computed.

## From E_DFT to G

The energies you extract from these 4 OUTCARs are electronic `E_DFT` only.
Full Gibbs free energy G needs ASE Thermochemistry post-processing:

    G(T, p) = E_DFT + ZPE + ∫ Cv dT − T S(vib, rot, trans) + kT ln(p/p°)

For a liquid reservoir (e.g. methanol as pure liquid at activity 1):

    μ(liq, T) = μ(g, T, p°) + RT ln(p_sat(T)/p°)
    methanol @ 298.15 K: p_sat ≈ 0.169 bar → correction ≈ -4.40 kJ/mol

**Do not add +ΔG_vap** — the correction is NEGATIVE at 298 K because
methanol's saturation vapor pressure is below 1 bar.

## Cubic box

Each molecule is centered in a **cubic box of edge 15 Å** (volume ≈ 3375 Å³).
Not "15 Å³" — that would be an atomic-scale volume. Explicit Γ-only KPOINTS.
Polar-molecule finite-cell error will be checked with paired-static 20 Å
runs later (NOT included in this tar).

## H200 submit

    sbatch submit_vasp.sh   (in each dir)

Each dir: partition=h200q, 16 CPU, Singularity+srun. Runtime per dir ≈ 5–30 min.
''')

write_sha256_manifest(gas_dir)
make_tar(gas_dir, STAGE/'stage2_gas_vacuum_primary.tar.gz')

# =====================================================================
# TAR 3: stage2_vaspsol_pilot
# =====================================================================
print('\n== stage2_vaspsol_pilot ==')
pilot_dir = STAGE/'stage2_vaspsol_pilot'
pilot_dir.mkdir()
src_pilot = ROOT/'calculations/T1_17_VASPsol/S1_clean'
dst_pilot = pilot_dir/'S1_clean_pilot'
shutil.copytree(src_pilot, dst_pilot)

# Convert to STATIC (IBRION=-1, NSW=0) per reviewer
incar_p = dst_pilot/'INCAR'
lines = incar_p.read_text().splitlines()
new_lines = []
for ln in lines:
    if ln.strip().startswith('IBRION'): new_lines.append('IBRION = -1')
    elif ln.strip().startswith('NSW'): new_lines.append('NSW = 0')
    elif ln.strip().startswith('EDIFFG'): pass  # not meaningful for static
    else: new_lines.append(ln)
incar_p.write_text('\n'.join(new_lines) + '\n')

# PILOT_ONLY marker
(dst_pilot/'PILOT_ONLY').write_text(
    "This directory is a VASPsol functional pilot.\n"
    "IBRION = -1, NSW = 0 (STATIC single-point).\n"
    "Purpose: verify LSOL/EB_K/TAU are recognized, solvation output present,\n"
    "no unknown-INCAR-tag warnings. The energy from this calc is NOT used\n"
    "in any final analysis.\n"
)

update_metadata(dst_pilot, {
    'stage': 'Stage 2 — VASPsol functional pilot',
    'purpose': ('Verify VASPsol build support (LSOL/EB_K/TAU recognized). '
                'Static single-point in methanol implicit solvent on the '
                'S1 clean slab. No ionic relaxation.'),
    'analysis_use': ('PILOT ONLY — do not use this energy in any final analysis. '
                     'Success criteria: OUTCAR contains solvation output '
                     '(cavity, dielectric response), no unknown INCAR tag '
                     'warnings, electronic convergence reached. If pilot passes, '
                     'proceed to submit molecular VASPsol calcs; if it fails, '
                     'rebuild VASP with solvation source.'),
    'ionic_relaxation': False,
})

(pilot_dir/'README.md').write_text('''# Stage 2 — VASPsol functional pilot

**1 dir, STATIC single-point.** Purpose: verify the H200 VASP build actually
implements VASPsol (i.e. that LSOL/EB_K/TAU are recognized and solvation is
applied) before submitting the 8 molecular VASPsol calcs later.

## Submit

    cd S1_clean_pilot/
    sbatch submit_vasp.sh

Expected walltime: minutes (single SCF).

## Success criteria

Grep OUTCAR after completion:

    grep -E "LSOL|EB_K|TAU|VASPsol|cavity|solvation" OUTCAR

Must show:
- LSOL flag read
- EB_K = 32.6 accepted
- TAU = 0 accepted
- Solvation output (dielectric cavity, contribution to energy)
- No "unknown INCAR tag" warnings for these flags
- Electronic SCF converged

## What if it fails

If unknown-INCAR-tag warnings appear or no solvation output, the VASP build
does not include VASPsol. Rebuild VASP with VASPsol source patch. Do NOT
submit any `*_vaspsol` molecule dir until this pilot passes.

## What NOT to do

- Do not use the energy from this pilot in any ΔG or descriptor map.
  This is a functional test only.
- Do not extend NSW>0 to make it a "real" T1.17 calc — the T1.17 top-1
  candidates will be chosen in Stage 3 after magnetic seed + mask
  sensitivity results.
''')

write_sha256_manifest(pilot_dir)
make_tar(pilot_dir, STAGE/'stage2_vaspsol_pilot.tar.gz')

print('\n=== all 3 slim tars built ===')
for t in ['stage2_magnetic_seed.tar.gz','stage2_gas_vacuum_primary.tar.gz','stage2_vaspsol_pilot.tar.gz']:
    p = STAGE/t
    print(f'  {p}  {p.stat().st_size:>10} bytes  ({sha256_of(p)[:16]}...)')
