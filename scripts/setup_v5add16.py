"""Build T1_16_DFT_L2_v5add16/ — 16 new v5-approved candidates.

Same INCAR/POTCAR/submit convention as T1_16_DFT_L2 (setup_dft_v3 with the
H-element POTCAR bug fixed). No files under T1_16_DFT_L2/ are touched.
No jobs submitted.

Deliverables:
  T1_16_DFT_L2_v5add16/{sid}/{ads}/{rank:02d}_{ads}_idx{idx:05d}/
      POSCAR, INCAR, POTCAR, submit_vasp_gpu.sh, metadata.json
  T1_16_DFT_L2_v5add16/manifest.csv    all 70 v4 + 16 v5-new
  T1_16_DFT_L2_v5add16/verification.md  POSCAR/POTCAR species check + bash -n
"""
import json, shutil, subprocess, csv
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT/'calculations/G2_slab'
G3 = ROOT/'calculations/G3_adsorption'
L2 = ROOT/'calculations/T1_16_DFT_L2'
OUT = ROOT/'calculations/T1_16_DFT_L2_v5add16'
if OUT.exists(): shutil.rmtree(OUT)
OUT.mkdir()

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

# ---------- 16 new picks per researcher approval (2026-07-17) ----------
NEW_PICKS = [
    # coads (8)
    dict(sid='S1',  ads='coads', idx=487,  priority='RECOMMENDED',
         reason='reactive-close br_PdPd+br_PdPd (ΔE=0.032)'),
    dict(sid='S1',  ads='coads', idx=4021, priority='RECOMMENDED',
         reason='thermodynamic br_PdPd+br_PdPd rep (ΔE=0.011)'),
    dict(sid='S2',  ads='coads', idx=3633, priority='RECOMMENDED',
         reason='thermodynamic br_PdPd+atop_Pd rep (ΔE=0.229)'),
    dict(sid='S2',  ads='coads', idx=8079, priority='RECOMMENDED',
         reason='reactive-close atop_Pd+atop_Pd distinct-site (ΔE=0.441)'),
    dict(sid='S3',  ads='coads', idx=3481, priority='MUST',
         reason='reactive-close atop_Pd+atop_Pd, bin gap not covered by v4 (ΔE=0.041)'),
    dict(sid='S3',  ads='coads', idx=5161, priority='MUST',
         reason='reactive-loose atop_Pd+atop_Pd, bin gap not covered by v4 (ΔE=0.021)'),
    dict(sid='S3b', ads='coads', idx=2051, priority='MUST-diagnostic',
         reason='reactive-close br_PdPd+atop_Pd, only reactive-close in pool. '
                'HIGH-ENERGY DIAGNOSTIC ΔE=0.554 eV — kept per workplan §9 rule 9'),
    dict(sid='S3b', ads='coads', idx=2754, priority='MUST',
         reason='reactive-loose br_PdPd+br_PdPd (ΔE=0.122)'),
    # single-ads OPTIONAL (8, per user's explicit list)
    dict(sid='S1',  ads='CO',   idx=43,  priority='OPTIONAL',
         reason='4-fold hollow site (4f_4Pd0O) diversity (ΔE=0.046)'),
    dict(sid='S1',  ads='CH3O', idx=283, priority='OPTIONAL',
         reason='atop_Pd/metal site (ΔE=0.002)'),
    dict(sid='S2',  ads='CH3O', idx=217, priority='OPTIONAL',
         reason='br_PdO interface site (ΔE=0.288)'),
    dict(sid='S2',  ads='CH3O', idx=496, priority='OPTIONAL',
         reason='atop_O oxide site (ΔE=0.322)'),
    dict(sid='S3',  ads='CH3O', idx=395, priority='OPTIONAL',
         reason='physisorbed reference (ΔE=0.029)'),
    dict(sid='S3b', ads='CO',   idx=6,   priority='OPTIONAL',
         reason='br_PdPd metal site (ΔE=0.134)'),
    dict(sid='S4',  ads='CH3O', idx=329, priority='OPTIONAL',
         reason='br_PdO interface site (ΔE=0.042)'),
    dict(sid='S4',  ads='CH3O', idx=383, priority='OPTIONAL',
         reason='physisorbed reference (ΔE=0.006)'),
]

# ---------- INCAR templates (identical to setup_dft_v3.py) ----------
INCAR_METAL = """SYSTEM = pddmc
ENCUT = 520
PREC = Accurate
LASPH = .TRUE.
ADDGRID = .TRUE.
ISPIN = 2
IVDW = 12
EDIFF = 1e-06
NELM = 500
NELMIN = 5
ALGO = Normal
NCORE = 1
LREAL = Auto
LWAVE = .FALSE.
LCHARG = .FALSE.
LORBIT = 11
ISYM = 0
IBRION = 2
NSW = 300
ISIF = 2
ISMEAR = 1
SIGMA = 0.1
EDIFFG = -0.03
LDIPOL = .TRUE.
IDIPOL = 3
KSPACING = 0.25
"""
INCAR_OXIDE = INCAR_METAL.replace("ISMEAR = 1\nSIGMA = 0.1", "ISMEAR = 0\nSIGMA = 0.05")

# ---------- POTCAR library ----------
POT_LIB = Path('/home/hyunjin/POTENTIAL/potpaw_PBE')
POT_FOLDER = {'C':'C', 'H':'H', 'O':'O', 'Pd':'Pd_pv'}

# ---------- validity + descriptor (reused) ----------
CO2_CUT = 2.0
def valid_CO(a):
    syms=a.get_chemical_symbols(); c=[i for i,s in enumerate(syms) if s=='C'][0]
    o_i=[i for i,s in enumerate(syms) if s=='O']
    d_o=sorted([(a.get_distance(c,oi,mic=True),oi) for oi in o_i])
    return (c, d_o[0][1])
def valid_CH3O(a):
    syms=a.get_chemical_symbols(); c=[i for i,s in enumerate(syms) if s=='C'][0]
    o_i=[i for i,s in enumerate(syms) if s=='O']
    d_o=sorted([(a.get_distance(c,oi,mic=True),oi) for oi in o_i])
    return (c, d_o[0][1])
def valid_coads(a):
    syms=a.get_chemical_symbols()
    c_i=[i for i,s in enumerate(syms) if s=='C']
    h_i=[i for i,s in enumerate(syms) if s=='H']
    o_i=[i for i,s in enumerate(syms) if s=='O']
    me_c,co_c=None,None
    for c in c_i:
        nh=sum(1 for h in h_i if a.get_distance(c,h,mic=True)<1.3)
        if nh==3: me_c=c
        else: co_c=c
    d_co=sorted([(a.get_distance(co_c,oi,mic=True),oi) for oi in o_i])
    co_o=d_co[0][1]
    d_me=sorted([(a.get_distance(me_c,oi,mic=True),oi) for oi in o_i if oi!=co_o])
    return (co_c, co_o, me_c, d_me[0][1])

def fix_bottom_half(atoms, n_sub):
    z = atoms.positions[:n_sub, 2]
    zm = np.median(z)
    fixed=[i for i in range(n_sub) if atoms.positions[i,2] < zm]
    atoms.set_constraint(FixAtoms(indices=fixed))

# ---------- reference submit script ----------
SUBMIT_REF = L2/'S1/CO/00_CO_idx00064/submit_vasp_gpu.sh'
SUBMIT_TXT = SUBMIT_REF.read_text()

# ---------- v5 descriptor CSV to fetch per-pick attributes ----------
v5_desc = list(csv.DictReader(open(G3/'DFT_shortlist_v5/proposed_additions_only.csv')))
desc_map = {(r['sid'], r['ads'], int(r['idx'])): r for r in v5_desc}

# ---------- Build ----------
built = []
issues = []
for p in NEW_PICKS:
    sid, ads, idx = p['sid'], p['ads'], p['idx']
    sdir = SDIRS[sid]
    slab = read(G2/sdir/'CONTCAR'); n_sub = len(slab)
    if ads == 'coads':
        uj = json.load(open(G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'))
        traj = list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
        n_ads = 7; vf = valid_coads
    else:
        uj = json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
        traj = list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
        n_ads = 2 if ads=='CO' else 5
        vf = valid_CO if ads=='CO' else valid_CH3O
    ur = next(r for r in uj if r['idx'] == idx)
    E = float(ur['E'])
    a = traj[idx]
    if len(a) != n_sub + n_ads:
        ads_at = a[-n_ads:]
        a = slab.copy()
        a += ads_at
    vf(a)  # runs but return not needed here
    fix_bottom_half(a, n_sub)

    dest = OUT/sid/ads/f'{ads}_idx{idx:05d}'
    dest.mkdir(parents=True)
    # POSCAR (sort=True groups species — matches setup_dft_v3)
    write(str(dest/'POSCAR'), a, format='vasp', direct=True, sort=True, vasp5=True)
    # INCAR
    incar = INCAR_METAL if sid == 'S1' else INCAR_OXIDE
    (dest/'INCAR').write_text(incar)
    # POTCAR built from library, species order = POSCAR 6th line
    species = (dest/'POSCAR').read_text().splitlines()[5].split()
    potcar = ''.join((POT_LIB/POT_FOLDER[s]/'POTCAR').read_text() for s in species)
    (dest/'POTCAR').write_text(potcar)
    # submit script (identical to L2)
    (dest/'submit_vasp_gpu.sh').write_text(SUBMIT_TXT)
    # metadata
    desc = desc_map.get((sid, ads, idx), {})
    meta = {
        'sid': sid, 'ads': ads, 'idx': idx,
        'E_MLIP': E,
        'dE_from_global': float(desc.get('dE_from_global','') or 0),
        'site': desc.get('site',''),
        'region': desc.get('region',''),
        'n_Pd_nbr': desc.get('n_Pd_nbr',''),
        'n_O_nbr': desc.get('n_O_nbr',''),
        'd_anchor_surf': desc.get('d_anchor_surf',''),
        'd_reactive': desc.get('d_reactive',''),
        'd_bin': desc.get('d_bin',''),
        'site_combo': desc.get('site_combo',''),
        'priority': p['priority'],
        'selection_reason': p['reason'],
        'neighbors_within_3.5A': desc.get('neighbors_within_3.5A',''),
        'source_traj': str(G3/sdir/(
            'MLIP_phase2/relaxed_SetA.traj' if ads=='coads'
            else f'MLIP_phase1/relaxed_{ads}.traj')),
        'source_unique_json': str(G3/sdir/(
            'MLIP_phase2_filtered/unique_SetA.json' if ads=='coads'
            else f'MLIP_phase1/unique_{ads}.json')),
        'notes': ('HIGH-ENERGY DIAGNOSTIC: only reactive-close candidate in pool; '
                  'ΔE=0.554 eV vs global min. DFT will show whether it truly binds '
                  'or relaxes back to a lower-E arrangement.'
                  if (sid,ads,idx)==('S3b','coads',2051) else ''),
    }
    json.dump(meta, open(dest/'metadata.json','w'), indent=2)

    # ---- verification: POSCAR species vs POTCAR TITEL order
    # POTCAR TITEL line: "TITEL  = PAW_PBE <species> <date>" → species is field [3]
    potcar_titles = [ln.split()[3].split('_')[0]
                     for ln in potcar.split('\n') if ln.strip().startswith('TITEL')]
    poscar_species = species
    if potcar_titles != poscar_species:
        issues.append(f'{dest}: POSCAR species {poscar_species} != POTCAR TITEL {potcar_titles}')
    # ---- verification: bash -n
    r = subprocess.run(['bash','-n', str(dest/'submit_vasp_gpu.sh')],
                       capture_output=True, text=True)
    if r.returncode != 0:
        issues.append(f'{dest}/submit_vasp_gpu.sh bash -n failed: {r.stderr[:120]}')

    built.append({'sid':sid,'ads':ads,'idx':idx,'dir':str(dest.relative_to(ROOT)),
                  'E_MLIP':E,'priority':p['priority'],
                  'poscar_species':' '.join(poscar_species),
                  'potcar_species':' '.join(potcar_titles),
                  'bash_n':'OK' if r.returncode==0 else 'FAIL'})

# ---------- manifest.csv: 70 v4 (existing) + 16 v5-new ----------
manifest_rows = []

# 70 v4 rows
v4 = json.load(open(G3/'DFT_shortlist_v3/summary.json'))
v4_desc = list(csv.DictReader(open(G3/'DFT_shortlist_v3/picks_analysis.csv')))
v4_desc_map = {(r['sid'], r['ads'], int(r['idx'])): r for r in v4_desc}
for r in v4:
    sid, ads, idx, rank = r['sid'], r['ads'], r['idx'], r['rank']
    d = v4_desc_map.get((sid, ads, idx), {})
    dir_rel = f"calculations/T1_16_DFT_L2/{sid}/{ads}/{rank:02d}_{ads}_idx{idx:05d}"
    # status
    if sid == 'S1' and ads == 'CO':
        status = 'DONE'
    else:
        status = 'PENDING'
    manifest_rows.append({
        'origin':'v4-existing',
        'sid':sid, 'ads':ads, 'idx':idx, 'rank':rank,
        'E_MLIP':r['E'], 'dE_from_global':d.get('dE_from_global',''),
        'site': (d.get('site','') if ads!='coads' else ''),
        'site_combo': (f"{d.get('site_CO','')}+{d.get('site_CH3O','')}"
                        if ads=='coads' and d.get('site_CO') else ''),
        'region': d.get('region',''),
        'd_reactive': d.get('d_reactive',''),
        'd_bin':'',   # v4 file has no d_bin
        'priority':'v4-baseline',
        'selection_reason':'top-5 MLIP E + xy dedup (v4 selector)',
        'status':status,
        'dir':dir_rel,
    })

# 16 v5-new rows
for r in built:
    d = desc_map.get((r['sid'], r['ads'], r['idx']), {})
    reason = next(p['reason'] for p in NEW_PICKS
                  if p['sid']==r['sid'] and p['ads']==r['ads'] and p['idx']==r['idx'])
    manifest_rows.append({
        'origin':'v5-new',
        'sid':r['sid'], 'ads':r['ads'], 'idx':r['idx'], 'rank':'',
        'E_MLIP':r['E_MLIP'], 'dE_from_global':d.get('dE_from_global',''),
        'site':d.get('site',''),
        'site_combo':d.get('site_combo',''),
        'region':d.get('region',''),
        'd_reactive':d.get('d_reactive',''),
        'd_bin':d.get('d_bin',''),
        'priority':r['priority'],
        'selection_reason':reason,
        'status':'PENDING (bundle only, not submitted)',
        'dir':r['dir'],
    })

keys=['origin','sid','ads','idx','rank','E_MLIP','dE_from_global','site',
      'site_combo','region','d_reactive','d_bin','priority','selection_reason',
      'status','dir']
with open(OUT/'manifest.csv','w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
    for r in manifest_rows: w.writerow(r)

# ---------- verification.md ----------
md=[]
md.append('# v5add16 bundle — verification report\n\n')
md.append(f'Built {len(built)} new candidate dirs under `{OUT.relative_to(ROOT)}/`\n\n')
md.append('## Per-candidate verification\n\n')
md.append('| # | sid | ads | idx | priority | POSCAR species | POTCAR TITEL | bash -n |\n')
md.append('|---|---|---|---|---|---|---|---|\n')
for i,b in enumerate(built,1):
    ok_species = '✅' if b['poscar_species']==b['potcar_species'] else '❌'
    md.append(f"| {i} | {b['sid']} | {b['ads']} | {b['idx']} | {b['priority']} | "
              f"`{b['poscar_species']}` | `{b['potcar_species']}` {ok_species} | {b['bash_n']} |\n")
md.append(f'\n## Issues found: {len(issues)}\n')
if issues:
    for e in issues: md.append(f'- {e}\n')
else:
    md.append('- None. All POSCAR/POTCAR species orders match; all submit scripts pass `bash -n`.\n')
md.append('\n## Guardrails observed\n\n')
md.append(f'- `T1_16_DFT_L2/` untouched (source of INCAR/POTCAR templates, read-only).\n')
md.append(f'- No jobs submitted. Each dir contains submit_vasp_gpu.sh identical to v4 baseline;\n')
md.append(f'  H200 environment paths must be adjusted before use, same as pending65 bundle.\n')
md.append(f'- `manifest.csv` distinguishes 70 v4-existing (5 DONE, 65 PENDING) vs 16 v5-new (all PENDING).\n')
md.append(f'\n## Special note\n\n')
md.append(f'- **S3b/coads idx=2051** flagged `MUST-diagnostic` in metadata (ΔE=0.554 eV).\n')
md.append(f'  Only reactive-close candidate available in the S3b coads MLIP pool. DFT will\n')
md.append(f'  determine whether it remains bound at this configuration or relaxes to a\n')
md.append(f'  lower-E arrangement.\n')
(OUT/'verification.md').write_text(''.join(md))

# ---------- report to stdout ----------
print(f'Built {len(built)} candidates under {OUT.relative_to(ROOT)}/')
print(f'manifest.csv: {len(manifest_rows)} rows (70 v4 + {len(built)} v5-new)')
print(f'verification issues: {len(issues)}')
if issues:
    for e in issues: print(f'  {e}')
else:
    print('  ✅ all POSCAR/POTCAR species match, all bash -n pass')
