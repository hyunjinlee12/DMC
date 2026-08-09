"""Consolidate raw data into paper-ready CSV tables.

Read-only aggregation of what already exists in the repo:
- G1 bulk OUTCARs → lattice + E_bulk table
- G2 slab OUTCAR + CONTCAR → slab reference table (rumpling, bond lengths)
- MACE-D3 references JSON → flat CSV
- MLIP unique JSON (Phase 1 singles + Phase 2 filtered coads) → row-per-candidate
- v5add16 combined manifest → 86-row DFT shortlist
- T1_16_DFT_L2 OUTCAR results (completed jobs only) → DFT results table

Output: paper_data/{01..07}_*.csv + README.md
"""
import json, csv, re
from pathlib import Path
from collections import Counter
import numpy as np
from ase.io import read

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
OUT = ROOT/'paper_data'
OUT.mkdir(exist_ok=True)

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

# ---------- helpers ----------
def parse_outcar_energy(p):
    """Return dict with final E_sigma_zero, E_free (F), max_force (free atoms only)."""
    text = Path(p).read_text()
    # 'energy(sigma->0)' — last occurrence = converged step
    e_sigma = None; e_free = None
    for m in re.finditer(r'energy\s*\(sigma->0\)\s*=\s*(-?\d+\.\d+)', text):
        e_sigma = float(m.group(1))
    for m in re.finditer(r'free\s+energy\s+TOTEN\s*=\s*(-?\d+\.\d+)', text):
        e_free = float(m.group(1))
    reached = 'reached required accuracy' in text
    # max force on non-fixed atoms
    fmax = None
    return {'E_sigma0': e_sigma, 'E_free': e_free, 'reached_accuracy': reached}

def parse_outcar_forces_max(outcar, contcar):
    """Get max |force| on free (not selective-dynamics-fixed) atoms."""
    atoms = read(contcar)
    text = Path(outcar).read_text()
    # take last force block: skip header + separator, then read numeric rows until next '-----'
    blocks = list(re.finditer(r'TOTAL-FORCE.*?\n\s*-----+\n([\s\S]+?)\n\s*-----', text))
    if not blocks: return None
    lines = blocks[-1].group(1).strip().split('\n')
    forces = []
    for ln in lines:
        parts = ln.split()
        if len(parts) >= 6:
            forces.append([float(x) for x in parts[3:6]])
    if not forces: return None
    forces = np.array(forces)
    if forces.ndim != 2 or forces.shape[1] != 3: return None
    if len(forces) != len(atoms): return float(np.linalg.norm(forces, axis=1).max())
    # apply constraint mask
    mask = np.ones(len(atoms), bool)
    if atoms.constraints:
        from ase.constraints import FixAtoms
        for c in atoms.constraints:
            if isinstance(c, FixAtoms):
                mask[c.index] = False
    return float(np.linalg.norm(forces[mask], axis=1).max()) if mask.any() else None

def lattice_from_contcar(p):
    a = read(p)
    cell = a.cell
    lengths = a.cell.lengths()
    angles = a.cell.angles()
    vol = a.get_volume()
    return {
        'a_A': round(float(lengths[0]),4),
        'b_A': round(float(lengths[1]),4),
        'c_A': round(float(lengths[2]),4),
        'alpha_deg': round(float(angles[0]),2),
        'beta_deg':  round(float(angles[1]),2),
        'gamma_deg': round(float(angles[2]),2),
        'volume_A3': round(float(vol),3),
        'natoms': len(a),
    }

# ---------- 01: bulk data ----------
EXP_LAT = {
    'Pd':   {'a':3.891, 'note':'fcc, expt (Owen 1927, refined many times)'},
    'PdO':  {'a':3.043, 'c':5.336, 'note':'tetragonal P4₂/mmc'},
    'PdO2': {'a':4.437, 'c':3.089, 'note':'rutile P4₂/mnm (calc. reference)'},
}
KMESH = {'Pd':'12x12x12','PdO':'8x8x6','PdO2':'6x6x8'}

def build_bulk():
    rows = []
    for name in ['Pd','PdO','PdO2','PdO2_PBE']:
        d = ROOT/'calculations/G1_bulk'/name
        if not (d/'OUTCAR').exists(): continue
        lat = lattice_from_contcar(d/'CONTCAR')
        e   = parse_outcar_energy(d/'OUTCAR')
        exp = EXP_LAT.get(name if name != 'PdO2_PBE' else 'PdO2', {})
        row = {
            'material': name,
            'functional': 'PBE-D3(BJ)' if name != 'PdO2_PBE' else 'PBE (no D3)',
            'ENCUT_eV': 520,
            'k_mesh': KMESH.get(name if name != 'PdO2_PBE' else 'PdO2',''),
            **lat,
            'E_bulk_eV_sigma0': e['E_sigma0'],
            'E_per_atom_eV': round(e['E_sigma0']/lat['natoms'], 4) if e['E_sigma0'] else None,
            'reached_accuracy': e['reached_accuracy'],
            'exp_a_A': exp.get('a',''),
            'exp_c_A': exp.get('c',''),
            'dev_a_percent': round(100*(lat['a_A']-exp.get('a',lat['a_A']))/exp['a'],3) if 'a' in exp else '',
            'dev_c_percent': round(100*(lat['c_A']-exp.get('c',lat['c_A']))/exp['c'],3) if 'c' in exp else '',
            'notes': exp.get('note',''),
        }
        rows.append(row)
    keys = list(rows[0].keys())
    with open(OUT/'01_bulk_data.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    print(f'01_bulk_data.csv: {len(rows)} rows')

# ---------- 02: slab data ----------
SLAB_MODEL = {
    'S1':  'Pd(100), p(4×4), 5 layers, bottom 2 fixed',
    'S2':  '1 ML PdO(101)/Pd(100), (√5×√5)R27°, 4 layers Pd + oxide top',
    'S3':  'PdO(100) O-rich, mixed termination, 128 atoms',
    'S3b': 'PdO(100) Pd-terminated, 104 atoms',
    'S4':  'PdO₂(110) stoichiometric, 144 atoms',
}
SLAB_OXIDATION = {'S1':'Pd⁰','S2':'Pd²⁺ interface','S3':'Pd²⁺ oxide',
                  'S3b':'Pd²⁺ oxide (Pd-term)','S4':'Pd⁴⁺'}

def rumpling_and_Pd_O(contcar):
    a = read(contcar)
    syms = a.get_chemical_symbols()
    z = a.positions[:,2]
    top_thr = z.max() - 1.5
    top_idx = [i for i,zi in enumerate(z) if zi>=top_thr]
    top_syms = Counter(syms[i] for i in top_idx)
    rumpling = float(z[top_idx].max() - z[top_idx].min()) if top_idx else 0.0
    # nearest Pd-O bond in the top region
    pd_idx = [i for i in top_idx if syms[i]=='Pd']
    o_idx  = [i for i,s in enumerate(syms) if s=='O']
    if pd_idx and o_idx:
        d = a.get_distances(pd_idx[0], o_idx, mic=True)
        d_min = float(d.min())
        # avg of Pd-O within 2.6 A across all top Pd
        pairs=[]
        for p in pd_idx:
            dd = a.get_distances(p, o_idx, mic=True)
            for x in dd:
                if x < 2.6: pairs.append(x)
        d_avg = float(np.mean(pairs)) if pairs else None
    else:
        d_min = None; d_avg = None
    return {
        'natoms': len(a),
        'n_Pd': syms.count('Pd'),
        'n_O':  syms.count('O'),
        'top_layer': ', '.join(f'{el}:{n}' for el,n in top_syms.most_common()),
        'rumpling_top_A': round(rumpling, 3),
        'd_PdO_min_top_A': round(d_min,3) if d_min else '',
        'd_PdO_avg_top_A': round(d_avg,3) if d_avg else '',
    }

def build_slab():
    rows = []
    for sid, sdir in SDIRS.items():
        d = ROOT/'calculations/G2_slab'/sdir
        outcar = d/'OUTCAR'; contcar = d/'CONTCAR'
        if not outcar.exists() or not contcar.exists(): continue
        e = parse_outcar_energy(outcar)
        fmax = parse_outcar_forces_max(outcar, contcar)
        geo = rumpling_and_Pd_O(contcar)
        row = {
            'sid': sid,
            'model': SLAB_MODEL.get(sid,''),
            'oxidation_state': SLAB_OXIDATION.get(sid,''),
            **geo,
            'E_slab_eV_sigma0': round(e['E_sigma0'],4) if e['E_sigma0'] else None,
            'F_max_free_eV_per_A': round(fmax,4) if fmax else None,
            'reached_accuracy': e['reached_accuracy'],
            'vacuum_A': 20,
            'ISMEAR': 1 if sid=='S1' else 0,
            'SIGMA_eV': 0.10 if sid=='S1' else 0.05,
            'LDIPOL': True,
            'IDIPOL': 3,
        }
        rows.append(row)
    keys = list(rows[0].keys())
    with open(OUT/'02_slab_data.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    print(f'02_slab_data.csv: {len(rows)} rows')

# ---------- 03: MACE-D3 references ----------
def build_references():
    refs = json.load(open(ROOT/'calculations/G3_adsorption/mace_d3_references.json'))
    rows = []
    for sid, E in refs['slab'].items():
        rows.append({'kind':'slab','name':sid,'E_MACE_D3_eV':round(E,6),
                     'note':'MACE-MH+D3(BJ)+cueq, MACE-relaxed slab, bottom-half fixed'})
    for gname, E in refs['gas'].items():
        note = {'CO':'MACE-relaxed CO in vacuum',
                'CH3OH':'MACE-relaxed methanol',
                'H2':'MACE-relaxed H2',
                'CH3O_ref':'derived: E(CH3OH) − ½ E(H2), MeOH(U)-style reference',
                'CH3O_radical':'MACE-relaxed CH3O radical (spin-polarised)'}.get(gname,'')
        rows.append({'kind':'gas','name':gname,'E_MACE_D3_eV':round(E,6),'note':note})
    with open(OUT/'03_mace_references.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['kind','name','E_MACE_D3_eV','note']); w.writeheader()
        for r in rows: w.writerow(r)
    print(f'03_mace_references.csv: {len(rows)} rows')

# ---------- 04: MLIP pool singles ----------
def E_bind_single(sid, ads, E, refs):
    ref = refs['slab'][sid] + (refs['gas']['CO'] if ads=='CO' else refs['gas']['CH3O_radical'])
    return E - ref

def build_mlip_singles():
    refs = json.load(open(ROOT/'calculations/G3_adsorption/mace_d3_references.json'))
    rows = []
    for sid, sdir in SDIRS.items():
        for ads in ['CO','CH3O']:
            p = ROOT/'calculations/G3_adsorption'/sdir/'MLIP_phase1'/f'unique_{ads}.json'
            if not p.exists(): continue
            data = json.load(open(p))
            for r in data:
                E = float(r['E'])
                rows.append({
                    'sid': sid, 'ads': ads, 'idx': r['idx'],
                    'E_MACE_D3_eV': round(E, 6),
                    'E_bind_MACE_D3_eV': round(E_bind_single(sid, ads, E, refs), 4),
                    'dE_rel_meV': round(float(r.get('dE_rel', 0))*1000, 1),
                    'converged': bool(r.get('converged', True)),
                    'n_steps': int(r.get('n_steps', -1)),
                    'd_min_ads_sub_A': round(float(r.get('d_min', 0)), 3),
                    'site_type_raw_MLIP': r.get('site_type',''),
                    'fingerprint': json.dumps(r.get('fingerprint',[])),
                })
    keys = list(rows[0].keys())
    with open(OUT/'04_mlip_pool_singles.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    print(f'04_mlip_pool_singles.csv: {len(rows)} rows')

# ---------- 05: MLIP pool coads ----------
def dist_bin(d):
    if d<2.1: return 'product-like'
    if d<3.0: return 'reactive-close'
    if d<4.0: return 'reactive-loose'
    if d<5.0: return 'separated'
    return 'thermodynamic'

def build_mlip_coads():
    refs = json.load(open(ROOT/'calculations/G3_adsorption/mace_d3_references.json'))
    rows = []
    for sid in ['S1','S2','S3','S3b']:   # S4 coads excluded per project decision
        sdir = SDIRS[sid]
        p = ROOT/'calculations/G3_adsorption'/sdir/'MLIP_phase2_filtered'/'unique_SetA.json'
        if not p.exists(): continue
        data = json.load(open(p))
        for r in data:
            E = float(r['E'])
            E_ref = refs['slab'][sid] + refs['gas']['CO'] + refs['gas']['CH3O_radical']
            d_react = float(r.get('d_reactive', 0))
            rows.append({
                'sid': sid, 'idx': r['idx'],
                'E_MACE_D3_eV': round(E, 6),
                'E_bind_MACE_D3_eV': round(E - E_ref, 4),
                'dE_rel_meV': round(float(r.get('dE_rel_meV', 0)), 1),
                'converged': bool(r.get('converged', True)),
                'n_steps': int(r.get('n_steps', -1)),
                'd_min_ads_sub_A': round(float(r.get('d_min', 0)), 3),
                'd_reactive_A': round(d_react, 3),
                'distance_bin': dist_bin(d_react),
                'fingerprint': json.dumps(r.get('fingerprint',[])),
            })
    keys = list(rows[0].keys())
    with open(OUT/'05_mlip_pool_coads.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    print(f'05_mlip_pool_coads.csv: {len(rows)} rows')

# ---------- 06: DFT shortlist (86 rows, reuse combined_summary + selection_candidate_rationale) ----------
def build_shortlist():
    src = ROOT/'calculations/T1_16_DFT_L2_v5add16/combined_summary.csv'
    rat = ROOT/'calculations/G3_adsorption/DFT_shortlist_v5/audit/selection_candidate_rationale.csv'
    rows_a = {(r['sid'],r['ads'],r['idx']):r for r in csv.DictReader(open(src))}
    rows_b = {(r['sid'],r['ads'],r['idx']):r for r in csv.DictReader(open(rat))}
    out = []
    for k, ra in rows_a.items():
        rb = rows_b.get(k, {})
        out.append({**ra, 'dft_hypothesis': rb.get('dft_hypothesis','')})
    keys = list(out[0].keys())
    with open(OUT/'06_dft_shortlist.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in out: w.writerow(r)
    print(f'06_dft_shortlist.csv: {len(out)} rows')

# ---------- 07: DFT results (parse OUTCARs where available) ----------
def build_dft_results():
    rows = []
    for r in csv.DictReader(open(OUT/'06_dft_shortlist.csv')):
        dir_rel = r['dir']
        outcar = ROOT/dir_rel/'OUTCAR'
        contcar = ROOT/dir_rel/'CONTCAR'
        result = {
            'sid': r['sid'], 'ads': r['ads'], 'idx': r['idx'],
            'origin': r['origin'], 'priority': r['priority'],
            'E_MLIP_eV': r['E_MLIP'], 'E_bind_MLIP_eV': r['E_bind_MLIP'],
            'DFT_status': 'PENDING',
            'E_DFT_sigma0_eV': '', 'E_bind_DFT_eV': '',
            'F_max_free_eV_per_A': '', 'reached_accuracy': '',
            'contcar_path': str(contcar.relative_to(ROOT)) if contcar.exists() else '',
        }
        if outcar.exists():
            e = parse_outcar_energy(outcar)
            if e['E_sigma0'] is not None:
                fmax = parse_outcar_forces_max(outcar, contcar) if contcar.exists() else None
                result['DFT_status'] = 'DONE' if e['reached_accuracy'] else 'RUNNING_OR_UNCONVERGED'
                result['E_DFT_sigma0_eV'] = round(e['E_sigma0'], 4)
                result['F_max_free_eV_per_A'] = round(fmax, 4) if fmax else ''
                result['reached_accuracy'] = e['reached_accuracy']
        rows.append(result)
    keys = list(rows[0].keys())
    with open(OUT/'07_dft_results.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    n_done = sum(1 for r in rows if r['DFT_status']=='DONE')
    print(f'07_dft_results.csv: {len(rows)} rows ({n_done} DONE, {len(rows)-n_done} PENDING)')

# ---------- README.md ----------
def build_readme():
    md = ['# Paper-ready raw data\n\n',
          f'Consolidated from repo state; regenerable via `scripts/build_paper_data.py`.\n',
          f'All energies in eV. All distances in Å. All angles in degrees.\n\n',
          '## Files\n\n',
          '- **01_bulk_data.csv** — G1 bulk relaxation results for Pd, PdO, PdO₂ '
              '(and one PdO₂ w/o D3 for comparison). Includes lattice a/b/c, angles, '
              'volume, total energy, per-atom energy, experimental lattice references, '
              'and deviation from experiment.\n',
          '- **02_slab_data.csv** — G2 clean-slab results for the 5 modeled surfaces '
              '(S1 Pd(100), S2 PdO(101)/Pd(100), S3 PdO(100) O-rich, S3b PdO(100) '
              'Pd-terminated, S4 PdO₂(110)). Includes atom counts, top-layer '
              'composition, rumpling, min/avg Pd–O bond length in the top region, '
              'total energy, max force on free atoms.\n',
          '- **03_mace_references.csv** — reference energies used to compute '
              'E_bind_MACE_D3 in downstream tables. Same MACE-MH+D3(BJ)+cueq calculator '
              'as the pool relaxations.\n',
          '- **04_mlip_pool_singles.csv** — every MLIP-relaxed unique CO*/CH₃O* '
              'candidate across the 5 surfaces (~1200 rows). Columns: sid, ads, idx, '
              'E_MACE_D3, E_bind_MACE_D3, dE_rel_meV, converged, n_steps, '
              'd_min_ads_sub, site_type (raw MLIP label), fingerprint.\n',
          '- **05_mlip_pool_coads.csv** — every MLIP-relaxed unique coadsorption '
              'candidate (Set A) for S1/S2/S3/S3b (~14600 rows). Adds d_reactive '
              '(C_CO ↔ O_CH3O, Å, MIC) and distance_bin classification.\n',
          '- **06_dft_shortlist.csv** — 86 candidates selected for DFT L2 (70 from v4 '
              'selector + 16 from v5 additions). Includes site/region/distance-bin '
              'labels, priority (v4-baseline / MUST / RECOMMENDED / OPTIONAL / '
              'MUST-diagnostic), selection reason, DFT hypothesis.\n',
          '- **07_dft_results.csv** — 86 rows in one-to-one correspondence with 06; '
              'DFT_status is one of DONE / RUNNING_OR_UNCONVERGED / PENDING. Once a '
              'given candidate finishes, `scripts/build_paper_data.py` re-run '
              'automatically fills E_DFT_sigma0, F_max_free, reached_accuracy.\n\n',
          '## Provenance / methods (short form)\n\n',
          '- Functional: PBE + D3(BJ), IVDW=12, ENCUT=520 eV, PREC=Accurate, LASPH, ADDGRID.\n',
          '- Bulk: ISIF=3 (full relaxation), ISMEAR=1 for Pd (σ=0.10), ISMEAR=0 for PdO/PdO₂ (σ=0.05).\n',
          '- Slab: ISIF=2 (ionic only, cell fixed from bulk), bottom-half atoms FixAtoms, '
              'vacuum 20 Å, LDIPOL=True, IDIPOL=3.\n',
          '- MLIP: MACE-MH mh-1 head oc20_usemppbe, float64, cueq enabled, '
              'D3(BJ, xc=pbe) dispersion active. LBFGS fmax=0.05 eV/Å, 200–300 steps.\n',
          '- DFT L2: same slab INCAR as L1 (vacuum); EDIFFG=-0.03 eV/Å; ISPIN=2.\n',
          '- Convention: **E_bind is defined against isolated gas monomers** '
              '(CO gas + CH₃O_radical gas). More negative = stronger binding.\n\n',
          '## Sign conventions & unit conventions\n\n',
          '| quantity | unit | sign |\n|---|---|---|\n',
          '| E_bulk, E_slab, E_MLIP, E_DFT | eV | absolute (arbitrary offset, only diffs meaningful) |\n',
          '| E_bind | eV | more negative = stronger binding |\n',
          '| ΔE (relative to global min) | eV | ≥0 |\n',
          '| lattice a/b/c | Å | positive |\n',
          '| force | eV/Å | positive magnitude |\n',
          '| d_reactive, d_PdO, d_anchor_surf | Å | positive |\n',
          '| rumpling | Å | positive |\n\n',
          '## Regenerate\n\n',
          '```bash\n',
          'python scripts/build_paper_data.py\n',
          '```\n',
          'The script is idempotent and re-parses OUTCARs each time so completed '
              'DFT jobs are automatically picked up in **07_dft_results.csv**.\n',
    ]
    (OUT/'README.md').write_text(''.join(md))
    print('README.md written')

if __name__=='__main__':
    build_bulk()
    build_slab()
    build_references()
    build_mlip_singles()
    build_mlip_coads()
    build_shortlist()
    build_dft_results()
    build_readme()
    print('done.')
