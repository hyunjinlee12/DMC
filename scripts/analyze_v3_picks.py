"""Per-pick analysis of the 70 v3 DFT candidates: site, key distances, adsorbate geometry."""
import json, csv
from pathlib import Path
import numpy as np
from ase.io import read

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
SHORT = ROOT/'calculations/G3_adsorption/DFT_shortlist_v3'
G2 = ROOT/'calculations/G2_slab'
SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

CUTOFF_NBR = 2.6

def site_type(a, anchor, slab_idx):
    syms = a.get_chemical_symbols()
    d = a.get_distances(anchor, slab_idx, mic=True)
    nbrs = [(slab_idx[i], d[i]) for i in range(len(slab_idx)) if d[i] < CUTOFF_NBR]
    n_pd = sum(1 for i,_ in nbrs if syms[i]=='Pd')
    n_o  = sum(1 for i,_ in nbrs if syms[i]=='O')
    tot = n_pd + n_o
    if tot == 0: return 'physi', 0, 0
    if tot == 1: return ('atop_Pd' if n_pd else 'atop_O'), n_pd, n_o
    if tot == 2:
        lab = 'br_PdPd' if n_pd==2 else ('br_OO' if n_o==2 else 'br_PdO')
        return lab, n_pd, n_o
    if tot == 3:
        lab = 'h3Pd' if n_pd==3 else ('h3O' if n_o==3 else f'h3({n_pd}Pd{n_o}O)')
        return lab, n_pd, n_o
    return f'{tot}f({n_pd}Pd{n_o}O)', n_pd, n_o

def z_surf_top(a, slab_idx):
    return float(a.positions[slab_idx, 2].max())

def nearest_surface_dist(a, anchor, slab_idx):
    d = a.get_distances(anchor, slab_idx, mic=True)
    i_min = int(np.argmin(d))
    return float(d[i_min]), slab_idx[i_min]

def analyze():
    picks = json.load(open(SHORT/'summary.json'))
    rows = []
    for p in picks:
        sid, ads, rank, idx = p['sid'], p['ads'], p['rank'], p['idx']
        E = p['E']
        slab = read(G2/SDIRS[sid]/'CONTCAR'); n_sub = len(slab)
        f = next((SHORT/sid/ads).glob(f"{rank:02d}_*_idx{idx:05d}.vasp"))
        a = read(f)
        syms = a.get_chemical_symbols()
        c_i = [i for i,s in enumerate(syms) if s=='C']
        o_i = [i for i,s in enumerate(syms) if s=='O']
        h_i = [i for i,s in enumerate(syms) if s=='H']
        row = {'sid':sid,'ads':ads,'rank':rank,'idx':idx,'E_MLIP':round(E,4),
               'natoms':len(a),'n_sub':n_sub}

        # Identify ads O(s): O closest to each C (only ads O bonds to C directly)
        if ads == 'CO':
            c = c_i[0]
            o_ads = min(o_i, key=lambda oi: a.get_distance(c,oi,mic=True))
            ads_atoms = {c, o_ads}
        elif ads == 'CH3O':
            c = c_i[0]
            o_bonded = min(o_i, key=lambda oi: a.get_distance(c,oi,mic=True))
            ads_atoms = {c, o_bonded, *h_i}
        else:  # coads: 2 C + 3 H + 2 O (one O per C)
            me_c, co_c = None, None
            for c in c_i:
                nh = sum(1 for h in h_i if a.get_distance(c,h,mic=True) < 1.3)
                if nh == 3: me_c = c
                else: co_c = c
            me_o = min(o_i, key=lambda oi: a.get_distance(me_c,oi,mic=True))
            co_o = min([oi for oi in o_i if oi != me_o], key=lambda oi: a.get_distance(co_c,oi,mic=True))
            ads_atoms = {me_c, co_c, me_o, co_o, *h_i}
        slab_idx = [i for i in range(len(a)) if i not in ads_atoms]
        z_top = z_surf_top(a, slab_idx)

        if ads == 'CO':
            d_co = a.get_distance(c, o_ads, mic=True)
            site, npd, no = site_type(a, c, slab_idx)
            dmin, imin = nearest_surface_dist(a, c, slab_idx)
            row.update({'anchor':'C', 'site':site, 'n_Pd_nbr':npd, 'n_O_nbr':no,
                        'd_CO':round(d_co,3),
                        'd_anchor_surf':round(dmin,3),
                        'anchor_nbr_species':syms[imin],
                        'h_anchor':round(a.positions[c,2]-z_top,3)})
        elif ads == 'CH3O':
            d_co = a.get_distance(c, o_bonded, mic=True)
            site, npd, no = site_type(a, o_bonded, slab_idx)
            dmin, imin = nearest_surface_dist(a, o_bonded, slab_idx)
            row.update({'anchor':'O', 'site':site, 'n_Pd_nbr':npd, 'n_O_nbr':no,
                        'd_CO':round(d_co,3),
                        'd_anchor_surf':round(dmin,3),
                        'anchor_nbr_species':syms[imin],
                        'h_anchor':round(a.positions[o_bonded,2]-z_top,3)})
        else:
            d_meC_meO = a.get_distance(me_c, me_o, mic=True)
            d_coC_coO = a.get_distance(co_c, co_o, mic=True)
            d_react = a.get_distance(co_c, me_o, mic=True)
            site_co, npd_co, no_co = site_type(a, co_c, slab_idx)
            site_me, npd_me, no_me = site_type(a, me_o, slab_idx)
            dmin_co, _ = nearest_surface_dist(a, co_c, slab_idx)
            dmin_me, _ = nearest_surface_dist(a, me_o, slab_idx)
            row.update({'anchor':'C_CO+O_CH3',
                        'site_CO':site_co, 'site_CH3O':site_me,
                        'd_CO':round(d_coC_coO,3),
                        'd_CH3O_bond':round(d_meC_meO,3),
                        'd_reactive':round(d_react,3),
                        'd_C_CO_surf':round(dmin_co,3),
                        'd_O_CH3_surf':round(dmin_me,3),
                        'h_C_CO':round(a.positions[co_c,2]-z_top,3),
                        'h_O_CH3':round(a.positions[me_o,2]-z_top,3)})
        rows.append(row)
    return rows

if __name__ == '__main__':
    rows = analyze()
    out = ROOT/'calculations/G3_adsorption/DFT_shortlist_v3/picks_analysis.json'
    json.dump(rows, open(out,'w'), indent=2)
    print(f'Wrote {out}')
    # Also CSV summary
    csv_out = out.with_suffix('.csv')
    keys = sorted({k for r in rows for k in r})
    with open(csv_out,'w') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows: w.writerow(r)
    print(f'Wrote {csv_out}')
    # Print compact table
    print()
    for kind in ['CO','CH3O','coads']:
        print(f'\n=== {kind} ===')
        subset = [r for r in rows if r['ads']==kind]
        if kind == 'coads':
            hdr = f"{'sid':<5}{'rank':<5}{'idx':>6}{'E':>10}  {'site_CO':<10}{'site_CH3O':<12}{'d(C-O)':>8}{'d(C-O)ch3':>10}{'d_react':>9}"
        else:
            hdr = f"{'sid':<5}{'rank':<5}{'idx':>6}{'E':>10}  {'site':<10}{'nPd':>4}{'nO':>4}{'d(C-O)':>8}{'d_anc_surf':>12}{'h_anc':>8}"
        print(hdr); print('-'*len(hdr))
        for r in subset:
            if kind == 'coads':
                print(f"{r['sid']:<5}{r['rank']:<5}{r['idx']:>6}{r['E_MLIP']:>10.3f}  {r['site_CO']:<10}{r['site_CH3O']:<12}{r['d_CO']:>8.3f}{r['d_CH3O_bond']:>10.3f}{r['d_reactive']:>9.3f}")
            else:
                print(f"{r['sid']:<5}{r['rank']:<5}{r['idx']:>6}{r['E_MLIP']:>10.3f}  {r['site']:<10}{r['n_Pd_nbr']:>4}{r['n_O_nbr']:>4}{r['d_CO']:>8.3f}{r['d_anchor_surf']:>12.3f}{r['h_anchor']:>8.3f}")
