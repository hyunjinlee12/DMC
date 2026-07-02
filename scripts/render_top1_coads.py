"""Top-1 lowest-E co-ads per surface, bc-plane side view, same style as CO*/CH3O*."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ase.io import read
from render_structure import render

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G3 = ROOT/'calculations/G3_adsorption'
G2 = ROOT/'calculations/G2_slab'
OUT = ROOT/'reports/predft_advisor_figures/top1_coads_sideviews'
OUT.mkdir(parents=True, exist_ok=True)

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}


def intramol_valid_coads(atoms):
    syms = atoms.get_chemical_symbols()
    c_idx = [i for i,s in enumerate(syms) if s=='C']
    h_idx = [i for i,s in enumerate(syms) if s=='H']
    o_all = [i for i,s in enumerate(syms) if s=='O']
    if len(c_idx) != 2 or len(h_idx) != 3: return False
    # For CO: find C with 1 close O
    # For methoxy: find C with 3 H + 1 O
    for c in c_idx:
        o_dist = sorted([(atoms.get_distance(c,oi,mic=True), oi) for oi in o_all])
        d_h = sorted([atoms.get_distance(c,hi,mic=True) for hi in h_idx])
        n_close_h = sum(1 for d in d_h if d<1.3)
        if n_close_h == 0:  # CO carbon
            if not (1.05 <= o_dist[0][0] <= 1.30): return False
            if len(o_dist)>=2 and o_dist[1][0]<1.5: return False  # exclude CO2 (2nd O too close)
        elif n_close_h == 3:  # methoxy carbon
            if not all(0.90<=d<=1.25 for d in d_h): return False
            if not (1.30 <= o_dist[0][0] <= 1.55): return False
            if len(o_dist)>=2 and o_dist[1][0]<1.5: return False
        else:
            return False
    return True


for sid, sdir in SDIRS.items():
    slab = read(G2/sdir/'CONTCAR'); n_sub = len(slab)
    unique = json.load(open(G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'))
    traj = list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
    # helper: classify anchor to substrate site type by neighbor count (2.6 A cutoff)
    def site_type(atoms, anchor):
        syms = atoms.get_chemical_symbols()
        n_sub_atoms = len(atoms) - 7  # slab atoms
        sub_indices = list(range(n_sub_atoms))
        d = atoms.get_distances(anchor, sub_indices, mic=True)
        nbrs = [(sub_indices[i], d[i]) for i in range(len(sub_indices)) if d[i] < 2.6]
        n_pd = sum(1 for i,_ in nbrs if syms[i]=='Pd')
        n_o = sum(1 for i,_ in nbrs if syms[i]=='O')
        total = n_pd + n_o
        if total == 0: return 'physi'
        if total == 1: return f'atop_{"Pd" if n_pd else "O"}'
        if total == 2:
            if n_pd == 2: return 'br_PdPd'
            if n_o == 2: return 'br_OO'
            return 'br_PdO'
        if total == 3:
            if n_pd == 3: return 'hollow_3Pd'
            if n_o == 3: return 'hollow_3O'
            return f'h3({n_pd}Pd{n_o}O)'
        return f'{total}f'

    # Sort by E, pick top-3 with DIFFERENT (CO_site, OMe_site) combo
    picks = []
    combos_seen = set()
    for r in sorted(unique, key=lambda r: r['E']):
        a = traj[r['idx']]
        if len(a) != n_sub + 7:
            ads = a[-7:]; a = slab.copy(); a += ads
        if not intramol_valid_coads(a): continue
        syms = a.get_chemical_symbols()
        c_idx = [i for i,s in enumerate(syms) if s=='C']
        h_idx = [i for i,s in enumerate(syms) if s=='H']
        me_c, co_c = None, None
        for c in c_idx:
            n_h = sum(1 for h in h_idx if a.get_distance(c,h,mic=True) < 1.3)
            if n_h == 3: me_c = c
            else: co_c = c
        if me_c is None or co_c is None: continue
        o_all = [i for i,s in enumerate(syms) if s=='O']
        me_o_cand = sorted([(a.get_distance(me_c,oi,mic=True), oi) for oi in o_all])
        me_o = me_o_cand[0][1]
        combo = (site_type(a, co_c), site_type(a, me_o))
        if combo in combos_seen: continue
        combos_seen.add(combo)
        picks.append((r, a, combo))
        if len(picks) >= 3: break

    for rank, (r, a, combo) in enumerate(picks, start=1):
        print(f'{sid} top-{rank}: idx={r["idx"]}, E={r["E"]:.3f}, sites=CO[{combo[0]}]/OMe[{combo[1]}]')
        a2 = a * (2,2,1)
        render(a2, OUT/f'{sid}_coads_top{rank}.png', rotation='-90x,-90y,0z',
               width=2400, show_cell=False)
