"""Read-only audit: compute sensitivity of site labels + group summary tables.

Does NOT change selection, does NOT create dirs, does NOT submit jobs.
Writes only to calculations/G3_adsorption/DFT_shortlist_v5/audit/.
"""
import json, csv, math
from pathlib import Path
from collections import Counter
import numpy as np
from ase.io import read

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
G2 = ROOT/'calculations/G2_slab'
G3 = ROOT/'calculations/G3_adsorption'
V5 = G3/'DFT_shortlist_v5'
V5A = ROOT/'calculations/T1_16_DFT_L2_v5add16'
OUT = V5/'audit'
OUT.mkdir(exist_ok=True)

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}
CO2_CUT = 2.0

def valid_CO(a):
    syms=a.get_chemical_symbols(); c_i=[i for i,s in enumerate(syms) if s=='C']
    o_i=[i for i,s in enumerate(syms) if s=='O']
    if len(c_i)!=1: return None
    c=c_i[0]; d_o=sorted([(a.get_distance(c,oi,mic=True),oi) for oi in o_i])
    if not(1.05<=d_o[0][0]<=1.30): return None
    if len(d_o)>=2 and d_o[1][0]<CO2_CUT: return None
    return (c, d_o[0][1])

def valid_CH3O(a):
    syms=a.get_chemical_symbols()
    c_i=[i for i,s in enumerate(syms) if s=='C']
    h_i=[i for i,s in enumerate(syms) if s=='H']
    o_i=[i for i,s in enumerate(syms) if s=='O']
    if len(c_i)!=1 or len(h_i)!=3: return None
    c=c_i[0]
    if not all(0.90<=a.get_distance(c,h,mic=True)<=1.25 for h in h_i): return None
    d_o=sorted([(a.get_distance(c,oi,mic=True),oi) for oi in o_i])
    if not(1.30<=d_o[0][0]<=1.55): return None
    return (c, d_o[0][1])

def valid_coads(a):
    syms=a.get_chemical_symbols()
    c_i=[i for i,s in enumerate(syms) if s=='C']
    h_i=[i for i,s in enumerate(syms) if s=='H']
    o_i=[i for i,s in enumerate(syms) if s=='O']
    if len(c_i)!=2 or len(h_i)!=3: return None
    me_c,co_c=None,None
    for c in c_i:
        nh=sum(1 for h in h_i if a.get_distance(c,h,mic=True)<1.3)
        if nh==3: me_c=c
        else: co_c=c
    if me_c is None or co_c is None: return None
    d_co=sorted([(a.get_distance(co_c,oi,mic=True),oi) for oi in o_i])
    co_o=d_co[0][1]
    d_me=sorted([(a.get_distance(me_c,oi,mic=True),oi) for oi in o_i if oi!=co_o])
    if not d_me: return None
    return (co_c, co_o, me_c, d_me[0][1])

def site_label(a, anchor, slab_idx, cutoff):
    syms=a.get_chemical_symbols()
    d=a.get_distances(anchor, slab_idx, mic=True)
    nbrs=[(slab_idx[i],d[i]) for i in range(len(slab_idx)) if d[i]<cutoff]
    n_pd=sum(1 for i,_ in nbrs if syms[i]=='Pd')
    n_o=sum(1 for i,_ in nbrs if syms[i]=='O')
    tot=n_pd+n_o
    if tot==0: return 'physi'
    if tot==1: return 'atop_Pd' if n_pd else 'atop_O'
    if tot==2: return 'br_PdPd' if n_pd==2 else ('br_OO' if n_o==2 else 'br_PdO')
    if tot==3: return 'h3Pd' if n_pd==3 else ('h3O' if n_o==3 else f'h3_{n_pd}Pd{n_o}O')
    return f'{tot}f_{n_pd}Pd{n_o}O'

def dist_bin(d, lo=3.0, hi=4.0):
    if d<2.1: return 'product-like'
    if d<lo:  return 'reactive-close'
    if d<hi:  return 'reactive-loose'
    if d<5.0: return 'separated'
    return 'thermodynamic'

# ---------- load MLIP references ----------
refs=json.load(open(G3/'mace_d3_references.json'))
E_SLAB=refs['slab']; E_CO=refs['gas']['CO']; E_CH3O=refs['gas']['CH3O_radical']
def E_bind(sid, ads, E):
    ref=E_CO+E_CH3O if ads=='coads' else (E_CO if ads=='CO' else E_CH3O)
    return E - E_SLAB[sid] - ref

# ---------- pool loader ----------
def load_pool(sid, ads):
    sdir=SDIRS[sid]
    if ads=='coads':
        uj=json.load(open(G3/sdir/'MLIP_phase2_filtered/unique_SetA.json'))
        traj=list(read(G3/sdir/'MLIP_phase2/relaxed_SetA.traj', index=':'))
        vf=valid_coads; n_ads=7
    else:
        uj=json.load(open(G3/sdir/f'MLIP_phase1/unique_{ads}.json'))
        traj=list(read(G3/sdir/f'MLIP_phase1/relaxed_{ads}.traj', index=':'))
        vf=valid_CO if ads=='CO' else valid_CH3O
        n_ads=2 if ads=='CO' else 5
    slab_ref=read(G2/sdir/'CONTCAR'); n_ref=len(slab_ref)
    return uj, traj, vf, n_ads, slab_ref, n_ref

# ---------- sensitivity for 16 v5-new ----------
new_picks=[]
for r in csv.DictReader(open(V5A/'manifest.csv')):
    if r['origin']=='v5-new':
        new_picks.append(r)

sensitivity_rows=[]
for p in new_picks:
    sid=p['sid']; ads=p['ads']; idx=int(p['idx']); priority=p['priority']
    uj, traj, vf, n_ads, slab_ref, n_ref = load_pool(sid, ads)
    ur=next(x for x in uj if x['idx']==idx)
    E=float(ur['E'])
    dEg=float(p['dE_from_global']) if p['dE_from_global'] else 0.0
    a=traj[idx]
    if len(a)!=n_ref+n_ads:
        ads_at=a[-n_ads:]; a=slab_ref.copy(); a+=ads_at
    v=vf(a)
    if ads=='coads':
        co_c, co_o, me_c, me_o = v
        h_i=[i for i,s in enumerate(a.get_chemical_symbols()) if s=='H']
        ads_set={co_c,co_o,me_c,me_o,*h_i}
        slab_idx=[i for i in range(len(a)) if i not in ads_set]
        d_react=float(a.get_distance(co_c, me_o, mic=True))
        labels={}
        for cut in [2.4, 2.6, 2.8]:
            s_co=site_label(a, co_c, slab_idx, cut)
            s_me=site_label(a, me_o, slab_idx, cut)
            labels[f'combo@{cut}A']=f'{s_co}+{s_me}'
        bins={}
        for lo,hi in [(2.8,3.8),(3.0,4.0),(3.2,4.2)]:
            bins[f'bin@{lo}-{hi}A']=dist_bin(d_react, lo, hi)
        row={'sid':sid,'ads':ads,'idx':idx,'priority':priority,
             'E_MLIP':round(E,4),'E_bind_MLIP':round(E_bind(sid,ads,E),4),
             'dE_from_global':round(dEg,4),
             'd_reactive':round(d_react,3),
             **labels, **bins,
             'label_stable': 'yes' if len(set(labels.values()))==1 else 'no',
             'bin_stable': 'yes' if len(set(bins.values()))==1 else 'no'}
    else:
        c, o_bonded = v
        anchor = c if ads=='CO' else o_bonded
        if ads=='CH3O':
            h_i=[i for i,s in enumerate(a.get_chemical_symbols()) if s=='H']
            ads_set={c, o_bonded, *h_i}
        else:
            ads_set={c, o_bonded}
        slab_idx=[i for i in range(len(a)) if i not in ads_set]
        labels={f'site@{cut}A': site_label(a, anchor, slab_idx, cut) for cut in [2.4, 2.6, 2.8]}
        row={'sid':sid,'ads':ads,'idx':idx,'priority':priority,
             'E_MLIP':round(E,4),'E_bind_MLIP':round(E_bind(sid,ads,E),4),
             'dE_from_global':round(dEg,4),
             'd_reactive':'',
             **labels,
             'bin@2.8-3.8A':'', 'bin@3.0-4.0A':'', 'bin@3.2-4.2A':'',
             'label_stable': 'yes' if len(set(labels.values()))==1 else 'no',
             'bin_stable': 'N/A'}
    # Classify robustness
    if row['label_stable']=='yes' and row.get('bin_stable') in ('yes','N/A'):
        row['robustness']='robust'
    elif row['label_stable']=='no' and (row.get('bin_stable') in ('yes','N/A') or ads=='coads'):
        row['robustness']='threshold-sensitive-label'
    elif row.get('bin_stable')=='no':
        row['robustness']='threshold-sensitive-bin'
    else:
        row['robustness']='mixed'
    sensitivity_rows.append(row)

keys=['sid','ads','idx','priority','E_MLIP','E_bind_MLIP','dE_from_global','d_reactive',
      'site@2.4A','site@2.6A','site@2.8A',
      'combo@2.4A','combo@2.6A','combo@2.8A',
      'bin@2.8-3.8A','bin@3.0-4.0A','bin@3.2-4.2A',
      'label_stable','bin_stable','robustness']
for r in sensitivity_rows:
    for k in keys:
        if k not in r: r[k]=''
with open(OUT/'selection_sensitivity.csv','w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
    for r in sensitivity_rows: w.writerow(r)
print(f'sensitivity: {len(sensitivity_rows)} rows')

# ---------- group summary ----------
pool_diag=list(csv.DictReader(open(V5/'pool_diagnostics.csv')))
diag_map={(r['sid'],r['ads']):r for r in pool_diag}
manifest=list(csv.DictReader(open(V5A/'manifest.csv')))

group_rows=[]
for sid in ['S1','S2','S3','S3b','S4']:
    for ads in ['CO','CH3O','coads']:
        if ads=='coads' and sid=='S4': continue
        uj, traj, vf, n_ads, slab_ref, n_ref = load_pool(sid, ads)
        E_vals=[]
        site_counter=Counter(); bin_counter=Counter()
        for r in uj:
            if not r.get('converged',True): continue
            a=traj[r['idx']]
            if len(a)!=n_ref+n_ads:
                ads_at=a[-n_ads:]; a=slab_ref.copy(); a+=ads_at
            v=vf(a)
            if v is None: continue
            E_vals.append(r['E'])
            if ads=='coads':
                co_c, co_o, me_c, me_o = v
                h_i=[i for i,s in enumerate(a.get_chemical_symbols()) if s=='H']
                ads_set={co_c,co_o,me_c,me_o,*h_i}
                slab_idx=[i for i in range(len(a)) if i not in ads_set]
                s_co=site_label(a, co_c, slab_idx, 2.6)
                s_me=site_label(a, me_o, slab_idx, 2.6)
                site_counter[f'{s_co}+{s_me}']+=1
                d_react=a.get_distance(co_c, me_o, mic=True)
                bin_counter[dist_bin(d_react)]+=1
            else:
                c, o = v
                anchor = c if ads=='CO' else o
                if ads=='CH3O':
                    h_i=[i for i,s in enumerate(a.get_chemical_symbols()) if s=='H']
                    ads_set={c,o,*h_i}
                else:
                    ads_set={c,o}
                slab_idx=[i for i in range(len(a)) if i not in ads_set]
                site_counter[site_label(a, anchor, slab_idx, 2.6)]+=1
        d = diag_map.get((sid,ads), {})
        n_v4=sum(1 for r in manifest if r['origin']=='v4-existing' and r['sid']==sid and r['ads']==ads)
        n_v5=sum(1 for r in manifest if r['origin']=='v5-new' and r['sid']==sid and r['ads']==ads)
        group_rows.append({
            'sid':sid, 'ads':ads,
            'raw':d.get('raw',''),
            'unique_MLIP':d.get('unique_json',''),
            'converged':d.get('converged',''),
            'valid_final':len(E_vals),
            'n_sites_distinct':len(site_counter),
            'top_sites':';'.join(f'{k}:{v}' for k,v in site_counter.most_common(5)),
            'n_distance_bins':len(bin_counter) if ads=='coads' else '',
            'distance_bins':(';'.join(f'{k}:{v}' for k,v in bin_counter.most_common())
                             if ads=='coads' else ''),
            'E_min':round(min(E_vals),4) if E_vals else '',
            'E_max':round(max(E_vals),4) if E_vals else '',
            'dE_span':round(max(E_vals)-min(E_vals),3) if E_vals else '',
            'E_bind_min':round(E_bind(sid,ads,min(E_vals)),4) if E_vals else '',
            'E_bind_max':round(E_bind(sid,ads,max(E_vals)),4) if E_vals else '',
            'n_v4_picks':n_v4,
            'n_v5_new_picks':n_v5,
            'total_dft_planned':n_v4+n_v5,
        })

keys=['sid','ads','raw','unique_MLIP','converged','valid_final',
      'n_sites_distinct','top_sites','n_distance_bins','distance_bins',
      'E_min','E_max','dE_span','E_bind_min','E_bind_max',
      'n_v4_picks','n_v5_new_picks','total_dft_planned']
with open(OUT/'selection_group_summary.csv','w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
    for r in group_rows: w.writerow(r)
print(f'group summary: {len(group_rows)} rows')

# ---------- candidate rationale (86 rows) ----------
# Reuse existing manifest + descriptors, add reasoning column with concrete
# scientific hypothesis. Read-only assembly.
picks_desc=list(csv.DictReader(open(V5/'picks_analysis.csv')))
v4_desc=list(csv.DictReader(open(G3/'DFT_shortlist_v3/picks_analysis.csv')))
v4_desc_map={(r['sid'],r['ads'],int(r['idx'])):r for r in v4_desc}
v5_addmap={(r['sid'],r['ads'],int(r['idx'])):r for r in
           csv.DictReader(open(V5/'proposed_additions_only.csv'))}

def dft_hypothesis(r):
    ads=r['ads']; origin=r['origin']; sid=r['sid']
    if origin=='v4-existing':
        # Compose hypothesis based on info available
        if ads=='CO':
            site=v4_desc_map.get((sid,ads,int(r['idx'])),{}).get('site','')
            if site=='physi':
                return 'Confirms MLIP physisorption (no chemical binding); '\
                       'DFT should give weak E_bind and long C-Pd/C-O_surf distance.'
            return f'DFT confirms {site} CO* chemisorption; check final E_bind vs neighbors.'
        if ads=='CH3O':
            site=v4_desc_map.get((sid,ads,int(r['idx'])),{}).get('site','')
            return f'DFT confirms {site} methoxide binding; validate C-O and O-surf distances.'
        # coads
        return 'DFT confirms coads pair energy + geometry; check whether CH3OCO* formation is accessible.'
    else:
        # v5-new
        d=v5_addmap.get((sid,ads,int(r['idx'])),{})
        prio=r['priority']
        if ads=='coads':
            db=d.get('d_bin','')
            if 'reactive' in db:
                return f'Confirm {db} pair is a real TS1 endpoint; DFT should keep d_reactive<4 Å.'
            if db=='thermodynamic':
                return 'Confirm thermodynamic reference E for coadsorption; long-separation baseline.'
            if db=='separated':
                return 'Cross-check separated-pair energy; may drift to thermo bin in DFT.'
        else:
            site=d.get('site',''); region=d.get('region','')
            return f'Confirm {site}/{region} single-ads is genuine minimum, not artifact.'
    return ''

rows=list(csv.DictReader(open(V5A/'combined_summary.csv')))
for r in rows:
    r['dft_hypothesis']=dft_hypothesis(r)

keys=list(rows[0].keys())
with open(OUT/'selection_candidate_rationale.csv','w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=keys); w.writeheader()
    for r in rows: w.writerow(r)
print(f'candidate_rationale: {len(rows)} rows')
