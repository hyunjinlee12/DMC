"""Compare MOST STABLE configurations between MACE and SevenNet (N=20).
For each surface×ads: top-1 (lowest E) from each MLIP — same idx? similar d_min?"""
import json
from pathlib import Path
from collections import defaultdict

OUT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/reports/mlip_compare/statistical/N20')
mace = json.load(open(OUT / 'mace_results.json'))
seve = json.load(open(OUT / 'sevenet_results.json'))

def group(recs):
    g = defaultdict(list)
    for r in recs:
        g[(r['surface'], r['ads'])].append(r)
    return g
mace_g = group(mace); seve_g = group(seve)

surfaces = ['S1','S2','S3','S3b','S4']
print(f"{'Surf':<5} {'ads':<6} | {'MACE top-1':<35} | {'SevenNet top-1':<35} | match?")
print('-'*120)
for sid in surfaces:
    for ads in ['CO','CH3O']:
        m = sorted(mace_g[(sid,ads)], key=lambda r: r['E'])
        s = sorted(seve_g[(sid,ads)], key=lambda r: r['E'])
        if not m or not s: continue
        m1, s1 = m[0], s[0]
        same_idx = m1['orig_idx'] == s1['orig_idx']
        d_diff = abs(m1['d_min'] - s1['d_min'])
        # Rank overlap top-3
        m_top3 = {r['orig_idx'] for r in m[:3]}
        s_top3 = {r['orig_idx'] for r in s[:3]}
        overlap = len(m_top3 & s_top3)
        marker = '✓' if same_idx else ('~' if d_diff < 0.3 else '✗')
        print(f"{sid:<5} {ads:<6} | idx={m1['orig_idx']:<4} d={m1['d_min']:.2f}Å E={m1['E']:+.2f}eV   "
              f"| idx={s1['orig_idx']:<4} d={s1['d_min']:.2f}Å E={s1['E']:+.2f}eV   "
              f"| {marker} (top-3 overlap: {overlap}/3)")
