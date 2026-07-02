"""Submit DFT jobs for one surface at a time. Usage: python submit_dft_L2_surface.py S1"""
import subprocess, sys
from pathlib import Path

if len(sys.argv) < 2:
    print("Usage: python submit_dft_L2_surface.py <S1|S2|S3|S3b|S4>")
    sys.exit(1)
sid = sys.argv[1]

L2 = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/T1_16_DFT_L2')
surf_dir = L2 / sid
if not surf_dir.exists():
    print(f"No dir {surf_dir}"); sys.exit(1)

submitted = []
for ads_dir in sorted(surf_dir.iterdir()):
    if not ads_dir.is_dir(): continue
    for cand_dir in sorted(ads_dir.iterdir()):
        if not cand_dir.is_dir(): continue
        submit = cand_dir/'submit_vasp_gpu.sh'
        if not submit.exists():
            print(f"MISSING submit: {cand_dir}"); continue
        # job name: sid_ads_rank
        jn = f"{sid}_{ads_dir.name}_{cand_dir.name.split('_')[0]}"
        r = subprocess.run(['sbatch','-J',jn,'--chdir',str(cand_dir),
                            str(submit)], capture_output=True, text=True)
        if r.returncode == 0:
            jid = r.stdout.strip().split()[-1]
            submitted.append((jid, jn))
            print(f"  {jid}  {jn}")
        else:
            print(f"  FAIL {jn}: {r.stderr[:100]}")

print(f'\nSubmitted {len(submitted)} jobs for {sid}')
