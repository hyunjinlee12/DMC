"""Prepare DFT input for 70 candidates in DFT_shortlist_v3.
Directory: calculations/T1_16_DFT_L2/{sid}/{ads}/{rank:02d}_{ads}_idx{orig_idx:05d}/
Files: INCAR, POSCAR, POTCAR, submit_vasp_gpu.sh
"""
import shutil, json
from pathlib import Path

ROOT = Path('/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc')
SHORT = ROOT/'calculations/G3_adsorption/DFT_shortlist_v3'
OUT = ROOT/'calculations/T1_16_DFT_L2'
if OUT.exists(): shutil.rmtree(OUT)
OUT.mkdir()

SDIRS = {'S1':'S1_Pd100','S2':'S2_PdO101_Pd100','S3':'S3_PdO100',
         'S3b':'S3b_PdO100_PdOterm','S4':'S4_PdO2_110'}

# INCAR templates (per surface: Pd metal (S1) uses ISMEAR=1, others oxide use ISMEAR=0)
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

# Template GPU submit
SUBMIT_TEMPLATE = ROOT/'calculations/T1_16_DFT_L1/S1/single_CO/00_single_CO_rank0_idx00064/submit_vasp_gpu.sh'
POTCAR_SRC = {'S1':'/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/T1_16_DFT_L1/S1/single_CO/00_single_CO_rank0_idx00064/POTCAR',
              'S2':'/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/T1_16_DFT_L1/S2/single_CO/00_single_CO_rank0_idx00059/POTCAR',
              'S3':'/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/T1_16_DFT_L1/S3/single_CO/00_single_CO_rank0_idx00123/POTCAR',
              'S3b':'/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/T1_16_DFT_L1/S3b/single_CO/00_single_CO_rank0_idx00056/POTCAR',
              'S4':'/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/calculations/T1_16_DFT_L1/S4/single_CO/00_single_CO_rank0_idx00007/POTCAR'}


picks = json.load(open(SHORT/'summary.json'))
count = 0
for p in picks:
    sid = p['sid']; ads = p['ads']; rank = p['rank']; idx = p['idx']
    # Source POSCAR
    src_dir = SHORT/sid/ads
    src_files = list(src_dir.glob(f"{rank:02d}_*_idx{idx:05d}.vasp"))
    if not src_files:
        print(f'MISSING: {sid}/{ads}/rank {rank} idx {idx}')
        continue
    src_poscar = src_files[0]
    # Dest dir
    dest = OUT/sid/ads/f"{rank:02d}_{ads}_idx{idx:05d}"
    dest.mkdir(parents=True)
    # POSCAR
    shutil.copy(src_poscar, dest/'POSCAR')
    # INCAR
    incar = INCAR_METAL if sid == 'S1' else INCAR_OXIDE
    (dest/'INCAR').write_text(incar)
    # POTCAR (sibling from T1_16 DFT L1)
    shutil.copy(POTCAR_SRC[sid], dest/'POTCAR')
    # submit script
    sub = SUBMIT_TEMPLATE.read_text()
    sub = sub.replace('vasp_%j.out', 'slurm.%j.out').replace('vasp_%j.err', 'slurm.%j.err')
    (dest/'submit_vasp_gpu.sh').write_text(sub)
    count += 1

print(f'Prepared {count} DFT dirs at {OUT}')
