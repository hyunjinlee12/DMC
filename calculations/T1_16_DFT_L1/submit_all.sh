#!/bin/bash
# T1.16 DFT batch submit
# Review then uncomment desired lines
# All paths relative to project root

cd /home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc

# sbatch -J S1_single_CO_00 --chdir=calculations/T1_16_DFT_L1/S1/single_CO/00_single_CO_rank0_idx00064 scripts/submit_vasp_gpu.sh
# sbatch -J S1_single_CO_01 --chdir=calculations/T1_16_DFT_L1/S1/single_CO/01_single_CO_rank1_idx00099 scripts/submit_vasp_gpu.sh
# sbatch -J S1_single_CO_02 --chdir=calculations/T1_16_DFT_L1/S1/single_CO/02_single_CO_rank2_idx00058 scripts/submit_vasp_gpu.sh
# sbatch -J S1_single_CH3O_00 --chdir=calculations/T1_16_DFT_L1/S1/single_CH3O/00_single_CH3O_rank0_idx00048 scripts/submit_vasp_gpu.sh
# sbatch -J S1_single_CH3O_01 --chdir=calculations/T1_16_DFT_L1/S1/single_CH3O/01_single_CH3O_rank1_idx00331 scripts/submit_vasp_gpu.sh
# sbatch -J S1_single_CH3O_02 --chdir=calculations/T1_16_DFT_L1/S1/single_CH3O/02_single_CH3O_rank2_idx00333 scripts/submit_vasp_gpu.sh
# sbatch -J S1_coads_SetA_00 --chdir=calculations/T1_16_DFT_L1/S1/coads_SetA/00_coads_SetA_rank0_idx00487 scripts/submit_vasp_gpu.sh
# sbatch -J S1_coads_SetA_01 --chdir=calculations/T1_16_DFT_L1/S1/coads_SetA/01_coads_SetA_rank1_idx00767 scripts/submit_vasp_gpu.sh
# sbatch -J S1_coads_SetA_02 --chdir=calculations/T1_16_DFT_L1/S1/coads_SetA/02_coads_SetA_rank2_idx04132 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CO_00 --chdir=calculations/T1_16_DFT_L1/S2/single_CO/00_single_CO_rank0_idx00059 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CO_01 --chdir=calculations/T1_16_DFT_L1/S2/single_CO/01_single_CO_rank1_idx00081 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CO_02 --chdir=calculations/T1_16_DFT_L1/S2/single_CO/02_single_CO_rank2_idx00130 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CO_03 --chdir=calculations/T1_16_DFT_L1/S2/single_CO/03_single_CO_rank3_idx00014 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CO_04 --chdir=calculations/T1_16_DFT_L1/S2/single_CO/04_single_CO_rank4_idx00020 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CH3O_00 --chdir=calculations/T1_16_DFT_L1/S2/single_CH3O/00_single_CH3O_rank0_idx00019 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CH3O_01 --chdir=calculations/T1_16_DFT_L1/S2/single_CH3O/01_single_CH3O_rank1_idx00297 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CH3O_02 --chdir=calculations/T1_16_DFT_L1/S2/single_CH3O/02_single_CH3O_rank2_idx00082 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CH3O_03 --chdir=calculations/T1_16_DFT_L1/S2/single_CH3O/03_single_CH3O_rank3_idx00048 scripts/submit_vasp_gpu.sh
# sbatch -J S2_single_CH3O_04 --chdir=calculations/T1_16_DFT_L1/S2/single_CH3O/04_single_CH3O_rank4_idx00047 scripts/submit_vasp_gpu.sh
# sbatch -J S2_coads_SetA_00 --chdir=calculations/T1_16_DFT_L1/S2/coads_SetA/00_coads_SetA_rank0_idx07672 scripts/submit_vasp_gpu.sh
# sbatch -J S2_coads_SetA_01 --chdir=calculations/T1_16_DFT_L1/S2/coads_SetA/01_coads_SetA_rank1_idx08079 scripts/submit_vasp_gpu.sh
# sbatch -J S2_coads_SetA_02 --chdir=calculations/T1_16_DFT_L1/S2/coads_SetA/02_coads_SetA_rank2_idx10461 scripts/submit_vasp_gpu.sh
# sbatch -J S2_coads_SetA_03 --chdir=calculations/T1_16_DFT_L1/S2/coads_SetA/03_coads_SetA_rank3_idx02538 scripts/submit_vasp_gpu.sh
# sbatch -J S3_single_CO_00 --chdir=calculations/T1_16_DFT_L1/S3/single_CO/00_single_CO_rank0_idx00123 scripts/submit_vasp_gpu.sh
# sbatch -J S3_single_CO_01 --chdir=calculations/T1_16_DFT_L1/S3/single_CO/01_single_CO_rank1_idx00040 scripts/submit_vasp_gpu.sh
# sbatch -J S3_single_CO_02 --chdir=calculations/T1_16_DFT_L1/S3/single_CO/02_single_CO_rank2_idx00020 scripts/submit_vasp_gpu.sh
# sbatch -J S3_single_CH3O_00 --chdir=calculations/T1_16_DFT_L1/S3/single_CH3O/00_single_CH3O_rank0_idx00315 scripts/submit_vasp_gpu.sh
# sbatch -J S3_single_CH3O_01 --chdir=calculations/T1_16_DFT_L1/S3/single_CH3O/01_single_CH3O_rank1_idx00369 scripts/submit_vasp_gpu.sh
# sbatch -J S3_single_CH3O_02 --chdir=calculations/T1_16_DFT_L1/S3/single_CH3O/02_single_CH3O_rank2_idx00126 scripts/submit_vasp_gpu.sh
# sbatch -J S3_coads_SetA_00 --chdir=calculations/T1_16_DFT_L1/S3/coads_SetA/00_coads_SetA_rank0_idx03458 scripts/submit_vasp_gpu.sh
# sbatch -J S3_coads_SetA_01 --chdir=calculations/T1_16_DFT_L1/S3/coads_SetA/01_coads_SetA_rank1_idx03481 scripts/submit_vasp_gpu.sh
# sbatch -J S3_coads_SetA_02 --chdir=calculations/T1_16_DFT_L1/S3/coads_SetA/02_coads_SetA_rank2_idx02834 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_single_CO_00 --chdir=calculations/T1_16_DFT_L1/S3b/single_CO/00_single_CO_rank0_idx00056 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_single_CO_01 --chdir=calculations/T1_16_DFT_L1/S3b/single_CO/01_single_CO_rank1_idx00084 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_single_CO_02 --chdir=calculations/T1_16_DFT_L1/S3b/single_CO/02_single_CO_rank2_idx00041 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_single_CH3O_00 --chdir=calculations/T1_16_DFT_L1/S3b/single_CH3O/00_single_CH3O_rank0_idx00198 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_single_CH3O_01 --chdir=calculations/T1_16_DFT_L1/S3b/single_CH3O/01_single_CH3O_rank1_idx00016 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_single_CH3O_02 --chdir=calculations/T1_16_DFT_L1/S3b/single_CH3O/02_single_CH3O_rank2_idx00260 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_coads_SetA_00 --chdir=calculations/T1_16_DFT_L1/S3b/coads_SetA/00_coads_SetA_rank0_idx02051 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_coads_SetA_01 --chdir=calculations/T1_16_DFT_L1/S3b/coads_SetA/01_coads_SetA_rank1_idx02603 scripts/submit_vasp_gpu.sh
# sbatch -J S3b_coads_SetA_02 --chdir=calculations/T1_16_DFT_L1/S3b/coads_SetA/02_coads_SetA_rank2_idx02754 scripts/submit_vasp_gpu.sh
# sbatch -J S4_single_CO_00 --chdir=calculations/T1_16_DFT_L1/S4/single_CO/00_single_CO_rank0_idx00007 scripts/submit_vasp_gpu.sh
# sbatch -J S4_single_CO_01 --chdir=calculations/T1_16_DFT_L1/S4/single_CO/01_single_CO_rank1_idx00080 scripts/submit_vasp_gpu.sh
# sbatch -J S4_single_CO_02 --chdir=calculations/T1_16_DFT_L1/S4/single_CO/02_single_CO_rank2_idx00033 scripts/submit_vasp_gpu.sh
# sbatch -J S4_single_CH3O_01 --chdir=calculations/T1_16_DFT_L1/S4/single_CH3O/01_single_CH3O_rank1_idx00241 scripts/submit_vasp_gpu.sh
# sbatch -J S4_single_CH3O_90 --chdir=calculations/T1_16_DFT_L1/S4/single_CH3O/90_single_CH3O_manual_atop_Pd scripts/submit_vasp_gpu.sh
# sbatch -J S4_single_CH3O_91 --chdir=calculations/T1_16_DFT_L1/S4/single_CH3O/91_single_CH3O_manual_bridge scripts/submit_vasp_gpu.sh

echo "Done queueing $(grep -c sbatch ${0}) jobs (commented). Uncomment to submit."