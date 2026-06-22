#!/bin/bash
#SBATCH --job-name=MACE_phase3
#SBATCH --partition=debug
#SBATCH --gres=gpu:rtx6000ada:1
#SBATCH --time=24:00:00
#SBATCH --output=calculations/G3_adsorption/MLIP_phase3_slurm_%j.log
#SBATCH --error=calculations/G3_adsorption/MLIP_phase3_slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# T1.14 Phase 3: MACE-MH + D3 ranking of co-ads SetTS + SetB subsample
# Expected: ~7.5 hr GPU on RTX 6000 Ada
# Inputs:  calculations/G3_adsorption/{surface}/coads_guide/{SetTS,SetB}.traj
# Outputs: calculations/G3_adsorption/{surface}/MLIP_phase3/

set -e

echo "=========================================="
echo "  T1.14 Phase 3 — MACE-MH + D3 ranking"
echo "  SLURM job:  $SLURM_JOB_ID"
echo "  Node:       $SLURMD_NODENAME"
echo "  GPU:        $SLURM_JOB_GPUS"
echo "  Start:      $(date)"
echo "=========================================="

cd /home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc

# SLURM provides CUDA_VISIBLE_DEVICES based on --gres; the script will use it.
# (No need to set manually like the nohup phase did.)

# Run via conda env pddmc (mace-torch 0.3.16 + torch-dftd 0.5.3 + cueq-cuda-12)
conda run -n pddmc python scripts/run_mace_phase3.py

echo "=========================================="
echo "  Phase 3 finished: $(date)"
echo "=========================================="
