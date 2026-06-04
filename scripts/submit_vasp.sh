#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=vasp_%j.out
#SBATCH --error=vasp_%j.err

# Usage: sbatch -J <jobname> --chdir=<calcdir> [--ntasks-per-node=N] scripts/submit_vasp.sh [vtst]
# Pass "vtst" as $1 to use VTST-compiled VASP (NEB/dimer)

source /home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc/scripts/env_config.sh

if [ "$1" = "vtst" ]; then
    VASP_BIN="${VASP_VTST}"
    echo "Using VTST VASP: ${VASP_BIN}"
else
    VASP_BIN="${VASP_STD}"
    echo "Using vanilla VASP: ${VASP_BIN}"
fi

echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | Cores: ${SLURM_NTASKS}"
echo "Start: $(date)"

mpirun -np ${SLURM_NTASKS} ${VASP_BIN}

echo "End: $(date)"
