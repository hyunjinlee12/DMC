#!/bin/bash
#SBATCH --partition=h200q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err
CONTAINER=${CONTAINER:-/scratch/taehun1/hyunjin/vasp_vaspsol.sif}
VASP_BIN=${VASP_BIN:-vasp_std}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
srun --mpi=pmix singularity exec --nv "${CONTAINER}" "${VASP_BIN}"
