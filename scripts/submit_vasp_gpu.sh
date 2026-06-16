#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=7-00:00:00
#SBATCH --output=vasp_%j.out
#SBATCH --error=vasp_%j.err

# Usage: sbatch -J <jobname> --chdir=<calcdir> scripts/submit_vasp_gpu.sh

# --- Clean conda/python from environment ---
unset PYTHONPATH PYTHONHOME CONDA_PREFIX CONDA_DEFAULT_ENV
unset CONDA_SHLVL CONDA_PROMPT_MODIFIER
export LD_LIBRARY_PATH=""

# --- NVHPC compiler + MPI ---
export NVHPC=$HOME/nvhpc
export NVARCH=Linux_x86_64
export NVVERSION=25.9
export PATH=$NVHPC/$NVARCH/$NVVERSION/compilers/bin:$PATH
export PATH=$NVHPC/$NVARCH/$NVVERSION/comm_libs/mpi/bin:$PATH

export LD_LIBRARY_PATH=$NVHPC/$NVARCH/$NVVERSION/comm_libs/mpi/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVHPC/$NVARCH/$NVVERSION/compilers/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVHPC/$NVARCH/$NVVERSION/compilers/extras/qd/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/fftw/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/cuda-13.0/lib64:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1

VASP_BIN=/home/hyunjin/vasp.6.4.3/bin/vasp_std

NPROCS=${SLURM_NTASKS:-1}
echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | MPI ranks: ${NPROCS} | GPUs: ${SLURM_GPUS_ON_NODE:-1}"
echo "Start: $(date)"

mpirun --bind-to none -np ${NPROCS} ${VASP_BIN}

echo "End: $(date)"
