#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=vasp_%j.out
#SBATCH --error=vasp_%j.err

# Usage: sbatch -J <jobname> --chdir=<calcdir> [--ntasks-per-node=N] scripts/submit_vasp.sh [vtst]
# Pass "vtst" as $1 to use VTST-compiled VASP (NEB/dimer)

# --- Clean conda/python from environment (conflicts with NVHPC MPI) ---
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

# --- Threads: 1 thread per MPI rank to avoid oversubscription ---
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# --- VASP binary selection ---
VASP_STD=/home/hyunjin/vasp.6.4.3/bin/vasp_std
VASP_VTST=/home/hyunjin/vasp.6.4.3_vtst/bin/vasp_std

if [ "${1:-}" = "vtst" ]; then
    VASP_BIN="${VASP_VTST}"
    echo "Using VTST VASP: ${VASP_BIN}"
else
    VASP_BIN="${VASP_STD}"
    echo "Using vanilla VASP: ${VASP_BIN}"
fi

NPROCS=${SLURM_NTASKS:-1}
echo "Job: ${SLURM_JOB_NAME} | Dir: $(pwd) | Cores: ${NPROCS} | OMP_THREADS: ${OMP_NUM_THREADS}"
echo "Start: $(date)"

mpirun --bind-to none --map-by slot -np ${NPROCS} ${VASP_BIN}

echo "End: $(date)"
