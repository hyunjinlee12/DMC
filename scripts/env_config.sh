#!/bin/bash
# Project environment config — source this in job scripts and interactive sessions
# Usage: source scripts/env_config.sh

export PROJECT_ROOT="/home/hyunjin/CLAUDE/Pd_DMC/research-pd-dmc"

# VASP
export VASP_STD="/home/hyunjin/vasp.6.4.3/bin/vasp_std"
export VASP_GAM="/home/hyunjin/vasp.6.4.3/bin/vasp_gam"
export VASP_NCL="/home/hyunjin/vasp.6.4.3/bin/vasp_ncl"
export VASP_VTST="/home/hyunjin/vasp.6.4.3_vtst/bin/vasp_std"

# POTCAR (ASE uses this)
export VASP_PP_PATH="/home/hyunjin/POTENTIAL"

# VTST scripts
export VTST_SCRIPTS="/home/hyunjin/VTST/vtstscripts-1040"
export PATH="${VTST_SCRIPTS}:${PATH}"

# Conda
export CONDA_ENV="pddmc"
