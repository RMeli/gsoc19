#!/bin/bash -l
#SBATCH --job-name=validation
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
#SBATCH --output=slurm/validation.out
#SBATCH --time=12:00:00

# --------------------------------------------------
# Run validation script to make sure that structures
# have been reconstructed correctly.
# --------------------------------------------------

# Python
conda activate gsoc
which python

python validation.py
