#!/bin/bash -l
#SBATCH --job-name=molcache
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
#SBATCH --output=slurm/molcache.out
#SBATCH --time=12:00:00

gscripts="$HOME/git/gninascripts"
typedir="carlos_cd"

# The union of the CV test sets is the full dataset
# Folds are obtained with GroupKFold
cat files/test_*.types > files/all.types

python ${gscripts}/create_caches2.py \
    "files/all.types" \
    -c 1 -d . \
    --recmolcache "rec.molcache2" \
    --ligmolcache "lig.molcache2"