#!/bin/bash -l
#SBATCH --job-name=rmsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --time=00:10:00
#SBATCH --time=12:00:00

# --------------------------------------------------
# Compute RMSD using obrms and standard GNINA 1.0
# script
# --------------------------------------------------

# Python
conda activate gsoc
which python

sroot="${HOME}/git/GNINA-1.0/analysis_scripts/"

#python3 ${sroot}/obrms_calc.py -i gnina_cmds.txt -d "AOFB"
python3 ${sroot}/obrms_calc.py -i gnina_cmds.txt -d ""