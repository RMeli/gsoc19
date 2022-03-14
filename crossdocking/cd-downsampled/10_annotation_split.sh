#!/bin/bash
#SBATCH --job-name=sannotation
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
#SBATCH --output=slurm/annotation/%x.%A_%a.out
#SBATCH --error=slurm/annotation/%x.%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --account=stat-ecr
#SBATCH --qos=standard

for annotation in "flex05" "flex1" "flex2" "max2"
do
    python splitannotation.py --receptor ${annotation}
    python splitannotation.py -c --receptor ${annotation}
done
