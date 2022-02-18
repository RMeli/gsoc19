#!/bin/bash
#SBATCH --job-name=annotation
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

# No receptor annotation
python annotation.py
python annotation.py -c

for annotation in "flex1" "flex2" "max2"
do
    python annotation.py --receptor ${annotation}
    python annotation.py -c --receptor ${annotation}
done
