#!/bin/bash
#SBATCH --job-name=docking
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --array=1-4999
#SBATCH --array=1-2971
#SBATCH --output=slurm/%x.%A_%a.out
#SBATCH --time=12:00:00

# --------------------------------------------------
# Perform docking with flexible residues using GNINA
# --------------------------------------------------

#cmds="gnina_cmds_1-4999.txt"
cmds="gnina_cmds_5000-7970.txt"

container="${DATA}/gsoc19/devtools/containers/gnina.sif"

echo ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}

cmd=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${cmds})

singularity run -B ${PDW}:${PWD} ${container} ${cmd}
