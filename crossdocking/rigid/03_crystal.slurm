#!/bin/bash
#SBATCH --job-name=crystal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
#SBATCH --array=1-4999
##SBATCH --array=1-2971
#SBATCH --output=slurm/crystal/%x.%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --account=stat-ecr

# --------------------------------------------------
# Extract residues treated as flexible from crystal
# (for RMSD calculation) and oprimise crysytal
# structure (with flexible residues) to provide the
# best possible pose.
# --------------------------------------------------

# Get flexible residues position from crystal
# Optimise crystal structure with GNINA

cmds="gnina_cmds_1-4999.txt"
#cmds="gnina_cmds_5000-7970.txt"

container="${DATA}/gsoc19/devtools/containers/gnina.sif"

echo ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}

cmd=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${cmds})

# Get flexible residues from crystal pose with --score-only
ccmd=$(echo ${cmd} | sed "s#flexdistNone\.sdf.\gz#flexdistNone_crys\.sdf.\gz#")
echo $ccmd
singularity run -B ${PDW}:${PWD} ${container} ${ccmd} --score_only

# Get minimised crystal pose
mcmd=$(echo ${cmd} | sed "s#flexdistNone\.sdf.\gz#flexdistNone_p0\.sdf.\gz#")
echo $mcmd
singularity run -B ${PDW}:${PWD} ${container} ${mcmd} --minimize
