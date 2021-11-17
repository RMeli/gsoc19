#!/bin/bash -l
#SBATCH --job-name=gtyper
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --array=1-4999
#SBATCH --array=1-2971
#SBATCH --output=slurm/gtyper/%x.%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --account=stat-ecr

# --------------------------------------------------
# Create GNINA type files for all systems
# --------------------------------------------------

#files="gnina_cmds_1-4999.txt"
files="gnina_cmds_5000-7970.txt"

container="${DATA}/gsoc19/devtools/containers/gnina.sif"

echo ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${files})

a=( ${line} )
protein=${a[2]}
refligand=${a[4]}
ligand=${a[-3]}

echo "Protein input: ${protein}"
echo "Reference ligand: ${refligand}"
echo "Ligand output: ${ligand}"

prot=$(basename ${protein} .pdb)
rlig=$(basename ${refligand} .sdf)
dir=$(dirname ${protein})

# ${prot} is needed because pockets contain multiple proteins
# Weed to distinguish them between different scripts
for lig in $(ls -d ${dir}/${prot}_${rlig}*_p*.sdf.gz)
do
    outfile=$(echo ${lig} | sed 's#\.sdf.\gz#.gninatypes#g')
    
    # Run GNINA typer
    singularity run -B ${PDW}:${PWD} ${container} gninatyper \
        ${lig} ${outfile}
done

# ${prot} is needed because pockets contain multiple proteins
# Weed to distinguish them between different scripts
for rec in $(ls -d ${dir}/${prot}.pdb)
do
    outfile=$(echo ${rec} | sed 's#\.pdb#.gninatypes#g')
    
    # Run GNINA typer
    singularity run -B ${PDW}:${PWD} ${container} gninatyper \
        ${rec} ${outfile}
done
