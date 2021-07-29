#!/bin/bash
#SBATCH --job-name=makeflex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=devel
#SBATCH --clusters=arc
#SBATCH --array=1-2
##SBATCH --output=slurm/%x.%A_%a.out
#SBATCH --time=00:05:00
##SBATCH --time=12:00:00

# Python
#conda activate gsoc
which python

pscripts="../../scripts/python"

files="gnina_cmds_1-4999.txt"

echo ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${files})

a=( ${line} )
protein=${a[2]}
ligand=${a[-3]}
flexible=${a[-1]}

echo "Protein input: ${protein}"
echo "Ligand output: ${ligand}"
echo "Flexible output: ${flexible}"

ligname=$(echo ${ligand} | sed "s/\.sdf\.gz//")
flexname=$(echo ${flexible} | sed "s/\.sdf\.gz//")

# Split ligand docking
obabel -m -isdf ${ligand} -osdf -O "${ligname}_p.sdf"

# Split flexible side chains docking
obabel -m -isdf ${flexible} -opdb -O "${flexname}_p.pdb"

# Combine flexible and rigid part of the receptor
prot=$(basename ${protein} .pdb)
dir=$(dirname ${protein})

# ${prot} is needed because pockets contain multiple proteins
# Weed to distinguish them between different scripts
for flex in $(ls -d ${dir}/${prot}*flex_p*.pdb)
do
    outfile=$(echo ${flex} | sed 's#\(.*\)flex#\1full#')
    echo ${protein} ${flex} ${outfile}
    python ${pscripts}/makeflex.py ${protein} ${flex} ${outfile}

    # compress files
    gzip -f ${outfile}
    gzip -f ${flex}
done