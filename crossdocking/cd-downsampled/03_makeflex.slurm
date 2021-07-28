#!/bin/bash -l
#SBATCH --job-name=makeflex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --array=1-4999
#SBATCH --array=1-2971
#SBATCH --output=slurm/%x.%A_%a.out
#SBATCH --time=01:00:00

# Python
conda activate gsoc
which python

pscripts="../../scripts/python"

#file="gnina_cmds_1-4999.txt"
file="gnina_cmds_5000-7970.txt"

echo ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${file})

a=( ${line} )
protein=${a[2]}
ligand=${a[-3]}
flexible=${a[-1]}

echo "Protein input: ${protein}"
echo "Ligand output: ${ligand}"
echo "Flexible output: ${flexible}"

ligname=$(echo ${ligand} | sed "s/\.sdf\.gz//")
flexname=$(echo ${flexible} | sed "s/\.pdb\.gz//")

# Split ligand docking
obabel -m -isdf ${ligand} -osdf -O "${ligname}_p.sdf"
gzip -f ${ligname}_p*.sdf

# Split flexible side chains docking
# Do not compress now, compress later
obabel -m -ipdb ${flexible} -opdb -O "${flexname}_p.pdb"

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
