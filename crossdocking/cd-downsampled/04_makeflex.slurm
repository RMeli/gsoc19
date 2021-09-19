#!/bin/bash -l
#SBATCH --job-name=makeflex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --array=1-4999
##SBATCH --array=1-2971
#SBATCH --array=86-89
#SBATCH --output=slurm/makeflex/%x.%A_%a.out
#SBATCH --time=12:00:00

# --------------------------------------------------
# Build whole protein from original structure and
# flexible residues.
# --------------------------------------------------

# Python
conda activate gsoc
which python

pscripts="../../scripts/python"

#cmds="gnina_cmds_1-4999.txt"
cmds="gnina_cmds_5000-7970.txt"

echo ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${cmds})

a=( ${line} )
protein=${a[2]}
refligand=${a[4]}
ligand=${a[-3]}
flexible=${a[-1]}

echo "Protein input: ${protein}"
echo "Reference ligand: ${refligand}"
echo "Ligand output: ${ligand}"
echo "Flexible output: ${flexible}"

ligname=$(echo ${ligand} | sed "s/\.sdf\.gz//")
flexname=$(echo ${flexible} | sed "s/\.pdb\.gz//")

# Split ligand docking and compress them
obabel -m -isdf ${ligand} -osdf -O "${ligname}_p.sdf"
gzip -f ${ligname}_p*.sdf


# Split flexible side chains docking
# Compression happens after makeflex
# _p0 (minimised crtystal) needs to be decompressed (re-compressed after makeflex)
obabel -m -ipdb ${flexible} -opdb -O "${flexname}_p.pdb"
gunzip -f ${flexname}_p0.pdb.gz


# Combine flexible and rigid part of the receptor
prot=$(basename ${protein} .pdb)
rlig=$(basename ${refligand} .sdf)
dir=$(dirname ${protein})

# ${prot} is needed because pockets contain multiple proteins
# Weed to distinguish them between different scripts
for flex in $(ls -d ${dir}/${prot}_${rlig}*flex_p*.pdb)
do
    outfile=$(echo ${flex} | sed 's#\(.*\)flex#\1full#')
    echo ${protein} ${flex} ${outfile}
    python ${pscripts}/makeflex.py ${protein} ${flex} ${outfile}

    # compress files
    gzip -f ${outfile}
    gzip -f ${flex}
done