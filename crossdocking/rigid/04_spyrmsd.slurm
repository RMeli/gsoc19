#!/bin/bash -l
#SBATCH --job-name=spyrmsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --array=1-4999
##SBATCH --array=1-2971
#SBATCH --array=1784
#SBATCH --output=slurm/spyrmsd/%x.%A_%a.out
#SBATCH --time=06:00:00
#SBATCH --account=stat-ecr

# --------------------------------------------------
# Compute RMSDs:
#   - Ligand RMSDs
# --------------------------------------------------

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
refligand=${a[4]}
ligand=${a[-1]}

echo "Protein input: ${protein}"
echo "Ligand output: ${ligand}"
echo "Reference ligand: ${refligand}"

# ${prot} and ${rlig} are needed because pockets contain multiple pairs
# Weed to distinguish them between different scripts
prot=$(basename ${protein} .pdb)
rlig=$(basename ${refligand} .sdf)
dir=$(dirname ${protein})

echo "Protein: ${prot}"
echo "Ligand: ${rlig}"

# Split ligand docking and compress them
ligname=$(echo ${ligand} | sed "s/\.sdf\.gz//")
obabel -m -isdf ${ligand} -osdf -O "${ligname}_p.sdf"
gzip -f ${ligname}_p*.sdf

pocket=$(basename ${dir})

outfname="${dir}/${prot}_${rlig}.rmsds"
rm -f $outfname

echo "pocket,protein,ligand,rank,rmsd,obrmsd,score" > ${outfname}

for lig in $(ls -d ${dir}/${prot}_${rlig}*_p*.sdf.gz)
do
    idx=$(basename ${lig} .sdf.gz | sed 's#.*_p##')

    # -------------------
    # Compute ligand RMSD
    # -------------------
    gunzip ${lig}
    reflig_nogz=$(echo ${refligand} | sed 's#.gz##')
    lig_nogz=$(echo ${lig} | sed 's#.gz##')
    echo "python -m spyrmsd ${reflig_nogz} ${lig_nogz}"
    rmsd=$(python -m spyrmsd ${reflig_nogz} ${lig_nogz})
    echo "obrms ${reflig_nogz} ${lig_nogz}"
    obrmsd=$(obrms ${reflig_nogz} ${lig_nogz} | awk '{print $3}')
    gzip ${lig_nogz}

    # grep score from (compressed) SDF file
    score=$(zgrep -A 1 "minimizedAffinity" ${lig} | tail -n 1)

    echo ${pocket},${prot:0:4},${rlig:0:4},${idx},${rmsd},${obrmsd},${score} >> ${outfname}
done