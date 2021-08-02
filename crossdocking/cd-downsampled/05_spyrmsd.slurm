#!/bin/bash -l
#SBATCH --job-name=spyrmsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --array=1-4999
#SBATCH --array=1-2971
##SBATCH --array=1-2
#SBATCH --output=slurm/%x.%A_%a.out
#SBATCH --time=12:00:00

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
ligand=${a[-3]}
flexible=${a[-1]}

echo "Protein input: ${protein}"
echo "Ligand output: ${ligand}"
echo "Flexible output: ${flexible}"

# ${prot} is needed because pockets contain multiple proteins
# Weed to distinguish them between different scripts
prot=$(basename ${protein} .pdb)
rlig=$(basename ${refligand} .sdf)
dir=$(dirname ${protein})

# Reference for flexible residues
refflex=$(echo ${flexible} | sed "s#_flex.pdb.gz#_flex_crys.pdb.gz#")
gunzip ${refflex}

pocket=$(basename $(dirname ${dir}))

outfname="${dir}/${prot}_${rlig}.rmsds"
rm -f $outfname

for flex in $(ls -d ${dir}/${prot}*flex_p*.pdb.gz)
do
    idx=$(basename ${flex} .pdb.gz | sed 's#.*flex_p##')

    # Change name and extension from flex to lig
    lig=$(echo ${flex} | sed 's#\(.*\)_flex#\1#')
    lig=$(echo ${lig} | sed 's#.pdb#.sdf#')

    gunzip ${lig}
    rmsd=$(python -m spyrmsd $(echo ${refligand} | sed 's#.gz##') $(echo ${lig} | sed 's#.gz##'))
    obrmsd=$(python -m spyrmsd $(echo ${refligand} | sed 's#.gz##') $(echo ${lig} | sed 's#.gz##'))
    gzip $(echo ${lig} | sed 's#.gz##')

    gunzip ${flex}
    flexrmsd=$(python -m spyrmsd $(echo ${refflex} | sed 's#.gz##') $(echo ${flex} | sed 's#.gz##'))
    flexobrmsd=$(python -m spyrmsd $(echo ${refflex} | sed 's#.gz##') $(echo ${flex} | sed 's#.gz##'))
    gzip $(echo ${flex} | sed 's#.gz##')
    

    # grep score from (compressed) SDF file
    score=$(zgrep -A 1 "minimizedAffinity" ${lig} | tail -n 1)

    echo ${pocket},${prot:0:4},${rlig:0:4},${idx},${rmsd},${obrmsd},${flexrmsd},${flexobrmsd},${score} >> ${outfname}
done

gzip $(echo ${refflex} | sed 's#.gz##')
