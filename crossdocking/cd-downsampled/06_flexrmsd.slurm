#!/bin/bash
#SBATCH --job-name=flexrmsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=long
#SBATCH --clusters=arc
##SBATCH --array=1-4999
#SBATCH --array=1-2971
##SBATCH --array=3,253,260,261,308,315,316,557-562,713,714,741-750,770,771,908,918,968,977,1116-1118,1161,1424-1428,1636,1644,1645,1647,1649,1651,1652,1654,1657,1661,1664,1696,1697,1722,1724,2089,2906-2914,2985,2999,3000,3490-3503,3511-3517,3534-3578,3582-3589,3806,3807,3809,4215,4529,4530
##SBATCH --array=3,253,260,261,308,315,316,557-562,713,714,741-750,770,771,908,918,968,977,1116-1118,1161,1424-1428,1636,1644,1645,1647,1649,1651,1652,1654,1657,1661,1664,1696,1697,1722,1724,2089,2906-2914
##SBATCH --array=3541,3542
#SBATCH --output=slurm/flexrmsd/%x.%A_%a.out
#SBATCH --error=slurm/flexrmsd/%x.%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --account=stat-ecr
#SBATCH --qos=standard

# --------------------------------------------------
# Compute RMSDs:
#   - Ligand RMSDs
#   - Flexible residues RMSDs (combined)
#   - Maximum flexible residues RMSDs
# --------------------------------------------------

# Python
#conda activate gsoc
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
echo "Reference ligand: ${refligand}"

# ${prot} and ${rlig} are needed because pockets contain multiple pairs
# Weed to distinguish them between different scripts
prot=$(basename ${protein} .pdb)
rlig=$(basename ${refligand} .sdf)
dir=$(dirname ${protein})

echo "Protein: ${prot}"
echo "Ligand: ${rlig}"

# Reference for flexible residues
refflex=$(echo ${flexible} | sed "s#_flex.pdb.gz#_flex_crys.pdb.gz#")
gunzip ${refflex}

pocket=$(basename $(dirname ${dir}))

outfname="${dir}/${prot}_${rlig}.rmsds2"
rm -f $outfname

echo "pocket,protein,ligand,rank,rmsd,obrmsd,flexrmsd,fmaxrmsd,score" > ${outfname}

for flex in $(ls -d ${dir}/${prot}_${rlig}*flex_p*.pdb.gz)
do
    idx=$(basename ${flex} .pdb.gz | sed 's#.*flex_p##')

    # Change name and extension from flex to lig
    lig=$(echo ${flex} | sed 's#\(.*\)_flex#\1#')
    lig=$(echo ${lig} | sed 's#.pdb#.sdf#')

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

    # ----------------------------------
    # Compute RMSD for flexible residues
    # ----------------------------------
    # Name of full receptor pose
    recpose=$(echo ${flex} | sed 's#_flex_p#_full_p#')
    # Name of full cognate receptor (same name as ligand)
    # Need to get receptors from full database (some missing from downsampled)
    # The full database has been renamed wierbowski_cd to avoid confusion
    reccognate="wierbowski_cd/${pocket}/${rlig:0:4}_PRO.pdb"
    echo "python ${pscripts}/flexrmsd.py ${recpose} ${reccognate} ${flex}"
    RMSDs=$(python ${pscripts}/flexrmsd.py ${recpose} ${reccognate} ${flex})
    flexrmsd=$(echo ${RMSDs} | awk '{print $1}')
    rmsd_fmax=$(echo ${RMSDs} | awk '{print $2}')

    # grep score from (compressed) SDF file
    score=$(zgrep -A 1 "minimizedAffinity" ${lig} | tail -n 1)

    echo ${pocket},${prot:0:4},${rlig:0:4},${idx},${rmsd},${obrmsd},${flexrmsd},${rmsd_fmax},${score}
    echo ${pocket},${prot:0:4},${rlig:0:4},${idx},${rmsd},${obrmsd},${flexrmsd},${rmsd_fmax},${score} >> ${outfname}
done

gzip $(echo ${refflex} | sed 's#.gz##')
