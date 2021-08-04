#!/bin/bash -l
#SBATCH --job-name=spyrmsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --array=1-4999
##SBATCH --array=1-2971
#SBATCH --array=33-34
##SBATCH --output=slurm/%x.%A_%a.out
#SBATCH --time=12:00:00
##SBATCH --time=00:10:00

# Python
conda activate gsoc
which python

pscripts="../../scripts/python"

file="gnina_cmds_1-4999.txt"
#file="gnina_cmds_5000-7970.txt"

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

outfname="${dir}/${prot}_${rlig}.rmsds"
rm -f $outfname

echo "pocket,protein,ligand,rank,rmsd,obrmsd,flexrmsd,flexobrmsd,fmaxrmsd,score" > ${outfname}

# Initialisation for per-residue RMSD computation
rmsdir=${dir}/${prot}_${rlig}.d
mkdir -p ${rmsdir}
rrfname="${rmsdir}/resrmsd.csv"
echo "rank,res,rmsd,rmsd_obrms" > ${rrfname}

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

    # -----------------------------------------
    # Compute global RMSD for flexible residues
    # -----------------------------------------
    gunzip ${flex}
    refflex_nogz=$(echo ${refflex} | sed 's#.gz##')
    flex_nogz=$(echo ${flex} | sed 's#.gz##')
    echo "python -m spyrmsd ${refflex_nogz} ${flex_nogz}"
    flexrmsd=$(python -m spyrmsd ${refflex_nogz} ${flex_nogz})
    echo "obrms ${refflex_nogz} ${flex_nogz}"
    flexobrmsd=$(obrms ${refflex_nogz} ${flex_nogz} | awk  '{print $3}')
    gzip $(echo ${flex} | sed 's#.gz##')

    # ------------------------
    # Compute per-residue RMSD
    # ------------------------

    # Split residues in different files
    python ${pscripts}/getflex-cd.py ${flex} ${refflex_nogz} --dir ${rmsdir} --rank ${idx}
    
    rmsd_fmax=-1
    for pfname in $(ls ${rmsdir}/pflex-p${idx}-*.pdb)
    do
        cfname=$( echo ${pfname} | sed "s#pflex#cflex#g" )

        # Single residue RMSD
        rmsd_res=$(python -m spyrmsd ${pfname} ${cfname})
        rmsd_res_obrms=$(obrms ${pfname} ${cfname} | awk '{print $3}')
            
        rmsdresname=$( echo $(basename $pfname) | sed "s#pflex-p.*-##g" | sed "s#.pdb##g" )
        echo "${idx},${rmsdresname},${rmsd_res},${rmsd_res_obrms}" >> ${rrfname}

        if (( $(echo "$rmsd_res > $rmsd_fmax" | bc -l) )); then
            rmsd_fmax=${rmsd_res}
        fi
    done

    # grep score from (compressed) SDF file
    score=$(zgrep -A 1 "minimizedAffinity" ${lig} | tail -n 1)

    echo ${pocket},${prot:0:4},${rlig:0:4},${idx},${rmsd},${obrmsd},${flexrmsd},${flexobrmsd},${rmsd_fmax},${score} >> ${outfname}
done

gzip $(echo ${refflex} | sed 's#.gz##')
