#!/bin/bash

source variables/paths

for dataset in "test"
do
    for dir in $(ls -d ${dataset}/????)
    do
        # PDB name
	    system=$(basename ${dir})

        # CSV file
        csvfile=${dir}/${system}_score.csv
        rm -f ${csvfile}

        echo "name,rmsd_lig,score" >> ${csvfile}

        ligand_crystal=${pdbbind}/${dataset}/${system}/${system}_ligand.mol2

        for ligand in $(ls ${dir}/${system}_ligand-*.pdb)
        do
            name=$(basename $ligand .pdb)
            rmsd=$(${obrms} ${ligand_crystal} ${ligand} | awk  '{print $2}')
            score=$(grep "minimizedAffinity" ${ligand} | awk '{print $3}')

            echo "${name},${rmsd},${score}" >> ${csvfile}
        done
    done
done