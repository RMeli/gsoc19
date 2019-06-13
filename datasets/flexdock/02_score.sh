#!/bin/bash

source variables/paths

csv_header="name,system,rank,rmsd_lig,score"

allscores="allscores.csv"
rm -f ${allscores}
echo ${csv_header} >> ${allscores}

for dataset in "refined"
do
    for dir in $(ls -d ${dataset}/????)
    do
        # PDB name
	    system=$(basename ${dir})

        # CSV file
        csvfile=${dir}/${system}_score.csv
        rm -f ${csvfile}

        echo ${csv_header} >> ${csvfile}

        ligand_crystal=${pdbbind}/${dataset}/${system}/${system}_ligand.mol2

        for ligand in $(ls ${dir}/${system}_ligand-*.pdb)
        do
            name=$(basename $ligand .pdb)
            rank=$( echo $name | sed "s#${system}_ligand-##g" )
            rmsd=$(${obrms} ${ligand_crystal} ${ligand} | awk  '{print $2}')
            score=$(grep "minimizedAffinity" ${ligand} | awk '{print $3}')

            info="${name},${system},${rank},${rmsd},${score}"
            echo ${info} >> ${csvfile}
            echo ${info} >> ${allscores}
        done
    done
done