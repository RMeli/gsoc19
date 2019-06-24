#!/bin/bash

source variables/paths

csv_header="system,rank,rmsd_lig,rmsd_flex,score"

allscores="analysis/allscores.csv"
rm -f ${allscores}
echo ${csv_header} > ${allscores}

for dataset in "test"
do
    for dir in $(ls -d ${dataset}/????)
    do
        echo ${dir}

        # PDB name
	    system=$(basename ${dir})

        # CSV file
        csvfile=${dir}/${system}_score.csv
        rm -f ${csvfile}

        echo ${csv_header} > ${csvfile}

        ligand_crystal=${pdbbind}/${dataset}/${system}/${system}_ligand.mol2
        protein_crystal=${pdbbind}/${dataset}/${system}/${system}_protein.pdb

        for ligand in $(ls ${dir}/${system}_ligand-*.pdb)
        do
            # Ligand name and rank
            name=$(basename $ligand .pdb)
            rank=$( echo $name | sed "s#${system}_ligand-##g" )

            # Ligand RMSD (with OpenBabel)
            rmsd_lig=$(${obrms} ${ligand_crystal} ${ligand} | awk  '{print $2}')

            # Flexible residues RMSD (with MDAnalysis)
            flex=${dir}/${system}_flex-${rank}.pdb
            protein=${dir}/${system}_protein-${rank}.pdb
            rmsd_flex=$(python3.6 ${pscripts}/flexrmsd.py ${flex} ${protein} ${protein_crystal})

            # Score (from ligand file)
            score=$(grep "minimizedAffinity" ${ligand} | awk '{print $3}')

            info="${name},${system},${rank},${rmsd_lig},${rmsd_flex},${score}"
            echo ${info} >> ${csvfile}
            echo ${info} >> ${allscores}
        done
    done
done