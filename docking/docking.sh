#!/bin/bash

source variables/paths

# Store working directory
wd=$PWD

for dataset in "refined"
do
    for dir in $(ls -d ${pdbbind}/${dataset}/[1]???)
    do

 	# PDB name
	system=$(basename ${dir})

        # Docking working directory
        ddir=${dataset}/${system}
        mkdir -p ${ddir}
        mkdir -p ${ddir}/logs

        # Ligand and receptor
        ligand=${dir}/${system}_ligand.mol2
        receptor=${dir}/${system}_protein.pdb

        # Docking parameters
        source variables/docking

        # Run docking
        ${smina} -r ${receptor} -l ${ligand} \
            --flexdist_ligand ${ligand} --flexdist 3 \
            --autobox_ligand ${ligand} --autobox_add ${autobox_add} \
	        --exhaustiveness ${exhaustiveness} --num_modes ${num_modes} --cpu ${cpu} \
            --out ${ddir}/dock.pdb --out_flex ${ddir}/flex.pdb \
            2>&1 | tee ${ddir}/logs/smina.log

    done
done
