#!/bin/bash

list=analysis/icodes.dat

source variables/paths

# Store working directory
wd=$PWD

for path in $(cat ${list})
do
 	# PDB name
	system=$(basename ${path})

    # PDBbind directory
    dir=${pdbbind}/${path}

    # Docking working directory
    ddir=${path}
    mkdir -p ${ddir}
    mkdir -p ${ddir}/logs

    # Ligand and receptor
    ligand=${dir}/${system}_ligand.mol2
    receptor=${dir}/${system}_protein.pdb

    # Docking parameters
    source variables/docking

    # Run docking
    echo "--receptor ${receptor}" | tee ${ddir}/logs/smina.log
    echo "--ligand ${ligand}" | tee ${ddir}/logs/smina.log
    echo "--autobox_add ${autobox_add}" | tee ${ddir}/logs/smina.log
    echo "--exaustiveness ${exhaustiveness}" | tee ${ddir}/logs/smina.log
    echo "--num_modes ${num_modes}" | tee ${ddir}/logs/smina.log
    echo "--cpu ${cpu}" | tee ${ddir}/logs/smina.log
    echo "--out ${ddir}/dock.pdb" | tee ${ddir}/logs/smina.log
    echo "--out_flex ${ddir}/flex.pdb" | tee ${ddir}/logs/smina.log

    ${smina} -r ${receptor} -l ${ligand} \
        --flexdist_ligand ${ligand} --flexdist 3 \
        --autobox_ligand ${ligand} --autobox_add ${autobox_add} \
	    --exhaustiveness ${exhaustiveness} --num_modes ${num_modes} --cpu ${cpu} \
        --out ${ddir}/dock.pdb --out_flex ${ddir}/flex.pdb \
        2>&1 | tee ${ddir}/logs/smina.log
done
