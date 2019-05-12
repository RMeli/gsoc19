#!/bin/bash

source path.var

# Store working directory
wd=$PWD

for dataset in "refined"
do
    for dir in $(ls -d ${pdbbind}/${dataset}/1[0-9]??)
    do

 	    # PDB name
	    system=$(basename -a ${dir})

        # Docking working directory
        ddir=${dataset}/${system}
        mkdir -p ${ddir}
        mkdir -p ${ddir}/logs

        # Ligand and receptor
        ligand=${dir}/${system}_ligand.mol2
        receptor=${dir}/${system}_protein.pdb

        # Compute ligand center
        python3.6 ${pscripts}/ligand_center.py ${ligand} \
            -o ${ddir}/center.dat

        # Get center
        c_x=$(cat ${ddir}/center.dat | awk '{print $1}')
        c_y=$(cat ${ddir}/center.dat | awk '{print $2}')
        c_z=$(cat ${ddir}/center.dat | awk '{print $3}')

        # Docking parameters
        source docking.var

        # Run docking
        ${smina} -r ${receptor} -l ${ligand} \
            --flexdist_ligand ${ligand} --flexdist 3 \
            --center_x ${c_x} --center_y ${c_y} --center_z ${c_z} \
            --size_x ${l_box} --size_y ${l_box} --size_z ${l_box} \
            --exhaustiveness ${exhaustiveness} --num_modes ${num_modes} --cpu ${cpu} \
            --out ${ddir}/smina.pdb --out_flex ${ddir}/flex.pdb \
            2>&1 | tee ${ddir}/logs/smina.log

        # Split docking results
        mkdir ${ddir}/dock ${ddir}/flex
        ${obabel} -m -ipdb ${ddir}/smina.pdb -opdb -O ${ddir}/dock/${system}_dock-.pdb \
            2>&1 | tee ${ddir}/logs/obabel_dock.log
        ${obabel} -m -ipdb ${ddir}/flex.pdb -opdb -O ${ddir}/flex/${system}_flex-.pdb \
            2>&1 | tee ${ddir}/logs/obabel_flex.log

    done
done