#!/bin/bash

# smina exec
smina=~/software/Autodock/Smina/smina.static

# scripts
scripts=../scripts

# Store working directory
wd=$PWD

for dataset in "refined"
do
    for dir in $(ls -d ${dataset}/????)
    do

 	# PDB name
	system=$(basename -a ${dir})

    # Compute ligand center
    python ${scripts}/python/ligand_center.py ${dir}/${system}_ligand.pdbqt \
        -o ${dir}/lig_center.dat

    # Get center
    c_x=$(cat ${dir}/lig_center.dat | awk '{print $1}')
    c_y=$(cat ${dir}/lig_center.dat | awk '{print $2}')
    c_z=$(cat ${dir}/lig_center.dat | awk '{print $3}')

    # Docking parameters
    l_box=23.5

    # Run docking
    ${smina} -r ${dir}/${system}_protein.pdbqt -l ${dir}/${system}_ligand.pdbqt \
        --flexdist_ligand ${dir}/${system}_ligand.pdbqt --flexdist 3 \
        --center_x ${c_x} --center_y ${c_y} --center_z ${c_z} \
        --size_x ${l_box} --size_y ${l_box} --size_z ${l_box} \
        --exhaustiveness 8 --num_modes 20 --cpu 6 \
        --log ${dir}/logs/smina.log
		
    done
done