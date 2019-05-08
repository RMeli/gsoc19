#!/bin/bash

# PDBBind database
pdbbind=../PDBBind18

# Store working directory
wd=$PWD

for dataset in "refined"
do
    for dir in $(ls -d ${pdbbind}/${dataset}/????)
    do
 	# PDB name
	system=$(basename -a ${dir})
	
	newdir=${dataset}/${system}
	mkdir -p ${newdir}

	# Receptor .pdbqt
	receptor=${system}_protein
	obabel $dir/${receptor}.pdb -xr -O ${newdir}/${receptor}.pdbqt \
	    2>&1 | tee ${newdir}/obabel.log
    done
done
