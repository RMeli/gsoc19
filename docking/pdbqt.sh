#!/bin/bash

# PDBBind database
pdbbind=../PDBbind18

# Store working directory
wd=$PWD

for dataset in "refined"
do
    for dir in $(ls -d ${pdbbind}/${dataset}/1[0-9]??)
    do
 	# PDB name
	system=$(basename -a ${dir})
	
	newdir=${dataset}/${system}
	mkdir -p ${newdir} && mkdir -p ${newdir}/logs

	# Receptor .pdbqt
	ligand=${system}_ligand
	obabel $dir/${ligand}.mol2 -xr -O ${newdir}/${ligand}.pdbqt \
	    2>&1 | tee ${newdir}/logs/obabel_lig.log

	# Receptor .pdbqt
	receptor=${system}_protein
	obabel $dir/${receptor}.pdb -xr -O ${newdir}/${receptor}.pdbqt \
	    2>&1 | tee ${newdir}/logs/obabel_rec.log
		
    done
done