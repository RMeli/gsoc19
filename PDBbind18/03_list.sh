#!/bin/bash

# List all PDBbind18 systems
# List all PDBbind18 systems with valid SMILES and molecular weight lower than 1000
# Create PDBFILES file used by gnina clustering.py script (for cross-validation)

listall=pdbbind18-all # Complete list
list=pdbbind18 # Filtered list
PDBFILES=pdbfiles # PDBFILES list

max_mass=1000 # Max allowed ligand mass

rm -f ${listall}.lst
rm -f ${list}.lst
rm -f ${PDBFILES}.lst

for dataset in "refined" "other"
do
    for dir in $(ls -d ${dataset}/????) 
    do
	system=$(basename ${dir})

    # Append to complete list
	echo "${dataset}/${system}" >> ${listall}.lst

    # Absolute path
    path=$PWD

    # Input
    recname=${path}/${dir}/${system}_protein.pdb
    sminame=${path}/${dir}/${system}_ligand.smi

    # PDBFILES 
    echo ${system} ${recname} ${sminame} >> ${PDBFILES}.lst

    # Heavy atoms molecular weight
    hamw=$(python ../scripts/python/molweight.py ${sminame} | awk '{print $2}')

    # Append to filtered list
    if (( $(echo "${hamw} < ${max_mass}" | bc -l) )) 
    then
        echo ${dataset}/${system} >> ${list}.lst
    else
        echo "Discarding ${dataset}/${system} (mw = ${hamw})..."
    fi
    
    done
done