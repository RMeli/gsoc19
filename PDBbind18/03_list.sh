#!/bin/bash

# List all paths to PDBbind18 systems

lname=pdbbind18-test # Complete list
list=pdbbind18-t # Filtered list

rm -f ${lname}.lst
rm -f ${list}.lst

for dataset in "test"
do
    for dir in $(ls -d ${dataset}/????) 
    do
	system=$(basename ${dir})

    # Append to complete list
	echo "${dataset}/${system}" >> ${lname}.lst

    # Absolute path
    path=$PWD

    # Input
    sminame=${path}/${dir}/${system}_ligand.smi

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