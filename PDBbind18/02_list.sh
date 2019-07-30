#!/bin/bash

# List all paths to PDBbind18 systems

lname=pdbbind18-all

rm -f ${lname}.list

for dataset in "refined" "other"
do
    for dir in $(ls -d ${dataset}/????) 
    do
	system=$(basename ${dir})

	echo "${dataset}/${system}" >> ${lname}.list
    
    done
done