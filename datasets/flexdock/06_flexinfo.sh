#!/bin/bash

datasets="test"

n_cpus=8

source variables/paths

flexinfo(){
    dir=$1
    dataset=$2

    # PDB name
    system=$(basename ${dir})

    # Get string of flexible residues from SMINA log
    flexstr=$(grep "Flexible residues:" ${dir}/logs/smina.log)
    flexstr=$(echo $flexstr | sed "s#Flexible residues:##g")

    # Concatenate flexible residues in comma-separated list
    flexlist=""
    for res in ${flexstr}
    do
        flexlist="${flexlist}${res},"
    done

    # Remove final comma
    flexlist=$(echo ${flexlist} | sed 's#,$##g')

    # Store flexible residue list to file for later use
    echo ${flexlist} > ${dataset}/${system}/flex.info
}

export -f flexinfo

for dataset in ${datasets}
do
    dirs=$(ls -d ${ddir}/${dataset}/????)

    parallel -j ${n_cpus} flexinfo ::: ${dirs} ::: ${dataset}
done