#!/bin/bash

datasets="refined other"

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

    # TODO: Remove waters beforehand and fix this (there might be PDB with empty chains)
    # Remove empty chains (usually water!)
    flexstr=$(echo $flexstr | sed "s# :[0-9]*##g")
    flexstr=$(echo $flexstr | sed "s#^:[0-9]*##g")

    # Concatenate flexible residues in comma-separated list
    flexlist=""
    for res in ${flexstr}
    do
        flexlist="${flexlist}${res},"
    done

    # Remove final comma
    flexlist=$(echo ${flexlist} | sed 's#,$##g')

    # Store flexible residue list to file for later use
    echo ${flexlist} > ${dataset}/${system}/${system}_flex.info
}

export -f flexinfo

for dataset in ${datasets}
do
    dirs=$(ls -d ${ddir}/${dataset}/????)

    parallel -j ${n_cpus} flexinfo ::: ${dirs} ::: ${dataset}
done