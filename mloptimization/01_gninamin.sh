#!/bin/bash

list=lists/test.lst


source variables/paths

mkdir -p optimal

for dir in $(cat ${list})
do
    system=$(basename ${dir})

    outdir=optimal/${system}
    mkdir -p ${outdir}

    ligand=${dataroot}/${dir}/${system}_ligand.mol2
    receptor=${dataroot}/${dir}/${system}_protein.pdb

    ${gnina} -l ${ligand} -r ${receptor} \
    --cnn default2017 \
    --cnn_scoring --minimize  \
    --gpu \
    --out ${outdir}/${system}_ligand.pdb
done