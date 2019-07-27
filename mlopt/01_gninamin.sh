#!/bin/bash

list=lists/test.lst


source variables/paths

mkdir -p test

for dir in $(cat ${list})
do
    system=$(basename ${dir})

    outdir=test/${system}
    mkdir -p ${outdir}

    ligand=${dataroot}/${dir}/${system}_ligand.mol2
    receptor=${dataroot}/${dir}/${system}_protein.pdb

    ${gnina} -l ${ligand} -r ${receptor} \
    --flexdist 3 --flexdist_ligand ${ligand} \
    --cnn default2017 \
    --cnn_scoring --minimize  \
    --gpu \
    --out ${outdir}/${system}_ligand-min.pdb\
    --out_flex ${outdir}/${system}_flex-min.pdb
done