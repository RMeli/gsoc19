#!/bin/bash

list=lists/mintest.lst

optdir=$1

if [[ $optdir == "" ]]
then
  echo "OPTDIR must be specified."
  exit
fi

cd ${optdir}

source ../variables/paths
source opt.vars

IFS='
'
for line in $(cat ../${list})
do
    dir=$(echo $line | awk '{print $1}')
    rank=$(echo $line | awk '{print $2}')

    system=$(basename ${dir})

    outdir=minimized/${system}
    mkdir -p ${outdir}

    ligand=${dataroot}/${dir}/${system}_ligand-${rank}.pdb
    receptor=${dataroot}/${dir}/${system}_protein-${rank}.pdb
    flex=${dataroot}/${dir}/${system}_flex-${rank}.pdb

    echo "ligand = ${ligand}"
    echo "receptor = ${receptor}"

    # TODO: Remove this copy
    cp ${ligand} ${outdir}
    cp ${flex} ${outdir}

    ${gnina} -l ${ligand} -r ${receptor} \
    --flexres $(cat ${dataroot}/${dir}/${system}_flex.info) \
    --cnn_model *.model \
    --cnn_weights ${weights} \
    --cnn_scoring --minimize  \
    --gpu \
    --out ${outdir}/${system}_ligand-${rank}-min.pdb\
    --out_flex ${outdir}/${system}_flex-${rank}-min.pdb
done