#!/bin/bash

root=${PWD}
datadir=NOflexcd
cddat=cd.dat

source bash/paths
source bash/docking

datapath=${root}/${datadir}

mkdir -p ${datapath}

cat ${cddat} | while read line
do
    fnames=($line)

    ligand=${fnames[0]}
    receptor=${fnames[1]}

    ligname=$(basename ${ligand})
    recname=$(basename ${receptor})

    lignamea=(${ligname//_/ }) # Split ligand name into array
    ligpdb=${lignamea[0]} # Ligand PDB code
    ligid=${lignamea[1]} # Ligand PDB id code

    recnamea=(${recname//_/ }) # Split receptor name into array
    recpdb=${recnamea[0]} # Receptor PDB code
    recchain=${recnamea[1]} # Recheptor chain

    #echo $ligname $recname
    #echo $ligpdb $ligid
    #echo $recpdb $recchain

    pocket=$(basename $(dirname ${ligand}))
    
    ddir=${datapath}/${pocket}

    mkdir -p ${ddir}/logs
    
    outname="${recpdb}_${recchain}_rec_${ligpdb}_${ligid}"

    ${smina} -r ${receptor} -l ${ligand} \
        --autobox_ligand ${ligand} --autobox_add ${autobox_add} \
	    --exhaustiveness ${exhaustiveness} --num_modes ${num_modes} --cpu ${cpu} \
        --out ${ddir}/${outname}_lig.pdb --out_flex ${ddir}/${outname}_flex.pdb \
        2>&1 | tee ${ddir}/logs/${outname}.log
done