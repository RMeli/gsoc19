#!/bin/bash

optdir=$1

# Datasets
# ! Specify as list with newline since IFS is modified later !
datasets="refined
other"

if [[ $optdir == "" ]]
then
  echo "OPTDIR must be specified."
  exit
fi

cd ${optdir}

source ../variables/paths

for fold in 0
do
# IFS='\n' is needed in order to loop over one line
# A line contains the PDB code of the system and its rank
IFS='
'
  for line in $(cat *${fold}.lst)
  do
    system=$(echo $line | awk '{print $1}')
    rank=$(echo $line | awk '{print $2}')

    # Check in which dataset (refined or other) the system is
    for dataset in ${datasets}
    do
      # Check if file exists
      if [ -f ${dataroot}/${dataset}/${system}/${system}_flex.info ] 
      then
        # If file exist the dataset is found and stored
        dir=${dataset}/${system}
        break # No need to look further
      fi
    done

    outdir=minimized/${system}
    mkdir -p ${outdir}

    ligand=${dataroot}/${dir}/${system}_ligand-${rank}.pdb
    receptor=${dataroot}/${dir}/${system}_protein-${rank}.pdb
    flex=${dataroot}/${dir}/${system}_flex-${rank}.pdb

    echo "ligand = ${ligand}"
    echo "receptor = ${receptor}"

    # TODO: Remove this copy
    #cp ${ligand} ${outdir}
    #cp ${flex} ${outdir}

    ${gnina} -l ${ligand} -r ${receptor} \
    --flexres $(cat ${dataroot}/${dir}/${system}_flex.info) \
    --cnn_model *.model \
    --cnn_weights *.${fold}_iter_*.caffemodel \
    --cnn_scoring --minimize  \
    --gpu \
    --out ${outdir}/${system}_ligand-${rank}-min.pdb\
    --out_flex ${outdir}/${system}_flex-${rank}-min.pdb
  done
done