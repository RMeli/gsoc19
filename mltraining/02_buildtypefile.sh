#!/bin/bash

box_size=28.5

source variables/paths
source variables/annotation

datasets="redined other"

traindir=$1

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

python ${pscripts}/buildtypefile.py \
    ${database} \
    ${typedir} \
    --lmin ${min} --lmax ${max} --fmin ${fmin} --lmax ${fmax} \
    -L ${box_size} -d ${datasets} \
    -o ${traindir}