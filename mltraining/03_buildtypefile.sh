#!/bin/bash

source variables/paths
source variables/annotation

datasets="refined other"

traindir=$1
box_size=$2

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

if [[ $box_size == "" ]]
then
  echo "BOXSIZE must be specified."
  exit
fi

python ${pscripts}/buildtypefile.py \
    ${database} \
    ${typedir} \
    --lmin ${min} --lmax ${max} --fmin ${fmin} --fmax ${fmax} \
    -L ${box_size} -d ${datasets} \
    --verbose -o ${traindir}/all.types \
    | tee ${traindir}/buildtypefile.log