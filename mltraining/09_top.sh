#!/bin/bash

traindir=$1
caffemodelnum=$2

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

if [[ $caffemodelnum == "" ]]
then
  echo "CAFFEMODELNUM must be specified."
  exit
fi

# This scrit is designed to run within a Singularity container
# The correct Python interpreter depend on the container

source ${traindir}/training.vars
source variables/paths

cd ${traindir}

# Caffe
export PYTHONPATH=${gcaffe}:${PYTHONPATH}

# Python 3.6 for Singularity container U18-C101
python3.6 ${gscripts}/calctop.py \
    -m ../${modeldir}/${model} \
    -p ${prefix} \
    -c def2017-pose-33.5.${caffemodelnum} \
    -d ../${dataroot} \
    -t 10 \
    -o top