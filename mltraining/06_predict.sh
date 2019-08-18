#!/bin/bash

traindir=$1

# This scrit is designed to run within a Singularity container
# The correct Python interpreter depend on the container 

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

source ${traindir}/training.vars
source variables/paths

cd ${traindir}

# Caffe
export PYTHONPATH=${gcaffe}:${PYTHONPATH}

# Python 3.6 for Singularity container U18-C101
for fold in 0 1 2
do
    python3.6 ${gscripts}/predict.py \
        -m ../${modeldir}/${model} \
        -d ../${dataroot} \
        -w *.${fold}_iter_*.caffemodel \
        -i alltest${fold}.types \
        -g ${gpu} \
        -s 0 \
        -o test${fold}.out
done