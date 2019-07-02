#!/bin/bash

traindir=training

# This scrit is designed to run within a Singularity container
# The correct Python interpreter is python3.6

source variables/training
source variables/paths

mkdir ${traindir} && cd ${traindir}

cp ../${dataroot}/all*.types .
cp ../variables/complete* .

# Caffe
export PYTHONPATH=${gcaffe}:${PYTHONPATH}

python3.6 ${gscripts}/train.py \
    -m ../${modeldir}/${model} \
    -p ${prefix} \
    -d ../${dataroot} \
    -i ${iters} \
    -t ${testinterval} \
    -g ${gpu}