#!/bin/bash

traindir=training

# This scrit is designed to run within a Singularity container
# The correct Python interpreter depend on the container 

source variables/training
source variables/paths

mkdir ${traindir} && cd ${traindir}

cp ../${dataroot}/all*.types .
cp ../variables/complete* .

# Caffe
export PYTHONPATH=${gcaffe}:${PYTHONPATH}

# Python 3.6 for Singularity container U18-C101
python3.6 ${gscripts}/train.py \
    -m ../${modeldir}/${model} \
    -p ${prefix} \
    -d ../${dataroot} \
    -i ${iters} \
    -t ${testinterval} \
    -g ${gpu}