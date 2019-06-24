#!/bin/bash

traindir=training

# This scrit is designed to run within a Singularity container
# The correct Python interpreter is python3.6

source variables/training
source variables/paths

mkdir ${traindir} && cd ${traindir}

cp ../${dataroot}/*.types .
cp ../${dataroot}/alltrain0.types alltest0.types
cp ../${dataroot}/alltrain1.types alltest1.types
cp ../${dataroot}/alltrain2.types alltest2.types
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