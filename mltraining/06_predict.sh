#!/bin/bash

traindir=$1

# This scrit is designed to run within a Singularity container
# The correct Python interpreter depend on the container 

source ${traindir}/training.vars
source variables/paths

cd ${traindir}

cp ../variables/complete* .

# Caffe
export PYTHONPATH=${gcaffe}:${PYTHONPATH}

# Python 3.6 for Singularity container U18-C101
for fold in 0 1 2
do
    python3.6 ${gscripts}/predict.py \
        -m ../${modeldir}/${model} \
        -d ../${dataroot} \
        -w default2017-pose.*.${fold}_iter_*.caffemodel \
        -i alltest${fold}.types \
        -g ${gpu} \
        -s 0 \
        -o test${fold}.out
done