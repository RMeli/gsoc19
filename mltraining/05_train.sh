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
python3.6 ${gscripts}/train.py \
    -m ../${modeldir}/${model} \
    -p ${prefix} \
    -d ../${dataroot} \
    -i ${iters} \
    -t ${testinterval} \
    --base_lr ${base_lr} \
    -g ${gpu} \
    --percent_reduced ${percent_reduced} \
    --dynamic --lr_policy fixed