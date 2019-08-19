#!/bin/bash

traindir=$1
caffemodelnum=$2

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
    -c default2017-pose.${caffemodelnum} \
    -d ../${dataroot} \
    -o top