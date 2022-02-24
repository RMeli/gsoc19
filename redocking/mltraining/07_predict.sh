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

# Test on systems contained in original alltestX.types
# Test on systems 
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

# Test on gull list of systems (in a given fold) including poses discarded for training
# Poses with ligand/flex RMSD in [min,max]/[fmin,fmax] are included here
# See 06_fulltypefiles.sh for more details
# Python 3.6 for Singularity container U18-C101
for fold in 0 1 2
do
    python3.6 ${gscripts}/predict.py \
        -m ../${modeldir}/${model} \
        -d ../${dataroot} \
        -w *.${fold}_iter_*.caffemodel \
        -i fulltest${fold}.types \
        -g ${gpu} \
        -s 0 \
        -o fulltest${fold}.out
done