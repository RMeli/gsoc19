#!/bin/bash

source variables/paths

traindir=$1

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

# Create folds
python ${gscripts}/clustering.py \
    --cpickle ${clusterdir}/dist-lsim.pickle \
    --input ${traindir}/all.types \
    -s2 0.4 \
    -v --output ${traindir}/all \
    | tee ${traindir}/clustering.log