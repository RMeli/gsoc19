#!/bin/bash

source variables/paths

outdir=$1

if [[ $outdir == "" ]]
then
  echo "OUTDIR must be specified."
  exit
fi

# Create folds
python ${gscripts}/clustering.py \
    --cpickle ${clusterdir}/dist-lsim.pickle \
    --input ${outdir}/all.types \
    -s2 0.4 \
    -v --output ${outdir}/all \
    | tee ${outdir}/clustering.log