#!/bin/bash

source variables/paths

outdir=$1

# Copy all.types in training folder
cp ${typedir}/all.types ${outdir}

# Create folds
python ${gscripts}/clustering.py \
    --cpickle ${clusterdir}/dist-lsim.pickle \
    --input ${typedir}/all.types \
    -s2 0.4 \
    -v --output ${outdir}/all