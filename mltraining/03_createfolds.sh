#!/bin/bash

source variables/paths

# Create folds
python ${gscripts}/clustering.py \
    --cpickle ${clusterdir}/matrix.pickle \
    --input ${typedir}/all.types \
    -s2 0.40 \
    --output ${typedir}/all