#!/bin/bash

source variables/paths

# Create folds
python ${gscripts}/clustering.py \
    --cpickle  ${clusterdir}/matrix.pickle \
    --input ${typedir}/all.types \
    --output ${typedir}/all