#!/bin/bash

box_size=28.5

source variables/paths
source variables/annotation

python ${pscripts}/buildtypefiles.py \
    ${database} \
    ${typedir} \
    --lmin ${min} --lmax ${max} --fmin ${fmin} --lmax ${fmax} \
    -L ${box_size}

# Concatenate all .types files into one
find gninatypes -name "????.types" | xargs cat > ${typedir}/all.types