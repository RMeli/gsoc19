#!/bin/bash

source variables/paths
source variables/annotation

python ${pscripts}/buildtypefiles.py \
    ${database} \
    ${typedir}

# Concatenate all .types files into one
find gninatypes -name "????.types" | xargs cat > ${typedir}/all.types