#!/bin/bash

source variables/paths
source variables/annotation

python ${pscripts}/buildtypefiles.py \
    ${database} \
    ${typedir}

# Concatenate all .types files into one
cat ${typedir}/????/????.types > ${typedir}/all.types