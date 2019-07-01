#!/bin/bash

source variables/paths
source variables/annotation

python3.6 ${pscripts}/buildtypefiles.py \
    ${database} \
    ${typedir}