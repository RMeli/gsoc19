#!/bin/bash

source variables/paths

python ${gscripts}/create_caches2.py \
    "${typedir}/alltrain0.types" \
    "${typedir}/alltrain1.types" \
    "${typedir}/alltrain2.types" \
    "${typedir}/alltest0.types" \
    "${typedir}/alltest1.types" \
    "${typedir}/alltest2.types" \
    -c 1 -d ${typedir} \
    --recmolcache "${typedir}/rec.molcache2" \
    --ligmolcache "${typedir}/lig.molcache2"