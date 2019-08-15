#!/bin/bash

source variables/paths

traindir=$1

python ${gscripts}/create_caches2.py \
    "${traindir}/alltrain0.types" \
    "${traindir}/alltrain1.types" \
    "${traindir}/alltrain2.types" \
    "${traindir}/alltest0.types" \
    "${traindir}/alltest1.types" \
    "${traindir}/alltest2.types" \
    -c 1 -d ${typedir} \
    --recmolcache "${traindir}/rec.molcache2" \
    --ligmolcache "${traindir}/lig.molcache2"