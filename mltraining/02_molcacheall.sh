#!/bin/bash

source variables/paths

echo -n "Building typefile..."

python ${pscripts}/buildtypefile.py \
    ${database} \
    ${typedir} \
    --lmin 2 --lmax 4 --fmin 1 --fmax 1.5 \
    -d refined other \
    --verbose --all \
    -o ${typedir}/all.types \
    | tee ${typedir}/buildtypefile.log

echo "done"
echo -n "Creating caches..."

python ${gscripts}/create_caches2.py \
    "${typedir}/all.types" \
    -c 1 -d ${typedir} \
    --recmolcache "${typedir}/rec.molcache2" \
    --ligmolcache "${typedir}/lig.molcache2"

echo "done"