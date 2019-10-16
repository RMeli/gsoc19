#!/bin/bash

source variables/paths

echo -n "Building typefile..."

# Build type file without checking box size
# This is useful to caching all systems (the cache can be a superset of the ones used)
# See README.md for more infos
python ${pscripts}/buildtypefile.py \
    ${database} \
    ${typedir} \
    --lmin ${min} --lmax ${max} --fmin ${fmin} --fmax ${fmax} \
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