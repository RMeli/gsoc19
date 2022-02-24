#!/bin/bash

traindir=$1

# Build full typefiles for testing
# The full typefiles contain all the poses of all systems contained in the test fold
# This allows to test on all poses and compare with docking

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

source variables/paths

cd ${traindir}

for fold in 0 1 2
do
    # Get list of systems in the training set
    cat alltest${fold}.types | awk '{print $2}' | sed "s#/.*##g" | sort -u > test${fold}-systems.lst

    # Construct fulltest.types file by picking systems from all.types in gninatypes
    # all.types in gninatypes contain ALL the systems
    rm -f fulltest${fold}.types
    for system in $(cat test${fold}-systems.lst)
    do
      grep ${system} ../${typedir}/all.types >> fulltest${fold}.types
    done

done