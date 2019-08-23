#!/bin/bash

traindir=$1
optdir=$2

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

if [[ $optdir == "" ]]
then
  echo "OPTDIR must be specified."
  exit
fi

cd ${optdir}

for fold in 0 1 2
do
  cat ../${traindir}/alltest${fold}.types | awk '{print $2}' | sed "s#..../##g" | \
    sed "s#_protein-# #g" | sed "s#.gninatypes##g" > test${fold}.lst
done