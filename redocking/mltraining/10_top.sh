#!/bin/bash

traindir=$1

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

source variables/paths

cd ${traindir}

python ../${pscripts}/top.py \
    fulltest -o analysis/top.pdf

python ../${pscripts}/top.py \
    fulltest --ligonly -o analysis/top_ligonly.pdf