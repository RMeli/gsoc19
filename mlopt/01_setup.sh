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

mkdir -p ${optdir}

cp ../mltraining/variables/completelig ../mltraining/variables/completerec ${optdir}
cp ${traindir}/*.molcache2 ${optdir}
cp ${traindir}/*.model ${optdir}
cp ${traindir}/*.caffemodel ${optdir}