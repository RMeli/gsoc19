#!/bin/bash

source variables/paths
source variables/annotation

num_good=0 # Number of good poses
num_bad=0 # Number of bad poses

for dataset in "test"
do
    for system in $(ls -d ${typedir}/????)
    do
        sname=$(basename ${system})

        fscore=${database}/${dataset}/${sname}/${sname}_score.csv # Score file (w/ RMSD)

        fgtypes=${typedir}/${sname}.types # gnina types file

        python ${pscripts}/buildtypefile.py \
            --output ${fgtypes} \
            --min ${min} --max ${max} \
            ${fscore} ${typedir}

        good=$(grep "^0" ${fgtypes} | wc | awk '{print $1}')
        bad=$(grep "^1" ${fgtypes} | wc | awk '{print $1}')

        num_good=$(expr ${num_good} + ${good})
        num_bad=$(expr ${num_bad} + ${bad})
    done
done

echo "Number of good poses: ${num_good}"
echo "Number of bad poses: ${num_bad}"