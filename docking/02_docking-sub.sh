#!/bin/bash

source variables/paths
source variables/misc

# Store working directory
wd=$PWD

# Job directory
jdir=jobs
mkdir -p ${jdir}

# Splits directory
tdir=splits

# Loop over lists of systems (splits, docked by a single job)
for s in $(ls ${tdir})
do
    sname=$(basename ${s})
    jobname=${jdir}/${sname}.job

    cp templates/${sn}.job ${jobname}

    sed -i "s#<SPLITNAME>#${s}#g" ${jobname}
    sed -i "s#<SPLITFILE>#${PWD}/${tdir}/${s}#g" ${jobname}

    ${sub} ${jobname}
 
    sleep 1 
done
