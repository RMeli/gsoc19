#!/bin/bash

source variables/paths
source variables/misc

# Store working directory
wd=$PWD

# Job directory
jdir=jobs
mkdir -p ${jdir}

jobname=${jdir}/docking.job

cp templates/${sn}.job ${jobname}

${sub} ${jobname}
 
sleep 1