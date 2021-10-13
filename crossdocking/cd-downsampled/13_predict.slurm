#!/bin/bash -l
#SBATCH --job-name=train
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --partition=short
#SBATCH --clusters=htc
#SBATCH --array=0-2
#SBATCH --output=slurm/predict/%x.%A_%a.out
#SBATCH --time=06:00:00
#SBATCH --account=stat-ecr

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
modelname="default2017-noaffinity"
prefix="cluster"
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

subdir=${PWD}
gscripts=${HOME}/git/gninascripts # Script on local filesystem

# GNINA within Singulartiy contianer
container=/data/biggin/lina3015/gsoc19/crossdocking/cd-downsampled/singularity/obabel.sif
groot=${HOME}/git/gnina/
gcaffe=${groot}/caffe/python

trainfolder=training/${modelname}
mkdir -p ${trainfolder}


cp files/${prefix}train* files/${prefix}test* ${trainfolder}
cp models/${modelname}.model ${trainfolder}
cp models/completerec models/completelig ${trainfolder}
cd ${trainfolder}

singularity exec --nv -B ${subdir}:${subdir} --env PYTHONPATH=${gcaffe}:/usr/local/python \
    ${container} \
    python ${gscripts}/predict.py \
    --model ${modelname}.model \
    --weights ${modelname}.*.${SLURM_ARRAY_TASK_ID}_iter_*.caffemodel \
    --data_root ${subdir} \
    --input ${prefix}test${SLURM_ARRAY_TASK_ID}.types \
    --output ${prefix}test${SLURM_ARRAY_TASK_ID}.out \
    -s 42 \
    --gpu 0
