#!/bin/bash -l
#SBATCH --job-name=predict
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:v100:1
#SBATCH --partition=short
#SBATCH --clusters=htc
#SBATCH --array=0-2
#SBATCH --output=slurm/predict/%x.%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --account=stat-ecr

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
modelname="dense"
prefix="nc"
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
modelname="${modelname}-noaffinity"
if [ ${prefix} = "nc" ]
then
    modelname="${modelname}-nostratified"
fi

subdir=${PWD}
gscripts=${HOME}/git/gninascripts # Script on local filesystem

# GNINA within Singulartiy contianer
container=/data/biggin/lina3015/gsoc19/crossdocking/cd-downsampled/singularity/obabel.sif
groot=${HOME}/git/gnina/
gcaffe=${groot}/caffe/python

trainfolder=training/${modelname}
cd ${trainfolder}

# Predictions of rigid docking model on rigid docking dataset
singularity exec --nv -B ${subdir}:${subdir} --env PYTHONPATH=${gcaffe}:/usr/local/python \
    ${container} \
    python ${gscripts}/predict.py \
    --model ${modelname}.model \
    --weights ${modelname}.*.${SLURM_ARRAY_TASK_ID}_iter_*.caffemodel \
    --data_root ${subdir} \
    --input ${prefix}test${SLURM_ARRAY_TASK_ID}.types \
    --output ${prefix}test${SLURM_ARRAY_TASK_ID}.out \
    --gpu 0

# Predictions of rigid docking model on flexible docking dataset
singularity exec --nv -B ${subdir}:${subdir} -B ${subdir}/../cd-downsampled/:${subdir}/../cd-downsampled/ \
    --env PYTHONPATH=${gcaffe}:/usr/local/python \
    ${container} \
    python ${gscripts}/predict.py \
    --model ${modelname}.model \
    --weights ${modelname}.*.${SLURM_ARRAY_TASK_ID}_iter_*.caffemodel \
    --data_root ${subdir}/../cd-downsampled/ \
    --input ${subdir}/../cd-downsampled/files/${prefix}test${SLURM_ARRAY_TASK_ID}.types \
    --output FLEX${prefix}test${SLURM_ARRAY_TASK_ID}.out \
    --gpu 0