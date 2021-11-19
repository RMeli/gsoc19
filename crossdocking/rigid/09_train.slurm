#!/bin/bash -l
#SBATCH --job-name=train
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:v100:1
#SBATCH --partition=medium
#SBATCH --clusters=htc
#SBATCH --array=0-2
#SBATCH --output=slurm/train/%x.%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --account=stat-ecr

hostname
nvidia-smi

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
modelname="default2018"
prefix="nc"
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
modelname="${modelname}-noaffinity"
if [ ${prefix} = "nc" ]
then
    modelname="${modelname}-nostratified"
fi

subdir=${PWD}
gscripts=${HOME}/git/gninascripts/ # Script on local filesystem

# GNINA within Singulartiy contianer
container=/data/biggin/lina3015/gsoc19/crossdocking/cd-downsampled/singularity/obabel.sif
groot=${HOME}/git/gnina/
gcaffe=${groot}/caffe/python

trainfolder=training/${modelname}
mkdir -p ${trainfolder}

cp files/${prefix}train* files/${prefix}test* ${trainfolder}
cp ../cd-downsampled/models/${modelname}.model ${trainfolder}
cp ../cd-downsampled/models/completerec ../cd-downsampled/models/completelig ${trainfolder}
cd ${trainfolder}

echo ${trainfolder}
echo ${subdir}
echo ${prefix}

singularity exec --nv -B ${subdir}:${subdir} --env PYTHONPATH=${gcaffe}:/usr/local/python \
    ${container} \
    python ${gscripts}/train.py \
    --model ${modelname}.model \
    --prefix ${prefix} \
    --data_root ${subdir} \
    --iterations 100000 \
    --test_interval 1000 \
    --foldnums ${SLURM_ARRAY_TASK_ID} \
    --dynamic --lr_policy fixed \
    --gpu 0
