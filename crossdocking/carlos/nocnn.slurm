#!/bin/bash
#SBATCH --array=1-8153
#SBATCH --job cd-carlos-nocnn
#SBATCH --partition=any_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm/slurm-%A-%a.out
#SBATCH --error=slurm/slurm-%A-%a.err

flexdist=3.5
exhaustiveness=8

seed=42

# Root output folder
EXPERIMENT="${SLURM_JOB_NAME}-d${flexdist}"

# Working directory and username
wdir=$(pwd)
user=$(whoami)

# Create output folders
outdir="${wdir}/${EXPERIMENT}"
mkdir -p ${outdir}/info

# Load CUDA
module load cuda/10.2

# Job directory on node
jdir="${user}_${SLURM_JOB_ID}.dcb.private.net"
mkdir /scr/$jdir
cd /scr/$jdir

# Starting time
date | tee $SLURM_ARRAY_TASK_ID.txt
hostname | tee -a $SLURM_ARRAY_TASK_ID.txt
echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}" | tee -a $SLURM_ARRAY_TASK_ID.txt
echo "SLURM_TASKS_PER_NODE ${SLURM_TASKS_PER_NODE}" | tee -a $SLURM_ARRAY_TASK_ID.txt

# Data root
dataroot=/net/pulsar/home/koes/paf46/Research/gnina1.0/
ifile=${dataroot}/carlos_cd_ds_cd.txt
echo "Data root: ${dataroot}" | tee -a $SLURM_ARRAY_TASK_ID.txt
echo "Input file:" ${ifile} | tee -a $SLURM_ARRAY_TASK_ID.txt

isystems=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ifile})
echo "Input systems:" ${isystems} | tee -a $SLURM_ARRAY_TASK_ID.txt

# Ligand and receptor names (with pocket)
ligand=$(echo ${isystems} | cut -f2 -d " ")
receptor=$(echo ${isystems} | cut -f1 -d " ")

echo "Ligand: ${ligand}" | tee -a $SLURM_ARRAY_TASK_ID.txt
echo "Receptor: ${receptor}" | tee -a $SLURM_ARRAY_TASK_ID.txt

# Path to ligand and receptor files
ligpath=${dataroot}/${ligand}
recpath=${dataroot}/${receptor}

# Pocket, ligand and receptor names
pocket=$(dirname ${ligand})
lig=$(basename ${ligand} .sdf)
rec=$(basename ${receptor} .pdb)

echo "Pocket: ${pocket}" | tee -a $SLURM_ARRAY_TASK_ID.txt
echo "Lig: ${lig}" | tee -a $SLURM_ARRAY_TASK_ID.txt
echo "Rec: ${rec}" | tee -a $SLURM_ARRAY_TASK_ID.txt

# Create dir for pocket
mkdir -p ${pocket}

# Copy ligand and receptor files on node
cp ${ligpath} ${ligand}
cp ${recpath} ${receptor}

# Ligand and receptor output names
ligout=${pocket}/flexlig-${lig}-${rec}.sdf
recout=${pocket}/flexrec-${lig}-${rec}.pdb

echo $outdir $pocket
echo $ligout $recout
echo ${outdir}/${pocket}
echo ${pocket}/$(basename ${lig} .sdf)-$(basename ${rec} .pdb)-gnina.log

# Run flexible docking
gnina -r ${receptor} -l ${ligand} \
    --flexdist_ligand ${ligand} --flexdist ${flexdist} \
    --autobox_ligand ${ligand} \
    --exhaustiveness ${exhaustiveness} --num_modes 10 --seed ${seed} \
    --cnn_scoring none \
    --no_gpu --cpu 1 \
    --out ${ligout} --out_flex ${recout} \
    2>&1 | tee ${pocket}/$(basename ${lig} .sdf)-$(basename ${rec} .pdb)-gnina.log

# Compress output
gzip ${ligout}
gzip ${recout}

# Remove input
rm ${ligand} ${receptor}

# Copy docking results on head node
cd ${pocket}
mkdir -p ${outdir}/${pocket}
rsync -a . ${outdir}/${pocket}

# End time
date | tee -a $SLURM_ARRAY_TASK_ID.txt

# Copy info on head node
rsync -a *.txt ${outdir}/info