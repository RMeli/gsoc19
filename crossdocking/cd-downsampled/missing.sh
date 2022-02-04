#/bin/bash

# Get SLURM ID of systems with missing RMSD values
# The string can be used directly to setup SLURM array jobs

# Creates missing.csv
python missing.py

missing=""
while read line; do
  pocket=$(echo ${line} | awk -F "," '{print $1}')
  rec=$(echo ${line} | awk -F "," '{print $2}')
  lig=$(echo ${line} | awk -F "," '{print $3}')

  ln=$(grep -n "${pocket}/PDB_Structures/${rec}_PRO_${lig}_LIG_.*" gnina_cmds_1-4999.txt | awk -F: '{print $1}')
  
  if [ -n "${ln}" ]; then
    missing="${missing},${ln}"
  fi
  
  ln=$(grep -n "${pocket}/PDB_Structures/${rec}_PRO_${lig}_LIG_.*" gnina_cmds_5000-7970.txt | awk -F: '{print $1}')
  
  if [ -n "${ln}" ]; then
    missing="${missing},${ln}"
  fi

done < missing.csv

echo ${missing:1}