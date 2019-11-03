#!/bin/bash

traindir=$1

if [[ $traindir == "" ]]
then
  echo "TRAINDIR must be specified."
  exit
fi

source variables/paths

cd ${traindir}

mkdir -p analysis

# scripts.plot.roc is available at https://github.com/RMeli/scripts/
python -m scripts.plot.roc \
    *.?.auc.finaltrain *.?.auc.finaltest \
     -g 0 0 0 1 1 1 \
     -l "Train 0" "Train 1" "Train 2" "Test 0" "Test 1" "Test 2" \
     -o analysis/auc.pdf


# scripts.plot.pr is available at https://github.com/RMeli/scripts/
python -m scripts.plot.pr \
    *.?.auc.finaltrain *.?.auc.finaltest \
     -g 0 0 0 1 1 1 \
     -l "Train 0" "Train 1" "Train 2" "Test 0" "Test 1" "Test 2" \
     -o analysis/pr.pdf


for fold in 0 1 2
do 
    echo "cnn_score,smina_score,lig_rmsd,fmax_rmsd,system,smina_rank" > fulltest${fold}.csv
    head -n -1 fulltest${fold}.out  | awk '{print $1","$9","$7","$8","$4}' | 
        sed "s#/...._protein-#,#g" | sed 's#\.gninatypes##g' >> fulltest${fold}.csv
done

python ../dist.py fulltest[0-2].csv -o analysis/dist_b.pdf -g boxen
python ../dist.py fulltest[0-2].csv -o analysis/dist_d.pdf -g distplot