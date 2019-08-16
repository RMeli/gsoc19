#!/bin/bash

traindir=$1

source variables/paths

cd ${traindir}

# scripts.plot.roc is available at https://github.com/RMeli/scripts/
python -m scripts.plot.roc \
    -i *.?.auc.finaltrain *.?.auc.finaltest \
     -g 0 0 0 1 1 1 \
     -l "Train 0" "Train 1" "Train 2" "Test 0" "Test 1" "Test 2" \
     -o auc.pdf