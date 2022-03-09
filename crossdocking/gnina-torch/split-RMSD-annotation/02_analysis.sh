#!/bin/bash

mkdir -p plots

for model in "default2017-nc" "default2018-nc" #"dense-nc"
do
    for prefix in "flex1" #"flex2" "max2"
    do
        python plot.py ${model} ${prefix}
    done
done
