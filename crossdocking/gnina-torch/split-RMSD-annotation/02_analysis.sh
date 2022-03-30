#!/bin/bash

mkdir -p plots

for model in "default2017-nc" "default2018-nc" "dense-nc"
do
    for prefix in "flex05" "flex1" #"flex2" "max2"
    do
        python plot.py ${model} ${prefix}
        python plot.py ${model} ${prefix} --suffix stratified
    done
done
