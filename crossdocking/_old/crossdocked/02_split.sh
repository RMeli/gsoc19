#!/bin/bash

datafile=data/crossdocking-m2.dat

mkdir -p lists

split ${datafile} -a 3 -d -l 50 lists/cd-m2_