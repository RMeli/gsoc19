#!/bin/bash

# Only for Linux
wget http://insilab.org/files/probis-algorithm/probis
chmod u+x probis

# Download full Carlos dataset
wget http://disco.csb.pitt.edu/disco.tar.gz
mkdir -p disco && cd disco
mv ../disco.tar.gz .
tar -xvf disco.tar.gz
