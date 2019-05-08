#!/bin/bash

# Dowlonad PDBbind databases
# - Refined Protein-Ligand (PL)
# - Other Protein-Ligand (PL)
# Extract databases and rename folders

echo "Downloading PDBbind PL (refined)..."
wget http://www.pdbbind.org.cn/download/pdbbind_v2018_refined.tar.gz

echo -n "Extracting PDBbind PL (refined)..."
tar -xzf pdbbind_v2018_refined.tar.gz
mv refined-set refined
echo "done"

echo "Downloading PDBbind PL (other)..."
wget http://www.pdbbind.org.cn/download/pdbbind_v2018_other_PL.tar.gz

echo -n "Extracting PDBbind PL (other)..."
tar -xzf pdbbind_v2018_other_PL.tar.gz
mv v2018-other-PL other
echo -n "done"
