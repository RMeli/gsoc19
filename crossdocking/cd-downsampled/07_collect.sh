#!/bin/bash

# --------------------------------------------------
# Collect all RMSDs and scores into a singl CSV file
# --------------------------------------------------

allscores="carlos_cd/rmsds2.csv"

# CSV Header
#csv_header="pocket,protein,ligand,rank,rmsd,obrmsd,flexrmsd,flexobrmsd,fmaxrmsd,score"
csv_header="pocket,protein,ligand,rank,rmsd,obrmsd,flexrmsd,fmaxrmsd,score"

# Create a file containing all the scores
rm -f ${allscores}
rm -f allscores.tmp # Remove tmp file
touch allscores.tmp # Create tmp file
for pocket in $(ls -d carlos_cd/*/)
do
    # Concatenate all score files
    cat allscores.tmp ${pocket}/PDB_Structures/*.rmsds2 > tmp
    mv tmp allscores.tmp
done
grep -v ${csv_header} allscores.tmp > ${allscores} # Remove all headers
echo ${csv_header} | cat - ${allscores} > allscores.tmp # Add top header
mv allscores.tmp ${allscores}