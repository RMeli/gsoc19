#!/bin/bash

n_cpus=12

datasets="refined"

allscores="analysis/allscores.csv"

# CSV Header
csv_header="system,rank,rmsd_lig,rmsd_flex,rmsd_tot,score"

# Create a file containing all the scores
rm -f ${allscores}
rm -f allscores.tmp # Remove tmp file
touch allscores.tmp # Create tmp file
for dataset in ${datasets}
do
    # Concatenate all score files
    cat allscores.tmp ${dataset}/????/????_score.csv > tmp
    mv tmp allscores.tmp
done
grep -v ${csv_header} allscores.tmp > ${allscores} # Remove all headers
echo ${csv_header} | cat - ${allscores} > allscores.tmp # Add top header
mv allscores.tmp ${allscores}
rm -f allscores.tmp # Remove tmp file