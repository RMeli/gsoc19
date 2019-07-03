#!/bin/bash

n_cpus=2

source variables/paths

datasets="test"

# File with all the scores
allscores="analysis/allscores.csv"

# CSV Header
csv_header="system,rank,rmsd_lig,rmsd_flex,rmsd_tot,score"

# Export variables for GNU parallel
export csv_header=${csv_header}
export pdbbind=${pdbbind}
export obrms=${obrms}
export pscripts=${pscripts}

# Score one system
score(){
    # Input variables from GNU parallel
    dir=$1
    dataset=$2

    echo ${dir} ${dataset}

    # PDB name
	system=$(basename ${dir})

    # CSV file
    csvfile=${dir}/${system}_score.csv
    rm -f ${csvfile}

    echo ${csv_header} > ${csvfile}

    ligand_crystal=${pdbbind}/${dataset}/${system}/${system}_ligand.mol2
    protein_crystal=${pdbbind}/${dataset}/${system}/${system}_protein.pdb

    for ligand in $(ls ${dir}/${system}_ligand-*.pdb)
    do
        # Ligand name and rank
        name=$(basename $ligand .pdb)
        rank=$( echo $name | sed "s#${system}_ligand-##g" )

        # Ligand RMSD (with OpenBabel)
        rmsd_lig=$(${obrms} ${ligand_crystal} ${ligand} | awk  '{print $2}')

        # Flexible residues RMSD (with MDAnalysis)
        flex=${dir}/${system}_flex-${rank}.pdb
        protein=${dir}/${system}_protein-${rank}.pdb
        rmsd_flex=$(python3.6 ${pscripts}/flexrmsd.py ${flex} ${protein} ${protein_crystal})

        # Combined RMSD
        rmsd_tot=$(echo ${rmsd_lig} + ${rmsd_flex} | bc)
            
        # Score (from ligand file)
        score=$(grep "minimizedAffinity" ${ligand} | awk '{print $3}')

        info="${system},${rank},${rmsd_lig},${rmsd_flex},${rmsd_tot},${score}"
        echo ${info} >> ${csvfile}
       
    done
}

# Export score for GNU parallel
export -f score

for dataset in ${datasets}
do
    dirs=$(ls -d ${dataset}/????)

    parallel -j ${n_cpus} score ::: ${dirs} ::: ${dataset}

done

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