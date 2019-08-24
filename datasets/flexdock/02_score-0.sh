#!/bin/bash

n_cpus=12

source variables/paths

datasets="test"

# Export variables for GNU parallel
export pdbbind=${pdbbind}
export obrms=${obrms}
export pscripts=${pscripts}

# Score one system
score(){
    # Input variables from GNU parallel
    dir=$1
    dataset=$2

    echo ${dir}

    # PDB name
	system=$(basename ${dir})

    # CSV file
    csvfile=${dir}/${system}_score.csv

    ligand_crystal=${pdbbind}/${dataset}/${system}/${system}_ligand.mol2
    protein_crystal=${pdbbind}/${dataset}/${system}/${system}_protein.pdb

    # File to store residues RMSD
    rrfname=${dir}"/resrmsd.csv"
    
    ligand=${dir}/${system}_ligand-0.pdb

    name=$(basename $ligand .pdb)
    rank=$( echo $name | sed "s#${system}_ligand-##g" )

    # Ligand RMSD (with OpenBabel)
    rmsd_lig=$(${obrms} ${ligand_crystal} ${ligand} | awk  '{print $3}')

    # Flexible  residues and protein filenames
    flex=${dir}/${system}_flex-${rank}.pdb
    protein=${dir}/${system}_protein-${rank}.pdb

    # Extract flexible residues from current protein and original crystal structure
    python3.6 ${pscripts}/getflex.py ${flex} ${protein} ${protein_crystal} --dir ${dir}
        
    # Compute RMSD for all flexible residues using OBRMS
    rmsd_flex=$(${obrms} ${dir}/pflex.pdb ${dir}/cflex.pdb | awk  '{print $3}')

    # Compute RMSD for single flexible residues using OBRMS (and store MAX)
    rmsd_fmax=-1
    for pfname in $(ls ${dir}/pflex-*.pdb)
    do
        cfname=$( echo $pfname | sed "s#pflex#cflex#g" )
            
        # Single residue RMSD
        rmsd_res=$(${obrms} ${pfname} ${cfname} | awk  '{print $3}')
            
        rmsdresname=$( echo $(basename $pfname) | sed "s#pflex-##g" | sed "s#.pdb##g" )

        echo "${rank},${rmsdresname},${rmsd_res}" >> ${rrfname}

        if (( $(echo "$rmsd_res > $rmsd_fmax" | bc -l) )); then
            rmsd_fmax=${rmsd_res}
        fi
    done

    # Remove temporary flex files
    rm ${dir}/pflex*.pdb ${dir}/cflex*.pdb

    # Combined RMSD
    rmsd_tot=$(echo ${rmsd_lig} + ${rmsd_flex} | bc)
            
    # Score (from ligand file)
    score=$(grep "minimizedAffinity" ${ligand} | awk '{print $3}')

    info="${system},${rank},${rmsd_lig},${rmsd_flex},${rmsd_fmax},${rmsd_tot},${score}"
    echo ${info} >> ${csvfile}
}

# Export score for GNU parallel
export -f score

for dataset in ${datasets}
do
    dirs=$(ls -d ${dataset}/???? )
    
    parallel -j ${n_cpus} score ::: ${dirs} ::: ${dataset}
done