#!/bin/bash

n_cpus=12

datasets="refined other"

optdir=$1

if [[ $optdir == "" ]]
then
  echo "OPTDIR must be specified."
  exit
fi

source variables/paths

cd ${optdir}

# CSV
csv_header="system,rank,rmsd_lig,rmsd_flex,rmsd_fmax,rmsd_tot,score"

# Export variables for GNU parallel
export csv_header=${csv_header}
export pdbbind=${pdbbind}
export obrms=${obrms}
export pscripts=${pscripts}
export datasets=${datasets}

# Score one system
score(){
    # Input variables from GNU parallel
    dir=$1

    echo ${dir}

    # PDB name
	system=$(basename ${dir})

    # CSV file
    csvfile=${dir}/${system}_score.csv
    echo ${csv_header} > ${csvfile}

    for dataset in ${datasets}
    do
        ligand_crystal=${pdbbind}/${dataset}/${system}/${system}_ligand.mol2
        protein_crystal=${pdbbind}/${dataset}/${system}/${system}_protein.pdb

        if [ -f ${ligand_crystal} ]
        then
            break
        else
            ligand_crystal=""
            protein_crystal=""
        fi
    done

    # File to store residues RMSD
    rrfname=${dir}"/resrmsd.csv"
    echo "rank,res,rmsd" > ${rrfname}
    
    for ligand in $(ls ${dir}/${system}_ligand-*.pdb)
    do
        # Ligand name and rank
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

            echo $pname $cfname
            
            # Single residue RMSD
            rmsd_res=$(${obrms} ${pfname} ${cfname} | awk  '{print $3}')
            
            rmsdresname=$( echo $(basename $pfname) | sed "s#pflex-##g" | sed "s#.pdb##g" )

            rank=$(echo ${rank} | sed "s#-min##g")
            echo "${rank},${rmsdresname},${rmsd_res}" >> ${rrfname}

            if (( $(echo "$rmsd_res > $rmsd_fmax" | bc -l) )); then
                rmsd_fmax=${rmsd_res}
            fi
        done

        # Remove temporary flex files
        #rm ${dir}/pflex*.pdb ${dir}/cflex*.pdb

        # Combined RMSD
        rmsd_tot=$(echo ${rmsd_lig} + ${rmsd_flex} | bc)
            
        # Score (from ligand file)
        score=$(grep "minimizedAffinity" ${ligand} | awk '{print $3}')

        # Strip rank of -min suffix
        rank=$(echo ${rank} | sed "s#-min##g")

        echo "${system},${rank},${rmsd_lig},${rmsd_fmax}"

        info="${system},${rank},${rmsd_lig},${rmsd_flex},${rmsd_fmax},${rmsd_tot},${score}"
        echo ${info} >> ${csvfile}

    done
}

# Export score for GNU parallel
export -f score

dirs=$(ls -d minimized/* | head -n 2)

echo $dirs

parallel -j ${n_cpus} score ::: ${dirs} ::: ${dataset}
