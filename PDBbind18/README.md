# PDBbind

[PDBbind Website](http://www.pdbbind.org.cn/)

> The aim of the PDBbind database is to provide a comprehensive collection of the experimentally measured binding affinity data for all types of biomolecular complexes deposited in the Protein Data Bank (PDB). It thus provides an essential linkage between energetic and structural information of these complexes, which is helpful for various computational and statistical studies on molecular recognition occurred in biological systems.

## PDBbind 2018

Version 2018 of PDBbind is based on the contents of [PDB](https://www.rcsb.org/) released by the 1st of January 2018.

### Obtaining the Dataset

The script `01_get.sh` downloads the two parts (`refined` and `other`) of the PDBbind 2018 database for ligand-protein complexes under `PDBbind18/refined` and `PDBbind18/other`, respectively. This consists in 16'126 protein-ligand complexes, of which 4'463 are refined.

The `refined` set is 2.5Gb, while the `other` set is 6.5Gb.

### Test

A small `test` dataset (only 6 structures) is extracted from the `refined` set. This test set allows to test different pipelines for docking and training.

## Pipeline

### 01 - Get

The script `01_get.sh` downloads the two parts (`refined` and `other`) of the PDBbind 2018 database for ligand-protein complexes under `PDBbind18/refined` and `PDBbind18/other`, respectively.

### 02 - SMILES

The script `02_smiles.sh` compute SMILES strings for all the ligands in the PDBbind18 dataset. This script uses [OpenBabel](http://openbabel.org/wiki/Main_Page) to perform the conversion from `.mol2` files to SMILES strings (saved on `.smi` files).

### 03 - List

The script `02_slist.sh` produces three lists:

* A list with all the PDB codes in the PDBbind18 dataset, prepended by the part dataset belong to (`refined` or `other`)
* A list where ligands with heavy-atoms molecular weight higher than 1000 are filtered out (using the `molweight.py` script)
* A list pof PDB codes, paths to the receptor PDB file and paths to the ligand SMILES file (separated by one space)

The first two lists are useful to create splits for docking multiple systems in parallel on HPC clusters.

The last list is used to generate cross-validation folds using `gnina`'s `clustering.py` script (or pipeline).

### 04 - NoWat

The script `04_nowat.sh` renames the original protein PDB files (as `SYSTEM_protein-wat.pdb`) and creates a PDB file with the original names where water molecules are stripped out.

## References

* Cheng T.J.; Li X.; Li Y.; Liu Z.H.; Wang R.X. *Comparative assessment of scoring functions on a diverse test set*, J. Chem. Inf. Model., 2009; 49(4); 1079-1093.
* Wang, R.; Fang, X.; Lu, Y.; Yang, C.-Y.; Wang, S. *The PDBbind Database: Methodologies and updates*, J. Med. Chem., 2005; 48(12); 4111-4119.
* Wang, R.; Fang, X.; Lu, Y.; Wang, S. *The PDBbind Database: Collection of Binding Affinities for Protein-Ligand Complexes with Known Three-Dimensional Structures*, J. Med. Chem., 2004; 47(12); 2977-2980.
