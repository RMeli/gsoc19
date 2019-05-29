# PDBbind

[PDBbind Website](http://www.pdbbind.org.cn/)

> The aim of the PDBbind database is to provide a comprehensive collection of the experimentally measured binding affinity data for all types of biomolecular complexes deposited in the Protein Data Bank (PDB). It thus provides an essential linkage between energetic and structural information of these complexes, which is helpful for various computational and statistical studies on molecular recognition occurred in biological systems.


## PDBbind 2018

Version 2018 of PDBbind is based on the contents of [PDB](https://www.rcsb.org/) released by the 1st of January 2018.

### Obtaining the Dataset

The script `get.sh` downloads the two parts (`refined` and `other`) of the PDBbind 2018 database for ligand-protein complexes under `PDBbind18/refined` and `PDBbind18/other`, respectively. This consists in 16'126 protein-ligand complexes, of which 4'463 are refined. 

The `refined` set is 2.5Gb, while the `other` set is 6.5Gb.

### Test

A small `test` dataset (only 6 structures) is extracted from the `refined` set. This test set allows to test different pipelines for docking and training.


## References

* Cheng T.J.; Li X.; Li Y.; Liu Z.H.; Wang R.X. *Comparative assessment of scoring functions on a diverse test set*, J. Chem. Inf. Model., 2009; 49(4); 1079-1093.
* Wang, R.; Fang, X.; Lu, Y.; Yang, C.-Y.; Wang, S. *The PDBbind Database: Methodologies and updates*, J. Med. Chem., 2005; 48(12); 4111-4119.
* Wang, R.; Fang, X.; Lu, Y.; Wang, S. *The PDBbind Database: Collection of Binding Affinities for Protein-Ligand Complexes with Known Three-Dimensional Structures*, J. Med. Chem., 2004; 47(12); 2977-2980.