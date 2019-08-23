# CNN Scoring for Flexible Docking

## Abstract

> Molecular docking—the prediction of binding modes and binding affinity of a molecule to a target of known structure—is a great computational tool for structure-based drug design. However, docking scoring functions are mostly empirical or knowledge-based and the flexibility of the receptor is completely neglected in most docking studies. Recent advances in the field showed that scoring functions can be effectively learnt by convolutional neural networks (CNNs). Here we want to build on top of these findings and develop a CNN scoring function for flexible docking by extending the capabilities of `gnina`—a state-of-the-art deep learning framework for molecular docking—and by building an high-quality training dataset for flexible docking.

## Project Description

This Google Summer of Code 2019 project aims to extend the capabilities of [gnina](https://github.com/gnina), the deep learning framework for molecular docking devloped in [David Koes](http://bits.csb.pitt.edu/)'s group, to build a CNN-based scoring function for docking with flexible side chains.

The main goals of the project are the following:

* Build a high-quality training dataset of docking with flexible side chains
  * Re-docking of [PDBbind18](http://www.pdbbind.org.cn/) using [smina](https://sourceforge.net/projects/smina/)
  * Receptor reconstruction with docking side chain poses
* Enable optimisation of flexible side chains (see [PR #73](https://github.com/gnina/gnina/pull/73))
  * Split ligand and receptor movable atoms in the correct channels
  * Combine ligand and receptor gradients for geometry optimisation
* Train a new CNN-based scoring function for docking with flexible side chains
  * Evaluate the performance of pose prediction
  * Evaluate the performance of pose optimisation
* Iterate training on datasets augmented with CNN-optimised poses

This repository collects the different pipelines built in order to achieve the project goals. A list of constributions and fixes to [openbabel](https://github.com/openbabel/openbabel), [smina](https://sourceforge.net/projects/smina/) and [gnina](https://github.com/gnina/gnina) (and [MDAnalysis](https://github.com/MDAnalysis/mdanalysis)) is given below.

The datasets built dusing this project will be released on [Zenodo](https://zenodo.org/) in due time.

### Pipelines

#### PDBbind18

The folder `PDBbind18` contains the pipeline to get the dataset, pre-process it, compute SMILES sequences and list the systems withing a given molecular weight threshold. 

See [`PDBbind18/README.md`](PDBbind18/README.md) for further information.

#### Docking

The folder `docking` contain the pipeline to dock the PDBbind18 dataset (both on a multi-core machine or on HPC clusters), analyse the results and cleaup the dataset (removing failures). 

See [`docking/README.md`](docking/README.md) for further information.

#### Datasets

The folder `datasets/flexdock` contains the pipeline to build the training dataset from flexible docking results. The dataset is built by re-constructing the targets with the flexible side chain poses and computing the RMSD for all the poses. After building the dataset and before computing the RMSD of all the poses, the dataset might be augmented with optimized crystal poses (see [opt](####opt)).  Validation and analysis scripts are also present, in order to check the quality of the dataset. 

See [`datasets/flexdock/README.md`](datasets/flexdock/README.md) for further information.

#### Optimization

#### ML Training

#### ML Opt

## Contributions

### GNINA

List of contributions to [gnina](https://github.com/gnina/gnina) and [gnina-scripts](https://github.com/gnina/scripts):

* Optimisation of flexible side chains ([PR #73](https://github.com/gnina/gnina/pull/73))
* Added option to `pymol_arrows.py` ([PR #31](https://github.com/gnina/scripts/pull/31))
* Low-memory and faster substitute `combine_rows.py` ([PR #30](https://github.com/gnina/scripts/pull/30))
* Attempt to decrease memory usage of `combine_rows.py` ([PR #29](https://github.com/gnina/scripts/pull/29))
* Added serialization of `struct residue` ([PR #74](https://github.com/gnina/gnina/pull/74))
* Small fixes to `gninavis` for gradients ([PR #72](https://github.com/gnina/gnina/pull/72))
* Fixed Python3 `pickle` in clustering pipeline ([PR #26](https://github.com/gnina/scripts/pull/26))
* Added insertion code support to `makeflex.py` ([PR #65](https://github.com/gnina/gnina/pull/65))
* Improved `makeflex.py` script to deal with PDB file without atom types ([PR #64](https://github.com/gnina/gnina/pull/64))
* Added test support for newer versions of Boost ([PR #62](https://github.com/gnina/gnina/pull/62))
* Provided documentation and PDB standardization for `makeflex.py` script ([PR #61](https://github.com/gnina/gnina/pull/61))
* Provided fixes for the `makeflex.py` script ([PR #60](https://github.com/gnina/gnina/pull/60))
* Raised issue about `gnina` parallel compilation without `libmolgrid` installed ([Issue #57](https://github.com/gnina/gnina/issues/57))
* Updated `PDBQTUtilities.cpp` to latest OpenBabel version ([PR #59](https://github.com/gnina/gnina/pull/59))

### LibMolGrid

List of contributions to [libmolgrid](https://github.com/gnina/gnina):

* Fixed issue with unsupported CUDA architecture ([PR #5](https://github.com/gnina/libmolgrid/pull/5))

### SMINA

List of contributions to [smina](https://sourceforge.net/projects/smina/):

* Fixed a problem with proline residues, broken by flexible docking ([MR #3](https://sourceforge.net/p/smina/code/merge-requests/3/))

### OpenBabel

List of contributions to [openbabel](https://github.com/openbabel/openbabel):

* Fixed various problems with PDB and PDBQT insertion codes ([PR #1998](https://github.com/openbabel/openbabel/pull/1998))
* Fixed CMake when compiling without RapidJSON ([PR #1988](https://github.com/openbabel/openbabel/pull/1988))

### MDAnalysis

List of contributions to [MDAnalysis](https://github.com/MDAnalysis/mdanalysis):

* Improved mass guess ([PR #2331](https://github.com/MDAnalysis/mdanalysis/pull/2331))
* Fixed issues with PDB `HEADER` field in `PDBReader` and `PDBWriter` ([PR #2325](https://github.com/MDAnalysis/mdanalysis/pull/2325))
* Allowed MOL2 parser to ignore status bit strings ([PR #2319](https://github.com/MDAnalysis/mdanalysis/pull/2319))

## Mentors

* [**Dr. David Ryan Koes**](http://bits.csb.pitt.edu/), Assistant Professor, Department of Computational and Systems Biology, University of Pittsburgh
* [**Jocelyn Sunseri**](http://pitt.edu/~jss97/), Computational Biology Doctoral Candidate, Carnegie Mellon and University of Pittsburgh
