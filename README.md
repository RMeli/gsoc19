# CNN Scoring for Flexible Docking

[![DOI](https://zenodo.org/badge/185625910.svg)](https://zenodo.org/badge/latestdoi/185625910)

[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![Powered by RDKit](https://img.shields.io/static/v1?label=Powered%20by&message=RDKit&color=3838ff&style=flat&logo=data:image/x-icon;base64,AAABAAEAEBAQAAAAAABoAwAAFgAAACgAAAAQAAAAIAAAAAEAGAAAAAAAAAMAABILAAASCwAAAAAAAAAAAADc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/////+B////AP///gB///wAP//4AB//+AAf//gAH//4AB//+AAf//gAH//8AD///gB///8A////gf////////)](https://www.rdkit.org/)

## Abstract

> Molecular docking—the prediction of binding modes and binding affinity of a molecule to a target of known structure—is a great computational tool for structure-based drug design. However, docking scoring functions are mostly empirical or knowledge-based and the flexibility of the receptor is completely neglected in most docking studies. Recent advances in the field showed that scoring functions can be effectively learnt by convolutional neural networks (CNNs). Here we want to build on top of these findings and develop a CNN scoring function for flexible docking by extending the capabilities of `gnina`—a state-of-the-art deep learning framework for molecular docking—and by building an high-quality training dataset for flexible docking.

## Project Description

This Google Summer of Code 2019 project aims to extend the capabilities of [gnina](https://github.com/gnina), the deep learning framework for molecular docking devloped in [David Koes](http://bits.csb.pitt.edu/)'s group, to build a CNN-based scoring function for docking with flexible side chains.

The main stages of the project are the following:

* Build a high-quality training dataset of docking with flexible side chains
  * Get and pre-process [PDBbind18](http://www.pdbbind.org.cn/) (see [`PDBbind18/README.md`](PDBbind18/README.md))
  * Re-docking with flexible side chains using [smina](https://sourceforge.net/projects/smina/) (see [`docking/README.md`](docking/README.md))
  * Optimize crystal poses using [smina](https://sourceforge.net/projects/smina/) (see [`optimisation/README.md`](optimisation/README.md))
  * Build training dataset (see [`datasets/flexdock/README.md`](datasets/flexdock/README.md))
* Enable optimisation of flexible side chains (see [PR #73](https://github.com/gnina/gnina/pull/73))
  * Split ligand and receptor movable atoms in the correct channels
  * Combine ligand and receptor gradients for geometry optimisation
* Train a new CNN-based scoring function for docking with flexible side chains (see [`mltraining/README.md`](mltraining/README.md))
  * Evaluate the performance of pose prediction
  * Evaluate the performance of pose optimisation
* Iterate training on datasets augmented with CNN-optimized poses
  * Optimize docking poses with the CNN (see [`mlopt/README.md`](mlopt/README.md))

This repository collects the different pipelines built in order to achieve the project goals. A list of constributions and fixes to [openbabel](https://github.com/openbabel/openbabel), [smina](https://sourceforge.net/projects/smina/) and [gnina](https://github.com/gnina/gnina) ([OpenChemistry](https://www.openchemistry.org) organisation) and [MDAnalysis](https://github.com/MDAnalysis/mdanalysis) ([NumFocus](https://numfocus.org) organisation) is given below.

The datasets related to this project will be released on [Zenodo](https://zenodo.org/) in due time.

[![Poster](https://img.shields.io/badge/Poster-10.5281%2Fzenodo.3385972-blue)](https://doi.org/10.5281/zenodo.3385972)

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
