# CNN Scoring for Flexible Docking

## Abstract

> Molecular docking—the prediction of binding modes and binding affinity of a molecule to a target of known structure—is a great computational tool for structure-based drug design. However, docking scoring functions are mostly empirical or knowledge-based and the flexibility of the receptor is completely neglected in most docking studies. Recent advances in the field showed that scoring functions can be effectively learnt by convolutional neural networks (CNNs). Here we want to build on top of these findings and develop a CNN scoring function for flexible docking by extending the capabilities of `gnina`—a state-of-the-art deep learning framework for molecular docking—and by building an high-quality training dataset for flexible docking.

## Contributions

### GNINA

List of contributions to [gnina](https://github.com/gnina/gnina) and [gnina-scripts](https://github.com/gnina/scripts):

* Considerably reduced memory usage of `combine_rows.py` ([PR #30](https://github.com/gnina/scripts/pull/30))
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

## Mentors

* [**Dr. David Ryan Koes**](http://bits.csb.pitt.edu/), Assistant Professor, Department of Computational and Systems Biology, University of Pittsburgh
* [**Jocelyn Sunseri**](http://pitt.edu/~jss97/), Computational Biology Doctoral Candidate, Carnegie Mellon and University of Pittsburgh
