# CNN Scoring for Flexible Docking

## Abstract

> Molecular docking—the prediction of binding modes and binding affinity of a molecule to a target of known structure—is a great computational tool for structure-based drug design. However, docking scoring functions are mostly empirical or knowledge-based and the flexibility of the receptor is completely neglected in most docking studies. Recent advances in the field showed that scoring functions can be effectively learnt by convolutional neural networks (CNNs). Here we want to build on top of these findings and develop a CNN scoring function for flexible docking by extending the capabilities of `gnina`—a state-of-the-art deep learning framework for molecular docking—and by building an high-quality training dataset for flexible docking.

## Contributions

### GNINA

List of contributions to [gnina](https://github.com/gnina/gnina):

* Provided fixes for the `makeflex.py` script ([PR #60](https://github.com/gnina/gnina/pull/60))
* Raised issue about `gnina` parallel compilation without `libmolgrid` installed ([Issue #57](https://github.com/gnina/gnina/issues/57))
* Updated `PDBQTUtilities.cpp` to latest OpenBabel version ([PR #59](https://github.com/gnina/gnina/pull/59))

### LibMolGrid

List of contributions to [libmolgrid](https://github.com/gnina/gnina):

* Fixed issue with unsupported CUDA architecture ([PR #5](https://github.com/gnina/libmolgrid/pull/5))

### OpenBabel

List of contributions to [openbabel](https://github.com/openbabel/openbabel):

* Fixed CMake when compiling without RapidJSON ([PR #1988](https://github.com/openbabel/openbabel/pull/1988))


## Mentors

* [**Dr. David Ryan Koes**](http://bits.csb.pitt.edu/), Assistant Professor, Department of Computational and Systems Biology, University of Pittsburgh
* [**Jocelyn Sunseri**](http://pitt.edu/~jss97/), Computational Biology Doctoral Candidate, Carnegie Mellon and University of Pittsburgh