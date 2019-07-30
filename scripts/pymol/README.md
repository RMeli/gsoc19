# PyMol

## Installation

### MacOS

Incentive PyMol:

```bash
conda create -n pymol python=3.6 ipython
conda install -n pymol -c schrodinger pymol
```

Open-source PyMol:

```bash
brew insall brewsci/bio/pymol
```

## Scripts

In order to make the functions defined by `pymol.cmd.extend(PYMOL_FUNCTION_NAME, FUNCTION)` available in PyMol, the script containing their definitions has to be run from PyMol's CLI:

```pymol
run PATH_TO_SCRIPT
```

### Load Complex

The script `loadcomplex.py` allows to load a PDBbind protein-ligand complex, given a `PDBID` and a `DATASET` (`refined`, `other` or `test`). The residues within 3Å from the ligand are outlined.

```pymol
loadcomplex PDBID, DATASET
```

### Load FlexDock

The script `loadflexdock.py` allows to load a ligand and flexible residues binding mode, given a `PDBID`, a `DATASET` (`refined`, `other` or `test`) and the binding mode index `BMIDX` (from 1 to the number of modes fount for `PDBID`). The ligand and flexible residues are outlined. The protein structure and the residues within 3Å from the ligand are also shown.

```pymol
loadflexdock PDBID, DATASET, BMIDX
```
