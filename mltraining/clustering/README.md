# Clustering

In order to obtain cross-validation folds, clustering is performed on the flexible docking dataset in order to split systems according to ligand similarity and protein sequence distance.

This pipeline is based on the following scripts from [gnina/scripts](https://github.com/gnina/scripts):

* `compute_row.py`: Compute ligand and sequence similarity of one system with all the others.
* `combine_rows.py`: Combine the results of `compute_row.py` in a single `pickle` object.
* `clustering.py`: Performs clustering for cross-validation according to the results of the data contained in the `pickle` object.

The `combine_rows.py` script was running out of memory on a machine with 64Gb of RAM, therefore it has ben rewritten as `combine_rows_lowmem.py`.

## SMILES and Protein Sequences

### 01 - Sequences

The script `01_seqs.sh` allows to extract the protein sequences from the PDBbind18 dataset and store them in the file `seqs.dat`, together with the corresponding ligand SMILES.

The format of the `seqs.dat` file is the following:
```
PDBCODE SMILES SEQUENCE
```

## Archer Pipeline

The `compute_row.py` takes some time to execute on the PDBbind18 dataset, since it computes ligand similarity and protein sequence distance of one system with all other ~16'000 systems.

Since each row is independent of the others, this can be efficiently parallelized on an HPC cluster. The following pipeline uses the [Archer UK National Supercomputing Sevice](http://www.archer.ac.uk/).

### 02 - Maps

Job sumbission on Archer is optimal when asking for 12 nodes (of 24 cores each). This means that with a single optimal job submission, 12 * 24 = 288 instances of `compute_row.py` can be run. The script `02_maps.sh` creates a mapping between a row number and one of the 288 MPI ranks.

### 03 - Compute Row

The script `03_computerow.sub` submits an array job, whose indices correspond to the maps computed above. Each job then uses 288 MPI processes to run `compute_row.py` in parallel, using the Python Task Farm (PTF).

The PTF passes the MPI rank as last argument of the python program. The beginning of the main function of the original `compute_row.py` script has to be modified as follow in order to read a map and the MPI rank from input:

```python
import sys
import numpy as np

pdbseqs=sys.argv[1]
map=sys.argv[2]
id=int(sys.argv[3])

maparray=np.loadtxt(map,dtype=int)

if maparray[id,0] == id:
    row = maparray[id,1]
else:
    raise Exception

input="--pdbseqs {pdbseqs} --row {row} --out rows/row-{row}".format(row=row, pdbseqs=pdbseqs)
```

The new custom input string can then be parsed as follows:

```python
args = parser.parse_args(input.split())
```

## Post-Processing

### 04 - Combine Rows

The scropt `04_combinerows.sh` allows to combine all the rows computed with `03_computerow.sub` in a single `matrix.pickle` binary file, containing protein sequence distances and ligand similarities.

This script calls `combine_rows_lowmem.py`, a lowe-memory version of the original `combine_rows.py` script (see [gnina/scripts](https://github.com/gnina/scripts)).
