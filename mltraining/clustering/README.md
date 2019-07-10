# Clustering

In order to obtain cross-validation folds, clustering is performed on the flexible docking dataset in order to split systems according to lignd and protein sequence similarity.

This pipeline is based on the following scripts from [gnina/scripts](https://github.com/gnina/scripts):

* `compute_row.py`: Compute ligand and sequence similarity of one system with all the others.
* `combine_rows.py`: Combine the results of `compute_row.py` in a single `pickle` object.
* `clustering.py`: Performs clustering for cross-validation according to the results of the data contained in the `pickle` object.

## Archer Pipeline

The `compute_row.py` takes some time to execute on the PDBbind18 dataset, since it computes ligand similariry and protein sequence similarity of one system with all othe ~16'000 systems.

Since each row is independent of the others, this can be efficiently parallelized on an HPC cluster. The following pipeline uses the [Archer UK National Supercomputing Sevice](http://www.archer.ac.uk/).

### 01 - Maps

Job sumbission on Archer is optimal when asking for 12 nodes (of 24 cores each). This means that with a single job summission, 12 * 24 = 288 instances of `compute_row.py` can be run. The script `01_maps.sh` creates a mapping between a row number and an MPI rank.

### 02 - Archer

The script `02_archer.sub` submits an array job, whose indices correspond to the maps computed above. Each job then uses 288 MPI processes to run `compute_row.py` in parallel, using the Python Task Farm (PTF).

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

The input strung can then be parsed as follows:

```python
args = parser.parse_args(input.split())
```
