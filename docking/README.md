## Flexible Docking

Flexible docking is performed using `smina` (Feb 12 2019) on the PDBbind18 dataset.

Flexible docking is performed using automatic flexible residue (`--flexdist_ligand`) and box size (`--autobox_ligand`) detection, with the following parameters (see `variables/docking`):

| Option            | Value | 
| :----------------:|:-----:|
| `--flexdist`      | 3     |
| `--autobox_add`   | 10    |
| `--exhaustiveness`| 8     |
| `--num_modes`     | 20    |


### Serial Docking

The script `docking.sh` performs docking in serial, i.e. going trough one system at a time. Given the large size of the PDBbind18 dataset and the increased computational cost, docking should be run in parallel.


### Parallel Docking

The script `docking.sub` allows to submit different docking jobs on a HPC cluster. The PDBbind18 dataset is split in chunks of 100 systems and each chunk is docked in parallel.

Only a job for the Sun Grid Engine (SGE) scheduler is provided (see `templates/sge.job`), but job scripts for other schedulers can be easily added.

### Analysis

The script `analysis.py` allows to automatically check how many systems have been docked successfully and provides some statistics. In addition, problematic systems where no conformations were found in the given search space are detected and listed.