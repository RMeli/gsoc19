# ML Traning

## Running Scripts

The scripts can be run with [Singularity](https://sylabs.io/singularity/) as follows:
```
singularity exec --nv PATH_TO_CONTAINER PATH_TO_SCRIPT
```

## Folds

### ProBiS

The [ProBiS](https://probis.nih.gov/) server can be used to split the systems based on the similarity of the binding site. The file `probis.csv` contains a 3-fold split of the PDBbind17 dataset, provided by David R. Koes and Jocelyn Sunseri.

### Receptor Sequence Similarity