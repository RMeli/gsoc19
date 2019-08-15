# ML Traning

## Running Scripts

The scripts can be run with [Singularity](https://sylabs.io/singularity/) as follows:

```bash
singularity exec --nv PATH_TO_CONTAINER PATH_TO_SCRIPT
```

The `--nv` flag is necessary to use `gnina` with a GPU.

## Training Preparation

### 01 - Run GNINA Typer

PDB files and other formats contain a lot of information that is not needed by `gnina`, which works only with atomic coordinates and `smina` atom types. In order to reduce I/O operations during training, receptor and ligand structures can be converted in `.gninatype` binary files.

The script `01_rungninatyper.sh` (to be run within Singularity) uses the `gninatyper` utility program to convert ligand and receptor PDB files into `.gninatype` binary files.

### 02 - Build Type File

The script `02_buildtypefile` creates a `.types` file via the `buildtypefiles.py` script. A `.types` file contains the training example annotation (or label, 0 for a good pose and 1 for a bad pose) and the paths to the receptor and ligand `.gninatype` files.

### 03 - Create Folds

The script `03_createfolds.sh` uses protein sequence distances and ligand similarities, contained in the `matrix.pickle` binary file (see [clustering/README.md](clustering/README.md) for more details), to split the content of the `all.types` file in different folds for cross-validation.

### 04 - Cache

The script `04_molcache.sh` cache all the the `.gninatype` files in monolithic receptor and ligand cache files for reduced I/O operations.

## Training

### 05 - Train

The script `05_train.sh` (to be run within Singularity) trains the CNN model.
