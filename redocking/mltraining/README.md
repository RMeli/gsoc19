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

### 02 - Molcache All

The script `02_molcacheall.sh` cache all the the `.gninatype` files in monolithic receptor and ligand cache files in order to firther reduced I/O operations.

Here we cache all the `.gninatype` files, even the ones discarded used for training. In order to achieve this, a general `all.types` file containing all the files available is created beforehand using the `buildtypefiles.py` script.

Note: `.molcache2` files can be a superset or even a subset of the files described in the `.types` file used for training. In case where the `.molcache2` is a subset, the missing `.gninatype` files are loaded from disk.

### 03 - Build Type File

The script `03_buildtypefile` creates a `.types` file via the `buildtypefiles.py` script. A `.types` file contains the training example annotation (or label, 0 for a good pose and 1 for a bad pose) and the paths to the receptor and ligand `.gninatype` files.

### 04 - Create Folds

The script `04_createfolds.sh` uses protein sequence distances and ligand similarities, contained in the `matrix.pickle` binary file (see [clustering/README.md](clustering/README.md) for more details), to split the content of the `all.types` file in different folds for cross-validation.

## Training

### 05 - Train

The script `05_train.sh` (to be run within Singularity) trains the CNN model.

### Testing

### 06 - Full Typefiles

The `.types` files used for training do not contain all the systems, since some ambiguous poses are discarded for training based on the RMSD annotation. The script `06_fulltypefiles.sh` build `.typefiles` containing all the systems, to be used to re-score all `smina` poses with the trained CNN scoring function. Cross-validation folds are conserved in this operation.

### 07 - Predict

The script `06_predict.sh` score all `smina` poses using the trained CNN scoring function, according to the cross-validation folds previously obtained.

### 08 - Analysis

The script `08_analysis.sh` provide some comparison analysis betewwn the CNN scoring function and `smina` scoring function.
