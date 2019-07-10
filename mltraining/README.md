# ML Traning

## Running Scripts

The scripts can be run with [Singularity](https://sylabs.io/singularity/) as follows:

```bash
singularity exec --nv PATH_TO_CONTAINER PATH_TO_SCRIPT
```

## Training Preparation

### Typing

PDB files and other formats contain a lot of information that is not needed by `gnina`, which works only with atomic coordinates and `smina` atom types. In order to reduce I/O operations during training, receptor and ligand structures can be converted in `.gninatype` binary files.

#### 01 - Run GNINA Typer

The script `01_rungninatyper.sh` (to be run within Singularity) uses the `gninatyper` utility program to convert ligand and receptor PDB files into `.gninatype` binary files.

#### 02 - Build Type File

The script `02_buildtypefile` creates a `.types` file via the `buildtypefiles.py` script. A `.types` file contains the training example annotation (or label, 0 for a good pose and 1 for a bad pose) and the paths to the receptor and ligand `.gninatype` files.

### Cross Validation

K-fold cross validation is employed in to evaluate the model performance. Folds for cross validation are based on protein sequence similarity and ligand similarity. See [clustering/README.md](clustering/README.md) for more details.
