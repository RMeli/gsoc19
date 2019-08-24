# Optimisation

Crystal poses optimised using `smina` can be included in the training/test dataset in order to ensure that there is at least one correct pose per system. This is particularly helpful when all docking poses have an high RMSD and they are all scored similarly by the CNN scoring function.

## 01 - Optimization

The script `01_opt.sh` performs the optimization of crystal poses using `smina`.

## 02 - Make Flex

The script `02_makeflex.sh` uses `makeflex.py` in order to reconstruct the full receptor with the minimised side chains.

## 03 - Augment

The script `03_augment.sh` augment the flexible docking dataset (in `../datasets/flexdock`) with the optimized crystal poses.
