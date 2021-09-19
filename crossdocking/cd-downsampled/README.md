# Cross-Docking (Downsampled)

## Data

Structures from [GNINA 1.0 Paper](http://bits.csb.pitt.edu/files/gnina1.0_paper/) cross-docking downsampled dataset (`crossdocked_ds_data.tar.gz`). 

### Changes

The following proteins have been maunally ameded to remove additional models in the PDB file:

* `FKB1A/1F40_PRO.pdb`
* `MMP13/1FM1_PRO.pdb`

The additional models were causing problems with the `makeflex.py` script, that assumes a single model for the receptor. The script has been amended to check the number of models before running and failing with a clear error message.
