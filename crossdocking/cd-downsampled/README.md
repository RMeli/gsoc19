# Cross-Docking (Downsampled)

## Data

Structures from [GNINA 1.0 Paper](http://bits.csb.pitt.edu/files/gnina1.0_paper/) cross-docking downsampled dataset (`crossdocked_ds_data.tar.gz`). 

### Changes

* Protein `FKB1A/1F40_PRO.pdb` has been manually modified to remove the additional (spuriout models) in the PDB 
    * Additional models were causing problems with `makeflex.py`
    * This problem was spotted with the `validation.py` script