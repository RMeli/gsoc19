# Cross-Docking (Downsampled)

## Data

Structures from [GNINA 1.0 Paper](http://bits.csb.pitt.edu/files/gnina1.0_paper/) cross-docking downsampled dataset (`crossdocked_ds_data.tar.gz`). 

### Changes

The following proteins have been maunally ameded to remove additional models in the PDB file:

* `FKB1A/1F40_PRO.pdb`
* `MMP13/1FM1_PRO.pdb`

The additional models were causing problems with the `makeflex.py` script, that assumes a single model for the receptor. The script has been amended to check the number of models before running and failing with a clear error message.

## Scripts

### `00_downloads.sh`

Download entire cross-docked dataset used in the GNINA 1.0 paper.

> Andrew T. McNutt, Paul Francoeur, Rishal Aggarwal, Tomohide Masuda, Rocco Meli, Matthew Ragoza, Jocelyn Sunseri & David Ryan Koes, *GNINA 1.0: molecular docking with deep learning*, Journal of Cheminformatics volume 13, Article number: 43 (2021).

### `01_make_gnina_commands.sh`

Create file with GNINA commands used to perform cross-docking with flexible residues.

### `02_docking.slurm`

Perform cross-docking with flexible residues using GNINA. Docking is performed with the AutoDock Vina scoring function. 

Docking is parallelised over each system using SLURM arrays.