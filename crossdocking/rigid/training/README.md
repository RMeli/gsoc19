# Training on Rigid Docking

Training of a CNN model on structures generated with rigid docking.

Only the `ligand` annotation (based on ligand RMSD) makes sense for rigid docking.

The subfolder `ligand` is however included here in order to have the same structure used for flexible docking (where there are different possible annotations, including flexible residues), which allows to use the same scripts.