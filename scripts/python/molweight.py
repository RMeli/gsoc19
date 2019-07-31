import rdkit.Chem as Chem
import rdkit.Chem.Descriptors as Descriptors

import argparse as ap

from typing import Optional

def molweight(fname: str):

    mols = Chem.SDMolSupplier(fname)

    mol = next(mols)

    print(mol)

    return Descriptors.ExactMolWt(mol)

def hvymolweight(fname: str):

    mol = Chem.MolFromSDFFile(fname)

    print(mol)

    return Descriptors.HeavyAtomMolWt(mol)

if __name__ == "__main__":

    def parse(args: Optional[str] = None) -> ap.Namespace:

        parser = ap.ArgumentParser(
            description="Compute RMSD distributions for docking poses."
        )

        parser.add_argument("input", type=str, help="SMILES file (.smi)")

        return parser.parse_args(args)

    args = parse()

    # Get molecular weight from SMILES
    mw = molweight(args.input)

    hamw = hvymolweight(args.input)

    # Print molecular weight
    print(f"{mw:.5f} {hamw:.5f}")