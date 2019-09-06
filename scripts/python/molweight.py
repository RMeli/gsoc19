import rdkit.Chem as Chem
import rdkit.Chem.Descriptors as Descriptors

import argparse as ap

from typing import Optional


def smiles(fname: str):

    with open(fname, "r") as f:
        smiles = f.read().strip()

    return smiles


def molweight(smiles: str):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Could not build molecule from smiles: {smiles}")
        exit(1)

    return Descriptors.ExactMolWt(mol)


def hvymolweight(smiles: str):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Could not build molecule from smiles: {smiles}")
        exit(1)

    return Descriptors.HeavyAtomMolWt(mol)


if __name__ == "__main__":

    def parse(args: Optional[str] = None) -> ap.Namespace:

        parser = ap.ArgumentParser(
            description="Compute molecular weight given a SMILES file."
        )

        parser.add_argument("input", type=str, help="SMILES file (.smi)")

        return parser.parse_args(args)

    args = parse()

    # Get SMILES from file
    smi = smiles(args.input)

    # Get molecular weight from SMILES
    mw = molweight(smi)
    hamw = hvymolweight(smi)

    # Print molecular weight
    print(f"{mw:.5f} {hamw:.5f}")
