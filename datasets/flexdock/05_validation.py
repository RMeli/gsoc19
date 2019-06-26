import MDAnalysis as mda # Needs dev version of MDAnalysis

import os
import re
import warnings

datasets = ["refined"]

pdbbindpath="../../PDBbind18"

def newline(yes : bool) -> str:
    return "\n" if yes else ""

def load(fpath: str, print_warnings=False) -> mda.Universe:
    """
    Load file with MDAnalysis, suppressing warnings by default.

    MDAnslysis complains about MOL2 atom types and some metal atoms.
    """

    if not print_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            return mda.Universe(fpath)

    return mda.Universe(fpath)


for dataset in datasets:

    # List all systems in dataset
    systems = [
        s
        for s in os.listdir(dataset)
        if re.match("^....$", s) and os.path.isdir(os.path.join(dataset, s))
    ]

    for system in systems:
        ok : bool = True
        print(f"Validating {dataset}/{system}...", end="")

        crecpath = os.path.join(pdbbindpath, dataset, system, f"{system}_protein.pdb")
        cligpath = os.path.join(pdbbindpath, dataset, system, f"{system}_ligand.mol2")

        crec = load(crecpath)  
        clig = load(cligpath)

        n_residues = len(crec.residues)
        n_hvy_atoms_rec = len(crec.select_atoms("protein and not type H"))
        n_hvy_atoms_lig = len(clig.select_atoms("not type H"))

        recnames = [
            rec for rec in os.listdir(os.path.join(dataset, system))
            if re.match(system + "_protein-?[0-9]{0,2}\.pdb", rec)
        ]
        for recname in recnames:
            recpath = os.path.join(dataset, system, recname)

            rec = load(recpath)

            try:
                assert len(rec.residues) == n_residues
            except AssertionError:
                print(f"{newline(ok)}\tWrong number of residues!")
                ok = False

            try:
                assert len(rec.select_atoms("protein and not type H")) == n_hvy_atoms_rec
            except AssertionError:
                print(f"{newline(ok)}\tWrong number of receptor heavy atoms!")
                ok = False
            
            if not ok:
                break

        lignames = [
            lig for lig in os.listdir(os.path.join(dataset, system))
            if re.match(system + "_ligand-?[0-9]{0,2}\.pdb", lig)
        ]
        for ligname in lignames:
            ligpath = os.path.join(dataset, system, ligname)

            lig = load(ligpath)

            try:
                assert len(lig.select_atoms("not type H")) == n_hvy_atoms_lig
            except:
                print(f"{newline(ok)}\tWrong number of ligand heavy atoms!")
                ok = False
            
            if not ok:
                break

        flexnames = [
            flex for flex in os.listdir(os.path.join(dataset, system))
            if re.match(system + "_flex-?[0-9]{0,2}\.pdb", flex)
        ]
        for flexname in flexnames:
            flexpath = os.path.join(dataset, system, flexname)

            flex = load(flexpath)

            # Check there are no PRO residues within the flexible residues
            try:
                PRO = flex.select_atoms("protein and resname PRO")
                assert len(PRO) == 0
            except AssertionError:
                print(f"{newline(ok)}\tFlexible PRO residues!")
                ok = False
            
            if not ok:
                break


        if ok:
            print("ok")