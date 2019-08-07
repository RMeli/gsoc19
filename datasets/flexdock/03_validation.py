import MDAnalysis as mda # Needs dev version of MDAnalysis

import numpy as np
import pandas as pd

import os
import re
import warnings

datasets = ["refined", "other"]

pdbbindpath="../../PDBbind18"

residue_selection = "protein and not (type H or name H*)"
ligand_selection = "not type H and not (type H or name H*)"
water_selection = "not protein and (resname WAT or resname HOH)"


def newline(yes : bool) -> str:
    return "\n" if yes else ""


def load(fpath: str, print_warnings=False) -> mda.Universe:
    """
    Load file with MDAnalysis, suppressing warnings by default.

    MDAnslysis complains about MOL2 atom types and some metal atoms in PDB files.
    """

    if not print_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            return mda.Universe(fpath)

    return mda.Universe(fpath)

with open("analysis/invalid.lst", "w", buffering=1) as finvalid, open("analysis/valid.lst", "w", buffering=1) as fvalid:

    # Loop over datasets
    for dataset in datasets:

        # List all systems in dataset
        systems = [
            s
            for s in os.listdir(dataset)
            if re.match("^....$", s) and os.path.isdir(os.path.join(dataset, s))
        ]

        # Loop over systems
        for system in systems:
            ok : bool = True
            print(f"Validating {dataset}/{system}...", end="")

            crecpath = os.path.join(pdbbindpath, dataset, system, f"{system}_protein.pdb")
            cligpath = os.path.join(pdbbindpath, dataset, system, f"{system}_ligand.mol2")

            crec = load(crecpath)
            clig = load(cligpath)

            n_residues = len(crec.residues)
            assert n_residues != 0

            n_hvy_atoms_rec = len(crec.select_atoms(residue_selection))
            assert n_hvy_atoms_rec != 0
            
            n_hvy_atoms_lig = len(clig.select_atoms(ligand_selection))
            assert n_hvy_atoms_lig != 0
            
            n_water = len(crec.select_atoms(water_selection))

            recnames = [
                rec for rec in os.listdir(os.path.join(dataset, system))
                if re.match(system + "_protein-[0-9]{0,2}\.pdb", rec)
            ]
            for recname in recnames:
                recpath = os.path.join(dataset, system, recname)

                try:
                    rec = load(recpath)
                except Exception:
                    print(f"{newline(ok)}    Failed loading {recpath}")
                    ok = False
                    continue

                try:
                    sel = rec.select_atoms(residue_selection)
                    assert len(sel) == n_hvy_atoms_rec


                    try:
                        assert len(rec.residues) == n_residues
                    except AssertionError:
                        print(f"{newline(ok)}    Wrong number of residues ({len(rec.residues)} vs {n_residues})!")
                        ok = False

                except AssertionError:
                    print(f"{newline(ok)}    Wrong number of receptor heavy atoms ({len(sel)} vs {n_hvy_atoms_rec})!")
                    ok = False

                try:
                    water = rec.select_atoms(water_selection)
                    assert len(water) == n_water
                except AssertionError:
                    print(f"{newline(ok)}    Wrong number of water molecules ({len(water)} vs {n_water})!")
                    ok = False

                
                if not ok:
                    break

            lignames = [
                lig for lig in os.listdir(os.path.join(dataset, system))
                if re.match(system + "_ligand-[0-9]{1,2}\.pdb", lig)
            ]
            for ligname in lignames:
                ligpath = os.path.join(dataset, system, ligname)

                try:
                    lig = load(ligpath)
                except Exception:
                    print(f"{newline(ok)}    Failed loading {recpath}")
                    ok = False
                    continue

                try:
                    sel = lig.select_atoms(ligand_selection)
                    assert len(sel) == n_hvy_atoms_lig
                except:
                    print(f"{newline(ok)}    Wrong number of ligand heavy atoms ({len(sel)} vs {n_hvy_atoms_lig})!")
                    ok = False
                
                if not ok:
                    break

            flexnames = [
                flex for flex in os.listdir(os.path.join(dataset, system))
                if re.match(system + "_flex-[0-9]{1,2}\.pdb", flex)
            ]
            for flexname in flexnames:
                flexpath = os.path.join(dataset, system, flexname)

                flex = load(flexpath)

                # Check there are no PRO residues within the flexible residues
                try:
                    PRO = flex.select_atoms("protein and resname PRO")
                    assert len(PRO) == 0
                except AssertionError:
                    print(f"{newline(ok)}    Flexible PRO residues!")
                    ok = False
                
                if not ok:
                    break

            # Load score
            scorepath = os.path.join(dataset, system, f"{system}_score.csv")
            try:
                df_score = pd.read_csv(scorepath)
            except:
                print(f"{newline(ok)}    Scores file not found ({scorepath})!")
                ok = False
        
            try:
                ranks = df_score["rank"].max()
                assert ranks == len(lignames)
            except AssertionError:
                print(f"{newline(ok)}    Number of scores mismatches number of poses ({ranks} vs {len(lignames)}!")
                ok = False
                
            rmsd = ["rmsd_lig", "rmsd_flex", "rmsd_tot"]
            for r in rmsd:
                if np.isnan(df_score[r].to_numpy()).any():
                    print(f"{newline(ok)}    {r.upper()} contains NaN!")
                    ok = False

                if np.isinf(df_score[r].to_numpy()).any():
                    print(f"{newline(ok)}    {r.upper()} contains Inf!")
                    ok = False

            if ok:
                print("ok")
                fvalid.write(f"{dataset}/{system}\n")
            else:
                finvalid.write(f"{dataset}/{system}\n")