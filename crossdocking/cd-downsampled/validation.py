import MDAnalysis as mda # Need version 2.0.0-beta or higher

from rdkit import Chem
from rdkit.Chem import Descriptors

import numpy as np
import pandas as pd

import os
import re
import warnings

import itertools
import gzip

cdpath = "carlos_cd"

residue_selection = "protein and not (type H or name H*)"
ligand_selection = "not type H and not (type H or name H*)"
water_selection = "not protein and (resname WAT or resname HOH)"


def newline(yes: bool) -> str:
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


with open("analysis/invalid.lst", "w", buffering=1) as finvalid, open(
    "analysis/valid.lst", "w", buffering=1
) as fvalid:

    with open("ds_cd_input_pairs.txt", "r") as fin:
        protein_ligand_pairs = fin.readlines()

    for line in protein_ligand_pairs:
        recpath, ligpath, _, _ = line.strip().split()

        rec = os.path.basename(recpath)
        lig = os.path.basename(ligpath)

        pocket = recpath.split("/")[1]
        
        ppath = os.path.join(cdpath, pocket, "PDB_Structures")

        ligid = lig[:4]
        recid = rec[:4]

        # Redocking is not performed
        if ligid == recid:
            continue

        ok: bool = True
        print(f"Validating {pocket}|{recid}_PRO_{ligid}_LIG...", end="")

        # Load protein crystal structure
        crecpath = os.path.join(ppath, rec)
        crec = load(crecpath)

        # Load ligand crystal structure
        cligpath = os.path.join(ppath, lig)
        clig = next(Chem.SDMolSupplier(cligpath, removeHs=False))
        if clig is None:
            print(f"{newline(ok)}    Failed loading {cligpath}")
            ok = False
            continue

        assert len(crec.atoms) > 0
        #n_residues = len(crec.residues)
        #assert n_residues != 0

        n_hvy_atoms_rec = len(crec.select_atoms(residue_selection))
        assert n_hvy_atoms_rec > 0

        n_hvy_atoms_lig = Descriptors.HeavyAtomCount(clig)
        assert n_hvy_atoms_lig > 0

        n_water = len(crec.select_atoms(water_selection))

        # !!! f-strings do not play well with regexp with {} !!!
        recnames = [
            r
            for r in os.listdir(ppath)
            if re.match(f"{recid}_PRO_{ligid}" + "_LIG_.*_full_p[0-9]{0,2}\.pdb\.gz", r)
        ]

        for recname in recnames:
            recpath = os.path.join(ppath, recname)

            try:
                rec = load(recpath)
            except Exception:
                print(f"{newline(ok)}    Failed loading {recpath}!")
                ok = False
                continue

            try:
                sel = rec.select_atoms(residue_selection)
                assert len(sel) == n_hvy_atoms_rec

                try:
                    # This does not really work when PDB files contain segment IDS
                    # For example: AKT1|3CQW_PRO_4EKK_LIG, AKT1|3CQW_PRO_4EKL_LIG, ...
                    # See: https://github.com/MDAnalysis/mdanalysis/issues/2874
                    #assert len(rec.residues) == n_residues

                    # Apply selection to crystal structure as well
                    # There might be a mistmatch on the number of hydrogens if the
                    # crystal structure had hidrogen on flexible residues, because
                    # GNINA only outputs polar hydrogens
                    csel = crec.select_atoms(residue_selection)
                    
                    # Compare selections (protein atoms but no hydrogens)
                    assert len(sel.atoms) == len(csel.atoms)
                    assert all(sel.atoms.resnums == csel.atoms.resnums)
                    assert all(sel.atoms.resids == csel.atoms.resids)
                    assert all(sel.atoms.chainIDs == csel.atoms.chainIDs)
                    assert all(sel.atoms.icodes == csel.atoms.icodes)

                except AssertionError:
                    print(
                        f"{newline(ok)}    Reference and reconstructed structures are different!\n    {recpath}"
                    )
                    ok = False

            except AssertionError:
                print(
                    f"{newline(ok)}    Wrong number of receptor heavy atoms ({len(sel)} vs {n_hvy_atoms_rec})!\n    {recpath}"
                )
                ok = False

            try:
                water = rec.select_atoms(water_selection)
                assert len(water) == n_water
            except AssertionError:
                print(
                    f"{newline(ok)}    Wrong number of water molecules ({len(water)} vs {n_water})!\n\t{recpath}"
                )
                ok = False

            if not ok:
                break

        lignames = [
            lig
            for lig in os.listdir(ppath)
            if re.match(f"{recid}_PRO_{ligid}" + "_LIG_.*_p[0-9]{0,2}\.sdf\.gz", lig)
        ]
        for ligname in lignames:
            ligpath = os.path.join(ppath, ligname)

            with gzip.open(ligpath, "r") as fgz:
                lig = next(Chem.ForwardSDMolSupplier(fgz, removeHs=False))
            if lig is None:
                print(f"{newline(ok)}    Failed loading {ligpath}")
                ok = False
                continue

            try:
                assert Descriptors.HeavyAtomCount(lig) == n_hvy_atoms_lig
            except:
                print(
                    f"{newline(ok)}    Wrong number of ligand heavy atoms ({Descriptors.HeavyAtomCount(lig)} vs {n_hvy_atoms_lig})!\n\t{ligpath}"
                )
                ok = False

            if not ok:
                break

        flexnames = [
            flex
            for flex in os.listdir(ppath)
            if re.match(f"{recid}_PRO_{ligid}" + "_LIG_.*_flex_p[0-9]{0,2}\.pdb\.gz", flex)
        ]
        for flexname in flexnames:
            flexpath = os.path.join(ppath, flexname)

            try:
                flex = load(flexpath)
            except:
                print(f"{newline(ok)}    Failed loading {flexpath}!\n\tNo flexible residues?")
                ok = False

            # Check there are no PRO residues within the flexible residues
            try:
                PRO = flex.select_atoms("protein and resname PRO")
                assert len(PRO) == 0
            except AssertionError:
                print(f"{newline(ok)}    Flexible PRO residues!\n\t{flexpath}")
                ok = False

            if not ok:
                break

        """         
        # Load score
        scorepath = os.path.join(dataset, system, f"{system}_score.csv")
        try:
            df_score = pd.read_csv(scorepath)
        except:
            print(f"{newline(ok)}    Scores file not found ({scorepath})!")
            ok = False

        try:
            ranks = df_score["rank"].max()

            # Take into account the possible presence of the crystal pose
            if df_score["rank"].min() == 0:
                ranks += 1

            assert ranks == len(lignames)
        except AssertionError:
            print(
                f"{newline(ok)}    Number of scores mismatches number of poses ({ranks} vs {len(lignames)})!"
            )
            ok = False

        rmsd = ["rmsd_lig", "rmsd_flex", "rmsd_tot"]
        for r in rmsd:
            if np.isnan(df_score[r].to_numpy()).any():
                print(f"{newline(ok)}    {r.upper()} contains NaN!")
                ok = False

            if np.isinf(df_score[r].to_numpy()).any():
                print(f"{newline(ok)}    {r.upper()} contains Inf!")
                ok = False
        """

        if ok:
            print("ok")
            fvalid.write(f"{pocket}\t{recid}\t{ligid}\n")
        else:
            finvalid.write(f"{pocket}\t{recid}\t{ligid}\n")