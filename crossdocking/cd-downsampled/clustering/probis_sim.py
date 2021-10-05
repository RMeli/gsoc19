"""
Based on calprobissim.py written by Paul Francoeur.

Please cite the following paper if you are using this script:

    Paul Francoeur, et al. Three-Dimensional Convolutional Neural Networks and a
    CrossDocked Data Set for Structure Based Drug Design. J. Chem. Inf. Model. (2020)
    DOI: https://doi.org/10.1021/acs.jcim.0c00411

Some functions are extracted from the code associated to the following paper:

    Rocco Meli, et al. Learning protein-ligand binding affinity with atomic environment
    vectors. Journal of Cheminformatics 13, 59 (2021)

"""

from typing import Union, List
import warnings
import os
import subprocess
import glob

import numpy as np
import qcelemental as qcel

import MDAnalysis as mda
from openbabel import pybel
import prody


def _universe_from_openbabel(obmol):
    """
    Create MDAnalysis universe from molecule parsed with OpenBabel.
    Parameters
    ----------
    obmol:
        Open Babel molecule
    Returns
    -------
    MDAnalysis universe

    Notes
    -----
    The molecule has resnum/resis set to 1, resname set to LIG and record type
    set to HETATM.

    Copyright 2020 Rocco Meli

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
    following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
    products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
    INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
    WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
    THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    """
    n_atoms = len(obmol.atoms)
    n_residues = 1  # LIG

    u = mda.Universe.empty(
        n_atoms,
        n_residues,
        atom_resindex=[0] * n_atoms,
        residue_segindex=[0] * n_residues,
        trajectory=True,
    )

    elements = []
    coordinates = np.zeros((n_atoms, 3))
    for idx, atom in enumerate(obmol):
        elements.append(qcel.periodictable.to_E(atom.atomicnum))
        coordinates[idx, :] = atom.coords

    # Complete records are needed for merging with protein PDB file
    u.add_TopologyAttr("elements", elements)
    u.add_TopologyAttr("type", elements)
    u.add_TopologyAttr("name", elements)
    u.add_TopologyAttr("resnum", [1] * n_residues)
    u.add_TopologyAttr("resid", [1] * n_residues)
    u.add_TopologyAttr("resname", ["LIG"] * n_residues)
    u.add_TopologyAttr("record_types", ["HETATM"] * n_atoms)
    u.add_TopologyAttr("segid", [""] * n_residues)

    u.atoms.positions = coordinates

    return u


def load_mols(
    ligand: str, receptor: str, datapaths: Union[str, List[str]]
) -> List[mda.Universe]:
    """
    Load ligand and receptor PDB files in a single mda.Universe
    Parameters
    ----------
    ligand: str
        Ligand file (SDF)
    receptor: str
        Receptor file (PDB)
    datapaths: Union[str, List[str]]
        Paths to root directory ligand and receptors are stored
    Returns
    -------
    Lit[mda.Universe]
        MDAnalysis universe for the protein-ligand complex
    Notes
    -----
    This function allows to load multiple ligands from SDF files for a single receptor.
    This is useful for docking and virtual screening, where multiple ligands are
    associated to a single target.
    The ligand is treated as a single entity named LIG. (:code:`resname LIG`).
    The folders containing ligand and receptor files data are defined by
    :code:`datapaths`.

    Copyright 2020 Rocco Meli

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
    following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
    products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
    INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
    WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
    THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    """

    ext = os.path.splitext(ligand)[-1].lower()[1:]

    # Assumes receptor is a PDB file
    # TODO: relax this assumption
    assert os.path.splitext(receptor)[-1].lower() == ".pdb"

    # Ensure list
    if isinstance(datapaths, str):
        datapaths = [datapaths]

    # TODO: Redirect warning instead of suppressing
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Try to load ligand
        for path in datapaths:
            ligfile = os.path.join(path, ligand)
            if os.path.isfile(ligfile):
                try:
                    uligs = [
                        _universe_from_openbabel(obmol)
                        for obmol in pybel.readfile(ext, ligfile)
                    ]
                except Exception:
                    print(f"Problems loading {ligfile}")
                    raise

                # Ligand file found in current path, no need to search further
                break
        else:
            raise RuntimeError(
                f"Could not find ligand file {ligfile} in {datapaths}..."
            )

        # Try to load receptor
        for path in datapaths:
            recfile = os.path.join(path, receptor)
            if os.path.isfile(recfile):
                try:
                    urec = mda.Universe(recfile)
                except Exception:
                    print(f"Problems loading {recfile}")
                    raise

                break
        else:
            raise RuntimeError(
                f"Could not find receptor file {recfile} in {datapaths}..."
            )

    ligs = [ulig.select_atoms("all") for ulig in uligs]
    rec = urec.select_atoms("protein")

    # Merge receptor and ligand in single universe
    systems = [mda.core.universe.Merge(lig, rec) for lig in ligs]

    return systems

def probis_sim(systems):
    """
    Compute similarity matrix between pockets using ProBiS.
    """
    n = len(systems)
    
    similarity_mtx  = -1 * np.ones((n,n))

    dist = 7.0

    for i in range(n):
        pocket1, pdbid1 = systems[i]

        # Load previously created file
        recname1 = f"pdbs/{pocket1}_{pdbid1}.pdb"

        p1 = prody.parsePDB(recname1)
        l1 = p1.select("resname LIG")

        lc1 = list(set(l1.getChids()))[0]
        rc1 = set(
            p1.select(
                f"within {dist} of (resname LIG and chain {lc1})"
            ).getChids()
        )
        ln1 = list(set(l1.getResnums()))

        for j in range(i, n):
            pocket2, pdbid2 = systems[j]

            recname2 = f"pdbs/{pocket2}_{pdbid2}.pdb"

            p2 = prody.parsePDB(recname2)

            l2 = p2.select("resname LIG")
            lc2 = list(set(l2.getChids()))[0]
            rc2  = set(
                p2.select(
                    f"within {dist} of (resname LIG and chain {lc2})"
                ).getChids()
            )
            ln2 = list(set(l2.getResnums()))

             # Start the string for the ProBIS command
            probis = f"./probis -compare -super -dist {dist} -f1 {recname1} -c1 "
            for c in rc1:
                probis += c
            probis += f" -bsite1 LIG.{ln1[0]}.{lc1} -f2 {recname2} -bsite2 LIG.{ln2[0]}.{lc2} -c2 "

            # Checking chains needed for probis for the second protein
            for c in rc2:
                probis += c

            # Make sure probis only runs on 1 CPU
            probis += " -ncpu 1"

            subprocess.call(probis, shell=True)

            # After running probis, assume proteins are not similar
            # Check output of probis & determine if they are similar
            # Similar if Z_SCORE > 2. *Asked developers if this will work
            if not glob.glob(f"{pocket1}*_{pocket2}*.0.rota.pdb"):
                similarity_mtx[i, j] = 0
                similarity_mtx[j, i] = 0
                print(f"Pockets {pocket1} and {pocket2} are not similar.")
                print(f"Pocket representatives: {recname1} and {recname2} are not similar.")
            else:
                esent = subprocess.check_output(
                    f"grep Z_SCORE {pocket1}*_{pocket2}*.0.rota.pdb",
                    shell=True,
                )
                esent = esent.split()

                Zscore = float(esent[3])
                if Zscore > 2:
                    similarity_mtx[i, j] = 1
                    similarity_mtx[j, i] = 1
                    print(f"Pockets {pocket1} and {pocket2} are similar.")
                    print(f"    Zscore: {Zscore}.")
                else:
                    similarity_mtx[i, j] = 0
                    similarity_mtx[j, i] = 0
                    print(f"Pockets {pocket1} and {pocket2} are not similar.")
                    print(f"    Pocket representatives: {recname1} and {recname2} are not similar.")
                    print(f"    Zscore: {Zscore}.")


    return similarity_mtx


if __name__ == "__main__":

    root = "../carlos_cd"

    # Loop over systems and create correct protein-ligand PDB file
    # This file will be loaded in with  ProDy in order to re-use most of Paul's code
    systems = []
    with open("representatives.list", "r") as fin:
        for line in fin:
            splitted_line = line.strip().split()
            pocket = splitted_line[0]
            pdbid = splitted_line[1]
            systems.append((pocket, pdbid))

            datapath = os.path.join(root, pocket, "PDB_Structures")
            ligname = f"{pdbid}_LIG_aligned.sdf"
            recname = f"{pdbid}_PRO.pdb"

            # Name of the output file
            outfile = f"pdbs/{pocket}_{pdbid}.pdb"

            # Create file if it does not exist
            if not os.path.exists(outfile):
                # Create unified PDB file wit ligand and protein
                # Ligand name is set to LIG
                # Some receptors have not been included in the GNINA1.0 paper because
                # of subsampling.
                # First look at the subsampled dataset, then look at the original dataset
                systems = load_mols(
                    ligname,
                    recname,
                    datapaths=[datapath, os.path.join("disco", pocket, "PDB_Structures")],
                )

                assert len(systems) == 1

                systems[0].select_atoms("all").write()

    sim = probis_sim(systems)

    print(sim)
