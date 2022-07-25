from pymol import cmd, stored

import os
import prody 

import prody
import numpy as np
import os
import warnings

import qcelemental as qcel

from spyrmsd.rmsd import symmrmsd
from spyrmsd.exceptions import NonIsomorphicGraphs

from collections import defaultdict


def _elements_to_atomicnums(elements: str):
    anums = np.zeros_like(elements, dtype=int)

    for i, e in enumerate(elements):
        anums[i] = qcel.periodictable.to_atomic_number(e)

    return anums


def _build_adjacency_matrix(selection):
    N = len(selection)
    A = np.zeros((N, N), dtype=int)
    for i, a1 in enumerate(selection.iterAtoms()):
        for j, a2 in enumerate(selection.iterAtoms()):
            for b1 in a1.iterBonded():
                if b1 == a2:
                    A[i, j] = 1
                    A[j, i] = 1
                    break

    return A

def calc_pocket_rmsd(pose, cognate, flex, root="", verbose=False, symm=True):
    """
    Calculate difference between the ligand reference receptor and
    the receptor it is being docked into.
    From original script by David Koes
    Parameters
    ----------
    Returns
    -------
    float
        RMSD
    """
    # Receptor pose (full receptor, reconstructed)
    pose = prody.parsePDB(os.path.join(root, pose))
    # Cognate receptor (full receptor, crystal structure)
    cognate = prody.parsePDB(os.path.join(root, cognate))

    # Automatically infer bonds for symmetry correction
    # Bonds are needed to build the adjacency matrix
    # BOTTLENECK!
    mbond = 2.1
    pose.inferBonds(max_bond=mbond)
    cognate.inferBonds(max_bond=mbond)

    # Some PDB files have the deprecates segment id column filled
    # This is not the case for the poses and therefore it can cause failures
    # Fill segment names from chain IDs (which should be correct)
    # See https://github.com/prody/ProDy/issues/1466
    pose.setSegnames(pose.getChids())
    cognate.setSegnames(cognate.getChids())

    flex = prody.parsePDB(os.path.join(root, flex))

    # Match POSE and COGNATE based on sequence identity
    matches = []
    for cutoff in range(90, 0, -10):
        # can't just set a low cutoff since we'll end up with bad alignments
        # try a whole bunch of alignments to maximize the likelihood we get the right one
        #print("overlap/seqid cutoff:", cutoff)
        m = prody.matchChains(
            pose, cognate, subset="all", overlap=cutoff, seqid=cutoff, pwalign=True
        )

        if m:
            matches += m

    # Print matches
    if verbose:
        for m in matches:
            print(m)

        print(flex.getResnums())
        print(flex.getIcodes())
        print(flex.getChids())

    flexResnums = flex.getResnums()
    flexIcodes = flex.getIcodes()
    flexChids = flex.getChids()

    atoms_per_residue = {}

    MATCHrmsd = np.inf  # RMSD for all residues

    # Need to distinguish residues in pose an cognate, they might not have same number/icode/chain
    matchingP = set() # Matching residues in pose
    matchingC = set() # Matching residues in cognate

    for Pmap, Cmap, _, _ in matches:
        Patoms = []
        Catoms = []
        Patoms_per_res = defaultdict(list)
        Catoms_per_res = defaultdict(list)

        for i, idx in enumerate(Pmap.getIndices()):
            atom = Pmap[i]

            # print(residx, flexIcodes[residx], atom.getIcode(), flexChids[residx], atom.getChid())

            if (
                (atom.getResnum() in flexResnums)
                and (atom.getIcode() in flexIcodes)
                and (atom.getChid() in flexChids)
            ):
                # Store index of current atom
                # This atom is in a flexible residue
                # The flexible residue is in the matched PMap
                Patoms.append(idx)

                # print(atom.getResnum(), atom.getIcode(), atom.getChid())

                Cidx = Cmap.getIndices()[i]
                Catoms.append(Cidx)

                # Both dictionaries are indexed by the same (resnum, icode, chain) as in the pose
                # This creates one-to-one correspondence between matching atoms/residues
                Patoms_per_res[
                    (atom.getResnum(), atom.getIcode(), atom.getChid())
                ].append(idx)
                Catoms_per_res[
                    (atom.getResnum(), atom.getIcode(), atom.getChid())
                ].append(Cidx)

                print(idx, Cidx, (pose[idx].getResnum(), pose[idx].getIcode(), pose[idx].getChid()), (cognate[Cidx].getResnum(), cognate[Cidx].getIcode(), cognate[Cidx].getChid()))

                matchingP.add((pose[idx].getResnum(), pose[idx].getIcode(), pose[idx].getChid()))
                matchingC.add((cognate[Cidx].getResnum(), cognate[Cidx].getIcode(), cognate[Cidx].getChid()))

        if len(Catoms) == 0:
            continue

        if not symm:
            rmsd = prody.calcRMSD(pose[Patoms], cognate[Catoms])
        else:
            Ap = _build_adjacency_matrix(pose[Patoms])
            Ac = _build_adjacency_matrix(cognate[Catoms])

            Panum = _elements_to_atomicnums(pose[Patoms].getElements())
            Canum = _elements_to_atomicnums(cognate[Catoms].getElements())

            try:
                rmsd = symmrmsd(
                    pose[Patoms].getCoords(),
                    cognate[Catoms].getCoords(),
                    Panum,
                    Canum,
                    Ap,
                    Ac,
                )
            except NonIsomorphicGraphs as e:  # Not isomorphic
                warnings.warn("NonIsomorphicGraphs | Computing standard RMSD")
                # Try computing RMSD without symmetry correction
                # This will possibly add some noise
                rmsd = prody.calcRMSD(pose[Patoms], cognate[Catoms]) 

        if rmsd < MATCHrmsd:
            MATCHrmsd = rmsd

            # Compute maximum per-redisude RMSD
            # Only when the RMSD is the lowest
            MATCHrmsd_mflex = -1  # Maximum RMSD between all residues
            for Pres, PatomsR in Patoms_per_res.items():
                CatomsR = Catoms_per_res[Pres]

                if not symm:
                    maxrmsd = prody.calcRMSD(pose[PatomsR], cognate[CatomsR])
                else:
                    ApR = _build_adjacency_matrix(pose[PatomsR])
                    AcR = _build_adjacency_matrix(cognate[CatomsR])

                    PanumR = _elements_to_atomicnums(pose[PatomsR].getElements())
                    CanumR = _elements_to_atomicnums(cognate[CatomsR].getElements())

                    try:
                        maxrmsd = symmrmsd(
                            pose[PatomsR].getCoords(),
                            cognate[CatomsR].getCoords(),
                            PanumR,
                            CanumR,
                            ApR,
                            AcR,
                        )
                    except NonIsomorphicGraphs as e:  # Not isomorphic
                        warnings.warn("NonIsomorphicGraphs | Computing standard RMSD")
                        # Try computing RMSD without symmetry correction
                        # This will possibly add some noise
                        maxrmsd = prody.calcRMSD(pose[PatomsR], cognate[CatomsR])

                if maxrmsd > MATCHrmsd_mflex:
                    MATCHrmsd_mflex = maxrmsd

    return MATCHrmsd, MATCHrmsd_mflex, matchingP, matchingC

def cdrmsd(
    ligand,
    protein,
    pocket,
    idx="0",
    flexdist="3.5", 
):
    # Clear everything
    cmd.reinitialize("everything")

    # Convert strings
    idx = int(idx)
    flexdist = float(flexdist)

    # Build paths
    path = ""
    ppath = os.path.join(path, pocket, "PDB_Structures")

    ligname = f"{protein}_PRO_{ligand}_LIG_aligned_v2_default_ensemble_none_flexdist3.5_p{idx}.sdf.gz"
    ligandpath = os.path.join(ppath, ligname)  # Docked ligand
    flexname = f"{protein}_PRO_{ligand}_LIG_aligned_v2_default_ensemble_none_flexdist3.5_flex_p{idx}.pdb.gz"
    flexpath = os.path.join(ppath, flexname)  # Flexible residues
    crystalname = f"{protein}_PRO.pdb"
    crystalpath = os.path.join(ppath, crystalname)  # Crystal receptor
    cognatename = f"{ligand}_PRO.pdb"
    cognatepath = os.path.join(ppath, cognatename) # Cognate receptor (crystal)
    cligname = f"{ligand}_LIG_aligned.sdf"
    cligpath = os.path.join(ppath, cligname)

    fRMSD, mfRMSD, matchingP, matchingC = calc_pocket_rmsd(flexpath.replace("_flex_", "_full_"), cognatepath, flexpath)
    print(f"flexrmsd = {fRMSD:.5f}\t fmaxrmsd = {mfRMSD:.5f}")
    print(matchingP, matchingC)

    # Load ligand and receptor
    cmd.load(ligandpath, "ligand")  # Selection name: ligand
    cmd.load(flexpath, "flex")  # Selection name: flex
    cmd.load(crystalpath, "crystal")  # Selection name: crystal
    cmd.load(cognatepath, "cognate") # selection name: cognate
    cmd.load(cligpath, "clig") # selection name: clig

    # Hide everything
    cmd.hide("all")

    # Show receptor
    cmd.show("cartoon", "crystal")
    cmd.color("lightorange", "crystal")

    # Show cognate receptro
    cmd.show("cartoon", "cognate")
    cmd.color("lightblue", "cognate")

    # Show docked ligand
    cmd.show("sticks", "ligand")
    cmd.show("spheres", "ligand")
    cmd.set("sphere_scale", 0.2, "ligand")
    cmd.color("limegreen", "ligand and name C*")

    # Show flexible residues
    cmd.show("licorice", "flex")
    cmd.set("stick_radius", 0.2, "flex")
    cmd.color("teal", "flex and name C*")

    cmd.show("sticks", "clig")
    cmd.color("grey", "clig")
            

    # Center and zoom to ligand
    cmd.center("ligand")
    cmd.zoom("ligand", 10)

    # Remove solvent
    cmd.remove("solvent")

    # Flexible residues in target
    for resi, icode, chid in matchingP:
        sel = f"flex and (resi {resi} in chain {chid})"
        cmd.show("sphere", sel)
        cmd.set("sphere_scale", 0.2, sel)
        cmd.color("marine", sel)
        cmd.remove(f"hydro in ({sel})")

    # Flexible residues in cognate
    for resi, icode, chid in matchingC:
        sel = f"cognate and (resi {resi} in chain {chid})"
        cmd.show("sticks", sel)
        cmd.show("sphere", sel)
        cmd.set("sphere_scale", 0.2, sel)
        cmd.color("orange", sel)
        cmd.remove(f"hydro in ({sel})")

cmd.extend("cdrmsd", cdrmsd)

