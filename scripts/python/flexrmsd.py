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
        print("overlap/seqid cutoff:", cutoff)
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

                Patoms_per_res[
                    (atom.getResnum(), atom.getIcode(), atom.getChid())
                ].append(idx)
                Catoms_per_res[
                    (atom.getResnum(), atom.getIcode(), atom.getChid())
                ].append(Cidx)

        if len(Catoms) == 0:
            continue

        # TODO: Symmetry correction?
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

    return MATCHrmsd, MATCHrmsd_mflex


if __name__ == "__main__":
    import sys

    pose = sys.argv[1]
    cognate = sys.argv[2]
    flex = sys.argv[3]

    rmsd, maxresrmsd = calc_pocket_rmsd(
        pose, cognate, flex, "", verbose=False, symm=True
    )

    print(f"{rmsd:.5f} {maxresrmsd:.5f}")
