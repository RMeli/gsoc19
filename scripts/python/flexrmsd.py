import prody
import numpy as np
import os

from collections import defaultdict

def calc_pocket_rmsd(pose, cognate, flex, root="", verbose=False):
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

    flex = prody.parsePDB(os.path.join(root, flex))

    # Match POSE and COGNATE based on sequence identity
    matches = []
    for cutoff in range(90, 0, -10):
        # can't just set a low cutoff since we'll end up with bad alignments
        # try a whole bunch of alignments to maximize the likelihood we get the right one
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

    MATCHrmsd = np.inf # RMSD for all residues
    for Pmap, Cmap, _, _ in matches:
        Patoms = []
        Catoms = []
        Patoms_per_res = defaultdict(list)
        Catoms_per_res = defaultdict(list)

        for i, idx in enumerate(Pmap.getIndices()):
            atom = Pmap[i]
            
            # Get index of residue number in flerx corresponding to current atom
            # Raises value error is no such atom is found
            #residx = flexResnums.index(atom.getResnum())
            residx = np.where(flexResnums == atom.getResnum())

            #print(residx, flexIcodes[residx], atom.getIcode(), flexChids[residx], atom.getChid())

            if (atom.getResnum() in flexResnums) and (atom.getIcode() in flexIcodes) and (atom.getChid() in flexChids):
                # Store index of current atom
                # This atom is in a flexible residue
                # The flexible residue is in the matched PMap
                Patoms.append(idx)

                #print(atom.getResnum(), atom.getIcode(), atom.getChid())

                Cidx = Cmap.getIndices()[i]
                Catoms.append(Cidx)

                Patoms_per_res[(atom.getResnum(), atom.getIcode(), atom.getChid())].append(idx)
                Catoms_per_res[(atom.getResnum(), atom.getIcode(), atom.getChid())].append(Cidx)

        if len(Catoms) == 0:
            continue

        # TODO: Symmetry correction?
        rmsd = prody.calcRMSD(pose[Patoms], cognate[Catoms])
        if rmsd < MATCHrmsd:
            MATCHrmsd = rmsd

            # Compute maximum per-redisude RMSD
            # Only when the RMSD is the lowest
            MATCHrmsd_mflex = -1 # Maximum RMSD between all residues
            for Pres, PatomsR in Patoms_per_res.items():
                CatomsR = Catoms_per_res[Pres]

                maxrmsd = prody.calcRMSD(pose[PatomsR], cognate[CatomsR])

                if maxrmsd > MATCHrmsd_mflex:
                    MATCHrmsd_mflex = maxrmsd

    return MATCHrmsd, MATCHrmsd_mflex

if __name__ == "__main__":
    import sys

    pose = sys.argv[1]
    cognate = sys.argv[2]
    flex = sys.argv[3]

    rmsd, maxresrmsd = calc_pocket_rmsd(pose, cognate, flex, "", verbose=False)

    print(f"{rmsd:.5f} {maxresrmsd:.5f}")
