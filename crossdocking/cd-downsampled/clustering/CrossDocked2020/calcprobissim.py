"""
Script originally written by Paul Francoeur.
"""

import os
import sys
import glob

# from multiprocessing import Pool
import numpy as np
from prody import *
import subprocess
import pickle

"""
This script takes in 2 arguments from the command-line: 1) Name of pocket dir, 2) largest ligand in the pocket
This script assumes that there exists a file dtup.pkl which is a pickeled python list of (pocket_dir, ligname)
This script then calculates the similarity of the pocket from the sys arg with every pocket in the list.
The similarity calculation returns a list starting with the pocket itself, & then the pockets it is similar to

Ex) if a,b,c,x,y,z are the pockets with a being similar to b,c, then this script when passed the things for a will write a file
    a.out --> which would be one line containing 'a b c'
"""


def helper(dtup_in):
    """
    Depreciated function for older version of script
    """
    rec, lig = dtup_in
    return probis_sim(rec, lig, d_tup)


def getzips():
    """
    Function to move & unzip all of the needed pdbs from the pdb folder used from PRODY
    """

    # move all the files in
    for rec, lig in d_tup:
        rname = lig.split("_")[0]
        subprocess.call("cp " + pdb_path + rname + "* .", shell=True)

    # unzip the files
    subprocess.call("gunzip *.gz", shell=True)


def probis_sim(pocket, lig, list2compare):
    """
    This function uses probis to calculate the similarity of given pocket with all the other pockets.
    Returns: list containing pocket itself, then other pockets it is similar to
    """
    pdb_path = "~/pdb/"
    ret = []
    ret.append(pocket)
    rname = lig.split("_")[0]
    p1 = parsePDB(rname)
    lname = lig.split("_")[1].upper()
    dist = "7"
    l1 = p1.select("resname " + lname)
    lc1 = list(set(l1.getChids()))[0]
    rc1 = set(
        p1.select(
            "within " + dist + " of (resname " + lname + " and chain " + lc1 + ")"
        ).getChids()
    )
    ln1 = list(set(l1.getResnums()))
    for pocket2, lig2 in list2compare:
        rname2 = lig2.split("_")[0]
        if rname2 == rname:
            continue
        p2 = parsePDB(rname2)
        lname2 = lig2.split("_")[1].upper()
        l2 = p2.select("resname " + lname2)
        lc2 = list(set(l2.getChids()))[0]
        rc2 = set(
            p2.select(
                "within " + dist + " of (resname " + lname2 + " and chain " + lc2 + ")"
            ).getChids()
        )
        ln2 = list(set(l2.getResnums()))
        # start the string for the ProBIS command
        probis = "probis -compare -super -dist 7.0 -f1 " + rname + ".pdb -c1 "
        # checking chains needed for probis for the first protein
        for c in rc1:
            probis += c
        probis += (
            " -bsite1 "
            + lname
            + "."
            + str(ln1[0])
            + "."
            + lc1
            + " -f2 "
            + rname2
            + ".pdb"
            + " -bsite2 "
            + lname2
            + "."
            + str(ln2[0])
            + "."
            + lc2
            + " -c2 "
        )
        # checking chains needed for probis for the second protein
        for c in rc2:
            probis += c
        # make sure probis only runs on 1 CPU
        probis += " -ncpu 1"
        subprocess.call(probis, shell=True)
        # after running probis, assume proteins are not similar
        # check output of probis & deetermine if they are if Z_SCORE > 2. **asked developers if this will work
        if not glob.glob(rname + "*" + "_" + rname2 + "*.0.rota.pdb"):
            pass
            # print "pockets not similar"
        else:
            esent = subprocess.check_output(
                "grep Z_SCORE " + rname + "*" + "_" + rname2 + "*.0.rota.pdb",
                shell=True,
            )
            esent = esent.split()
            if float(esent[3]) > 2:
                # print "pocket is similar"
                ret.append(pocket2)
    return ret


# START OF PROGRAM
d_tup = pickle.load(open("dtup.pkl", "rb"))
# getzips()
# pool=Pool(processes=63)
# probis_sim_list=pool.map(helper,d_tup)
# with open("probis_sim_list.pkl","wb") as f:
# 	pickle.dump(probis_sim_list,f,pickle.HIGHEST_PROTOCOL)
pocket_dir = sys.argv[1]
lig_pdb = sys.argv[2]
sim_list = probis_sim(pocket_dir, lig_pdb, d_tup)
with open(pocket_dir + ".out", "w") as f:
    for item in sim_list:
        f.write(item + " ")
