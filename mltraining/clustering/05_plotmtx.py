import pickle
import argparse as ap

from matplotlib import pyplot as plt

parser = ap.ArgumentParser(description="Protein sequence distance and ligand similarity visualisation.")
parser.add_argument("input", type=str, help="Input (pickle)")
args = parser.parse_args()

dist, targets, lsim = pickle.load(open(args.input, "rb"))

plt.figure(1)
plt.imshow(dist)
plt.colorbar()
plt.title("Protein Sequence Distance")
plt.savefig("analysis/protseqdist.png", dpi=600)
plt.close()

plt.figure(2)
plt.imshow(lsim)
plt.colorbar()
plt.title("Ligand Similarity")
plt.savefig("analysis/ligsim.png", dpi=600)
plt.close()