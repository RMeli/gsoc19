"""
Create GNINA commands for flexible docking.

Based on:
    https://github.com/dkoes/GNINA-1.0/blob/main/analysis_scripts/make_gnina_cmds.py
"""

import argparse

parser = argparse.ArgumentParser(description="Create file with GNINA commands")
parser.add_argument(
    "input",
    help="Space-delimited file: <Receptor file> <Ligand file> <autobox ligand file> <outfile prefix>",
)
parser.add_argument(
    "-o",
    "--output",
    default="gnina_cmds.txt",
    type=str,
    help="Output file with GNINA commands",
)
parser.add_argument(
    "--cnn_scoring",
    default="rescore",
    type=str,
    choices=["none", "rescore", "refinement", "all"],
    help="CNN scoring",
)
parser.add_argument("--exhaustiveness", default=8, type=int, help="Exhaustiveness")
parser.add_argument("--num_modes", default=20, type=int, help="Number of poses")
parser.add_argument("--seed", default=42, type=int, help="RNG seed")
parser.add_argument(
    "--flexdist",
    default=3.5,
    type=float,
    help="Flexible residues distance. Uses <autobox_ligand> file as <flexdist_ligand>",
)
parser.add_argument("--n_cpus", default=1, type=int, help="Number of CPUs")
args = parser.parse_args()

# List of tuples: (recfile, ligfile, autobox_ligand, docked_prefix)
todock = []
with open(args.input) as infile:
    for line in infile:
        rec, lig, box, docked_prefix = line.strip().split()
        todock.append((rec, lig, box, docked_prefix))

with open(args.output, "w") as outfile:
    for r, l, box, out_prefix in todock:
        cmd = (
            f"gnina -r {r} -l {l}"
            + f" --autobox_ligand {box} --flexdist_ligand {box} --flexdist {args.flexdist}"
            + f" --cnn_scoring {args.cnn_scoring} --cpu {args.n_cpus} --seed {args.seed}"
            + f" --num_modes {args.num_modes} --exhaustiveness {args.exhaustiveness}"
        )

        lig_out = (
            out_prefix
            + "default_ensemble_"
            + args.cnn_scoring
            + f"_flexdist{args.flexdist}.sdf.gz"
        )

        cmd += f" --out {lig_out}"

        flex_out = (
            out_prefix
            + "default_ensemble_"
            + args.cnn_scoring
            + f"_flexdist{args.flexdist}_flex.sdf.gz"
        )

        cmd += f" --out_flex {flex_out}"

        outfile.write(cmd + "\n")
