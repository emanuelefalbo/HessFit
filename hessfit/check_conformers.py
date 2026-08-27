#!/usr/bin/env python3

import argparse
import glob
import numpy as np

from ric_reference import (
    build_canonical_ric,
    validate_canonical_ric
)

from conformer_generation import read_xyz
from ric_geometry import evaluate_ric


parser = argparse.ArgumentParser(
    description="Check conformers against a canonical RIC reference."
)
parser.add_argument(
    "--log",
    default="template.log",
    help="Quantum chemistry log file to read (default: %(default)s)",
)
parser.add_argument(
    "--fchk",
    default="template.fchk",
    help="Formatted checkpoint file to read (default: %(default)s)",
)
args = parser.parse_args()

LOG = args.log
FCHK = args.fchk


ric = build_canonical_ric(
    LOG,
    FCHK
)

validate_canonical_ric(ric)


files = sorted(
    glob.glob("conformers/conf_*.xyz")
)

print(
    f"Found {len(files)} conformers."
)


for filename in files:

    atoms, coords = read_xyz(filename)

    q = evaluate_ric(
        coords,
        ric
    )

    print(
        f"{filename}: "
        f"{len(q)} RICs"
    )

    # First five bonds
    print(
        "  bonds:",
        " ".join(
            f"{x:.3f}"
            for x in q[:ric.n_bonds]
        )
    )

    # First five dihedrals
    start = ric.n_bonds + ric.n_angles

    print(
        "  dihedrals:",
        " ".join(
            f"{x:.2f}"
            for x in q[start:start + 5]
        )
    )