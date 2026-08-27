#!/usr/bin/env python3

import argparse
import os
import sys

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(__file__),
        "hessfit"
    )
)

from ric_reference import (
    build_canonical_ric,
    validate_canonical_ric
)

from conformer_generation import (
    read_xyz,
    write_xyz,
    generate_conformers
)


DEFAULT_OUTDIR = "conformers"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate conformers from Gaussian log/fchk/xyz inputs."
    )
    parser.add_argument("log", nargs="?")
    parser.add_argument("fchk", nargs="?")
    parser.add_argument("xyz", nargs="?")
    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help=f"Directory for generated conformers (default: {DEFAULT_OUTDIR})"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    log = args.log
    fchk = args.fchk
    xyz = args.xyz
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)

    print("Building canonical RIC...")

    ric = build_canonical_ric(
        log,
        fchk
    )

    validate_canonical_ric(ric)

    print("Canonical RIC:")
    print("  atoms      :", ric.natoms)
    print("  RICs       :", ric.n_ric)
    print("  bonds      :", ric.n_bonds)
    print("  angles     :", ric.n_angles)
    print("  dihedrals  :", ric.n_dihedrals)

    atoms, coords = read_xyz(xyz)

    if len(atoms) != ric.natoms:

        raise ValueError(
            f"XYZ contains {len(atoms)} atoms, "
            f"but Gaussian RIC contains {ric.natoms}"
        )

    conformers = generate_conformers(
        atoms,
        coords,
        ric,
        n_conformers=20,
        max_rotation=180.0,
        min_rmsd=0.20,
        seed=12345
    )

    for i, conf in enumerate(conformers):

        filename = os.path.join(
            outdir,
            f"conf_{i:03d}.xyz"
        )

        write_xyz(
            filename,
            atoms,
            conf,
            comment=f"HessFit conformer {i:03d}"
        )

    print(
        f"\nStructures written to: {outdir}/"
    )


if __name__ == "__main__":
    main()