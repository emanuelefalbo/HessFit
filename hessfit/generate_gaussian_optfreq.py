#!/usr/bin/env python3

"""
Generate Gaussian Opt + Freq input files from
xTB-selected conformers.

The XYZ geometries are used only as starting geometries.

Gaussian subsequently performs:
    geometry optimization
    +
    frequency calculation

The resulting optimized geometry and Hessian are
therefore both generated at the high-level QM method.

The atom ordering from the XYZ file is preserved.
"""

import os
import glob
import argparse


# ================================================================
# XYZ READER
# ================================================================

def read_xyz(filename):

    with open(filename) as f:
        lines = f.readlines()

    natoms = int(lines[0].strip())

    atoms = []

    for line in lines[2:2 + natoms]:

        fields = line.split()

        if len(fields) < 4:
            raise ValueError(
                f"Malformed XYZ line in {filename}: {line}"
            )

        atom = fields[0]

        x = float(fields[1])
        y = float(fields[2])
        z = float(fields[3])

        atoms.append(
            (atom, x, y, z)
        )

    if len(atoms) != natoms:
        raise ValueError(
            f"Expected {natoms} atoms in {filename}, "
            f"but found {len(atoms)}."
        )

    return atoms


# ================================================================
# GAUSSIAN INPUT
# ================================================================

def write_gaussian_input(
    xyz_file,
    output_file,
    method,
    basis,
    nproc,
    memory,
    charge,
    multiplicity,
    title
):

    atoms = read_xyz(
        xyz_file
    )

    with open(output_file, "w") as f:

        # --------------------------------------------------------
        # Gaussian resource directives
        # --------------------------------------------------------

        f.write(
            f"%mem={memory}\n"
        )

        f.write(
            f"%nprocshared={nproc}\n"
        )

        # Use the output filename without extension for the checkpoint file
        chk = os.path.splitext(output_file)[0] + ".chk"

        f.write(
            f"%chk={chk}\n"
        )

        # --------------------------------------------------------
        # Route section
        # --------------------------------------------------------

        f.write(
            f"#p {method}/{basis} "
            "opt(calcall,tight,maxstep=7) freq=intmodes\n"
        )

        f.write("\n")

        # --------------------------------------------------------
        # Title
        # --------------------------------------------------------

        f.write(
            f"{title}\n"
        )

        f.write("\n")

        # --------------------------------------------------------
        # Charge and multiplicity
        # --------------------------------------------------------

        f.write(
            f"{charge} {multiplicity}\n"
        )

        # --------------------------------------------------------
        # Geometry
        # --------------------------------------------------------

        for atom, x, y, z in atoms:

            f.write(
                f"{atom:<3s} "
                f"{x: .8f} "
                f"{y: .8f} "
                f"{z: .8f}\n"
            )

        f.write("\n")


# ================================================================
# MAIN
# ================================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate Gaussian Opt + Freq inputs "
            "from selected xTB conformers."
        )
    )

    parser.add_argument(
        "--input-dir",
        default="selected_conformers",
        help="Directory containing selected XYZ files."
    )

    parser.add_argument(
        "--output-dir",
        default="gaussian_optfreq",
        help="Directory for Gaussian inputs."
    )

    parser.add_argument(
        "--method",
        default="m052x",
        help="Gaussian method. Default: m052x"
    )

    parser.add_argument(
        "--basis",
        default="6-31+g(d,p)",
        help="Basis set."
    )

    parser.add_argument(
        "--nproc",
        type=int,
        default=48,
        help="Number of processors."
    )

    parser.add_argument(
        "--memory",
        default="96GB",
        help="Gaussian memory."
    )

    parser.add_argument(
        "--charge",
        type=int,
        default=0,
        help="Molecular charge."
    )

    parser.add_argument(
        "--multiplicity",
        type=int,
        default=1,
        help="Spin multiplicity."
    )

    args = parser.parse_args()

    os.makedirs(
        args.output_dir,
        exist_ok=True
    )

    xyz_files = sorted(
        glob.glob(
            os.path.join(
                args.input_dir,
                "*.xyz"
            )
        )
    )

    if not xyz_files:

        raise RuntimeError(
            f"No XYZ files found in {args.input_dir}"
        )

    print(
        f"Found {len(xyz_files)} selected conformers."
    )

    for xyz_file in xyz_files:

        basename = os.path.splitext(
            os.path.basename(xyz_file)
        )[0]

        output_file = os.path.join(
            args.output_dir,
            basename + ".com"
        )

        title = (
            f"HessFit high-level Opt+Freq "
            f"starting from {basename}"
        )

        write_gaussian_input(
            xyz_file=xyz_file,
            output_file=output_file,
            method=args.method,
            basis=args.basis,
            nproc=args.nproc,
            memory=args.memory,
            charge=args.charge,
            multiplicity=args.multiplicity,
            title=title
        )

        print(
            f"Written: {output_file}"
        )


if __name__ == "__main__":
    main()