#!/usr/bin/env python3

"""
Generate Gaussian QM optimization/frequency jobs for all
HessFit candidate conformers.

The original Gaussian input is used as the master template.
Only:
    - checkpoint filename
    - title
    - Cartesian coordinates

are changed.

The Gaussian calculation itself remains identical for every
candidate conformer.
"""

import os
import glob
import re
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate Gaussian QM optimization/frequency jobs for HessFit conformers."
    )
    parser.add_argument(
        "--template",
        default="../examples/but_qm.gjf",
        help="Path to Gaussian template file (default: %(default)s)"
    )
    parser.add_argument(
        "--conformer-dir",
        default="conformers",
        help="Directory containing conformer XYZ files (default: %(default)s)"
    )
    # parser.add_argument(
    #     "--output-dir",
    #     default="gaussian_qm_conformers",
    #     help="Output directory for Gaussian input files (default: %(default)s)"
    # )
    return parser.parse_args()


# TEMPLATE = "../examples/but_qm.gjf"
# CONFORMER_DIR = "conformers"
# OUTDIR = "gaussian_qm_conformers"


def read_xyz(filename):

    with open(filename) as f:
        lines = [x.rstrip() for x in f]

    natoms = int(lines[0].strip())

    atoms = []
    coords = []

    for line in lines[2:2 + natoms]:

        fields = line.split()

        atoms.append(fields[0])

        coords.append(
            (
                float(fields[1]),
                float(fields[2]),
                float(fields[3])
            )
        )

    return natoms, atoms, coords


def find_geometry_block(lines):
    """
    Identify the Gaussian charge/multiplicity line.

    Returns
    -------
    charge_line_index
    """

    for i, line in enumerate(lines):

        stripped = line.strip()

        if re.match(
            r"^[+-]?\d+\s+[+-]?\d+$",
            stripped
        ):
            return i

    raise RuntimeError(
        "Could not identify Gaussian charge/multiplicity line."
    )


def find_coordinate_end(lines, charge_idx):
    """
    Determine where the existing Cartesian geometry ends.

    The original but_qm.gjf contains an empty line after the
    coordinate section.
    """

    start = charge_idx + 1

    for i in range(start, len(lines)):

        if not lines[i].strip():

            return i

    raise RuntimeError(
        "Could not identify end of Gaussian geometry."
    )


def create_input(
    template,
    xyz_file,
    output_file,
    chk_file,
    conformer_id
):

    with open(template) as f:
        lines = f.readlines()

    natoms, atoms, coords = read_xyz(xyz_file)

    charge_idx = find_geometry_block(lines)

    geom_end = find_coordinate_end(
        lines,
        charge_idx
    )

    # ------------------------------------------------------------
    # Keep everything from template through charge/multiplicity.
    # ------------------------------------------------------------

    new_lines = lines[:charge_idx + 1]

    # ------------------------------------------------------------
    # Replace geometry
    # ------------------------------------------------------------

    for atom, xyz in zip(atoms, coords):

        new_lines.append(
            f"{atom:2s}"
            f" {xyz[0]: .8f}"
            f" {xyz[1]: .8f}"
            f" {xyz[2]: .8f}\n"
        )

    # ------------------------------------------------------------
    # Keep the remainder of the template.
    #
    # This preserves connectivity or other sections if they
    # are present in future Gaussian templates.
    # ------------------------------------------------------------

    new_lines.extend(
        lines[geom_end:]
    )

    text = "".join(new_lines)

    # ------------------------------------------------------------
    # Replace checkpoint filename
    # ------------------------------------------------------------

    text = re.sub(
        r"(%chk=)[^\s]+",
        rf"\1{chk_file}",
        text,
        count=1,
        flags=re.IGNORECASE
    )

    # ------------------------------------------------------------
    # Replace title
    # ------------------------------------------------------------

    # The title is the line after the blank line following
    # the route section.
    #
    # For the current template this is simply "title".
    text = re.sub(
        r"(?m)^title\s*$",
        f"HessFit conformer {conformer_id}",
        text,
        count=1
    )

    with open(output_file, "w") as f:
        f.write(text)


def main():

    args = parse_arguments()

    OUTDIR = "gaussian_qm_conformers"
    CONFORMER_DIR = args.conformer_dir
    TEMPLATE = args.template

    os.makedirs(
        OUTDIR,
        exist_ok=True
    )

    xyz_files = sorted(
        glob.glob(
            os.path.join(
                CONFORMER_DIR,
                "conf_*.xyz"
            )
        )
    )

    if not xyz_files:

        raise RuntimeError(
            f"No conformers found in {CONFORMER_DIR}"
        )

    print(
        f"Found {len(xyz_files)} conformers."
    )

    for xyz_file in xyz_files:

        basename = os.path.basename(
            xyz_file
        )

        conformer_id = os.path.splitext(
            basename
        )[0]

        gjf_file = os.path.join(
            OUTDIR,
            conformer_id + ".gjf"
        )

        chk_file = (
            conformer_id +
            ".chk"
        )

        create_input(
            TEMPLATE,
            xyz_file,
            gjf_file,
            chk_file,
            conformer_id
        )

        print(
            f"Created {gjf_file}"
        )


if __name__ == "__main__":
    main()