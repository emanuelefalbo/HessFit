#!/usr/bin/env python3

"""
xTB conformational prescreening for HessFit.

xTB is used only to optimize candidate conformers and obtain
relative energies.

The resulting structures are subsequently used as starting
geometries for the high-level Gaussian Opt + Freq calculations.

The canonical HessFit RIC definition is NOT modified.
"""

import os
import subprocess
import shutil
import numpy as np


class XTBResult:

    def __init__(
        self,
        conformer_id,
        input_xyz,
        optimized_xyz,
        energy,
        success,
        directory
    ):

        self.conformer_id = conformer_id
        self.input_xyz = input_xyz
        self.optimized_xyz = optimized_xyz
        self.energy = energy
        self.success = success
        self.directory = directory


def run_xtb(
    xyz_file,
    output_dir,
    charge=0,
    uhf=0,
    method="gfn2",
    opt_level="tight",
    max_iterations=250
):
    """
    Optimize one candidate using xTB.

    Parameters
    ----------
    xyz_file
        Input Cartesian geometry.

    output_dir
        Directory for xTB calculation.

    charge
        Molecular charge.

    uhf
        Number of unpaired electrons.

    method
        xTB method. Default = GFN2-xTB.

    opt_level
        Optimization level.

    max_iterations
        Maximum optimization iterations.
    """

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    input_xyz = os.path.abspath(
        xyz_file
    )

    local_xyz = os.path.join(
        output_dir,
        "input.xyz"
    )

    shutil.copy2(
        input_xyz,
        local_xyz
    )

    command = [
        "xtb",
        "input.xyz",
        f"--{method}",
        "--opt",
        opt_level,
        "--chrg",
        str(charge),
        "--uhf",
        str(uhf),
        "--iterations",
        str(max_iterations)
    ]

    output_file = os.path.join(
        output_dir,
        "xtb.out"
    )

    with open(
        output_file,
        "w"
    ) as fout:

        process = subprocess.run(
            command,
            cwd=output_dir,
            stdout=fout,
            stderr=subprocess.STDOUT
        )

    optimized_xyz = os.path.join(
        output_dir,
        "xtbopt.xyz"
    )

    success = (
        process.returncode == 0
        and os.path.isfile(optimized_xyz)
    )

    energy = None

    if success:

        energy = extract_xtb_energy(
            output_file
        )

    return XTBResult(
        conformer_id=os.path.basename(
            xyz_file
        ),
        input_xyz=xyz_file,
        optimized_xyz=(
            optimized_xyz
            if success
            else None
        ),
        energy=energy,
        success=success,
        directory=output_dir
    )


def extract_xtb_energy(output_file):
    """
    Extract final total energy from xTB output.

    Returns energy in Hartree.
    """

    energy = None

    with open(output_file) as f:

        for line in f:

            # Typical xTB final energy line:
            #
            # | TOTAL ENERGY       -... Eh
            #
            if "TOTAL ENERGY" in line:

                fields = line.split()

                for i, field in enumerate(fields):

                    try:

                        value = float(field)

                    except ValueError:

                        continue

                    # Avoid integer fields and identify
                    # the energy from the context.
                    if i > 0:

                        energy = value

    return energy