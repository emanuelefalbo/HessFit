#!/usr/bin/env python3

"""
Canonical RIC representation for HessFit.

The RIC definition is taken from the original Gaussian reference
calculation and is kept fixed for all conformers.

IMPORTANT:
    - RIC topology comes from Gaussian .log
    - RIC dimensions come from Gaussian .fchk
    - RIC Hessian is read directly from Gaussian .fchk
    - Cartesian Hessians are NOT converted to RIC Hessians here.
"""

from dataclasses import dataclass
from typing import List
import numpy as np

import parser_gau as pgau


@dataclass
class CanonicalRIC:
    """Immutable definition of the RIC coordinate system."""

    natoms: int

    # Atom indices are kept in Gaussian/HessFit 1-based convention
    bonds: List[List[int]]
    angles: List[List[int]]
    dihedrals: List[List[int]]

    n_ric: int
    n_bonds: int
    n_angles: int
    n_dihedrals: int

    # Reference QM information
    elements: List[str]
    coordinates: np.ndarray

    # Gaussian RIC Hessian read directly from FCHK
    hessian_ric: np.ndarray

    # Optional reference force/gradient information
    internal_forces: np.ndarray = None


def build_canonical_ric(log_file, fchk_file):
    """
    Build the canonical RIC representation from the original
    Gaussian reference calculation.

    Parameters
    ----------
    log_file : str
        Gaussian .log containing the RIC topology.

    fchk_file : str
        Gaussian .fchk containing RIC dimensions and
        internal force constants.

    Returns
    -------
    CanonicalRIC
    """

    # ------------------------------------------------------------
    # Read Gaussian files
    # ------------------------------------------------------------

    log_lines = pgau.store_any_file(log_file)
    fchk_lines = pgau.store_any_file(fchk_file)

    # ------------------------------------------------------------
    # Read molecular information from FCHK
    # ------------------------------------------------------------

    fchk = pgau.read_fchk(fchk_file)

    natoms = fchk["natoms"]
    coordinates = np.asarray(fchk["coords"], dtype=float)
    elements = list(fchk["elements"])

    # ------------------------------------------------------------
    # Read RIC dimensions from FCHK
    #
    # This is already how HessFit currently determines:
    # [total RIC, bonds, angles, dihedrals]
    # ------------------------------------------------------------

    ric_list, internal_forces = pgau.read_RicDim_Grad(fchk_lines)

    n_ric = ric_list[0]
    n_bonds = ric_list[1]
    n_angles = ric_list[2]
    n_dihedrals = ric_list[3]

    # ------------------------------------------------------------
    # Read RIC topology from Gaussian LOG
    # ------------------------------------------------------------

    bonds, angles, dihedrals = pgau.read_Top(
        log_lines,
        ric_list
    )

    # ------------------------------------------------------------
    # Basic consistency checks
    # ------------------------------------------------------------

    if len(bonds) != n_bonds:
        raise ValueError(
            f"RIC bond count mismatch: "
            f"log={len(bonds)}, fchk={n_bonds}"
        )

    if len(angles) != n_angles:
        raise ValueError(
            f"RIC angle count mismatch: "
            f"log={len(angles)}, fchk={n_angles}"
        )

    if len(dihedrals) != n_dihedrals:
        raise ValueError(
            f"RIC dihedral count mismatch: "
            f"log={len(dihedrals)}, fchk={n_dihedrals}"
        )

    # ------------------------------------------------------------
    # Read Gaussian INTERNAL RIC Hessian directly.
    #
    # DO NOT reconstruct this from the Cartesian Hessian.
    # ------------------------------------------------------------

    hessian_ric = pgau.read_HessRIC(
        fchk_lines,
        ric_list
    )

    expected_shape = (n_ric, n_ric)

    if hessian_ric.shape != expected_shape:
        raise ValueError(
            f"RIC Hessian shape mismatch: "
            f"{hessian_ric.shape} != {expected_shape}"
        )

    return CanonicalRIC(
        natoms=natoms,
        bonds=bonds,
        angles=angles,
        dihedrals=dihedrals,
        n_ric=n_ric,
        n_bonds=n_bonds,
        n_angles=n_angles,
        n_dihedrals=n_dihedrals,
        elements=elements,
        coordinates=coordinates,
        hessian_ric=hessian_ric,
        internal_forces=np.asarray(internal_forces, dtype=float)
    )

def validate_canonical_ric(ric):
    """
    Perform strict consistency checks on a CanonicalRIC object.
    """

    # Coordinate dimensions
    if ric.coordinates.shape != (ric.natoms, 3):
        raise ValueError(
            "Coordinate array has unexpected dimensions: "
            f"{ric.coordinates.shape}"
        )

    # Hessian dimensions
    if ric.hessian_ric.shape != (ric.n_ric, ric.n_ric):
        raise ValueError(
            "RIC Hessian has unexpected dimensions: "
            f"{ric.hessian_ric.shape}"
        )

    # Number of internal coordinates
    n_total = (
        ric.n_bonds +
        ric.n_angles +
        ric.n_dihedrals
    )

    if n_total != ric.n_ric:
        raise ValueError(
            "RIC dimension inconsistency: "
            f"{n_total} topology coordinates != "
            f"{ric.n_ric} total RICs"
        )

    # Atom index validation
    for group_name, group in [
        ("bonds", ric.bonds),
        ("angles", ric.angles),
        ("dihedrals", ric.dihedrals),
    ]:
        for indices in group:
            for idx in indices:
                if idx < 1 or idx > ric.natoms:
                    raise ValueError(
                        f"{group_name}: atom index {idx} "
                        f"outside 1..{ric.natoms}"
                    )

    # Hessian symmetry check
    asymmetry = np.max(
        np.abs(ric.hessian_ric - ric.hessian_ric.T)
    )

    if asymmetry > 1.0e-8:
        raise ValueError(
            f"RIC Hessian is not symmetric. "
            f"Maximum asymmetry = {asymmetry:.3e}"
        )

    return True