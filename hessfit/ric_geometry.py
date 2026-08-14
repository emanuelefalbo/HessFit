#!/usr/bin/env python3

"""
Geometry evaluation using a fixed HessFit canonical RIC definition.

The RIC topology is NEVER regenerated from the conformer.
All conformers are evaluated using the RIC ordering extracted from
the original Gaussian reference calculation.
"""

import numpy as np


def distance(coords, i, j):
    """Bond distance in Angstrom."""

    i -= 1
    j -= 1

    return np.linalg.norm(coords[i] - coords[j])


def angle(coords, i, j, k):
    """Angle i-j-k in degrees."""

    i -= 1
    j -= 1
    k -= 1

    v1 = coords[i] - coords[j]
    v2 = coords[k] - coords[j]

    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)

    cosang = np.clip(np.dot(v1, v2), -1.0, 1.0)

    return np.degrees(np.arccos(cosang))


def dihedral(coords, i, j, k, l):
    """
    Dihedral i-j-k-l in degrees.

    Convention follows the usual signed torsion definition.
    """

    i -= 1
    j -= 1
    k -= 1
    l -= 1

    p0 = coords[i]
    p1 = coords[j]
    p2 = coords[k]
    p3 = coords[l]

    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return np.degrees(np.arctan2(y, x))


def evaluate_ric(coords, ric):
    """
    Evaluate all canonical RICs for a Cartesian geometry.

    Returns
    -------
    q : ndarray
        RIC vector in exactly the same ordering as the Gaussian
        RIC Hessian.
    """

    coords = np.asarray(coords, dtype=float)

    if coords.shape != (ric.natoms, 3):
        raise ValueError(
            f"Expected coordinates with shape "
            f"({ric.natoms}, 3), got {coords.shape}"
        )

    values = []

    # ------------------------------------------------------------
    # Bonds
    # ------------------------------------------------------------

    for i, j in ric.bonds:
        values.append(distance(coords, i, j))

    # ------------------------------------------------------------
    # Angles
    # ------------------------------------------------------------

    for i, j, k in ric.angles:
        values.append(angle(coords, i, j, k))

    # ------------------------------------------------------------
    # Dihedrals
    # ------------------------------------------------------------

    for i, j, k, l in ric.dihedrals:
        values.append(dihedral(coords, i, j, k, l))

    q = np.asarray(values, dtype=float)

    if len(q) != ric.n_ric:
        raise ValueError(
            f"Generated {len(q)} RIC values but expected "
            f"{ric.n_ric}"
        )

    return q