#!/usr/bin/env python3
"""
evaluate_hessfit_conformers.py

Evaluate selected conformers after HessFit --test.

For each conformer:

    QM:
        <basename>.log
        <basename>.fchk

    HessFit:
        <basename>_hessfit.log
        <basename>_hessfit.fchk

Three quantities are evaluated:

    1. Normalized RIC Hessian error
    2. Geometry RMSD
    3. Vibrational-frequency error

The final merit function is:

    M = wH * H_norm
      + wR * RMSD_norm
      + wF * Freq_norm

The three quantities are min-max normalized over the selected
conformers before calculating the final merit.

The conformer with the lowest merit is selected.

IMPORTANT
---------
The RIC Hessians MUST be generated using exactly the same RIC
definition used by HessFit.

This script therefore contains an adapter function,
load_ric_hessian(), which should be connected to the existing
HessFit FCHK/RIC routines.

No artificial MM Hessian is reconstructed here.
"""

import argparse
import csv
import os
import re
import parser_gau as pgau
import numpy as np


# ================================================================
# Constants
# ================================================================

HARTREE_TO_KCAL = 627.509474
EPS = 1.0e-12


# ================================================================
# Generic file utilities
# ================================================================

def read_text(filename):

    with open(filename, "r") as f:
        return f.read()


# ================================================================
# Gaussian LOG
# ================================================================

def check_gaussian_termination(logfile):

    text = read_text(logfile)

    return "Normal termination of Gaussian" in text


def extract_fchk_geometry(fchk_file):
    """
    Extract atomic numbers and Cartesian coordinates from
    a Gaussian formatted checkpoint file.

    Coordinates are returned in Angstrom.
    """

    lines = read_text(fchk_file).splitlines()

    atomic_numbers = None
    coordinates = None

    # ------------------------------------------------------------
    # Atomic numbers
    # ------------------------------------------------------------

    for i, line in enumerate(lines):

        if line.startswith("Atomic numbers"):

            parts = line.split()

            # Number of atoms
            n_atoms = int(parts[-1])

            values = []

            j = i + 1

            while len(values) < n_atoms:

                values.extend(
                    int(x)
                    for x in lines[j].split()
                )

                j += 1

            atomic_numbers = np.asarray(
                values[:n_atoms],
                dtype=int
            )

            break

    if atomic_numbers is None:

        raise RuntimeError(
            "Could not find Atomic numbers in {}".format(
                fchk_file
            )
        )

    # ------------------------------------------------------------
    # Cartesian coordinates
    # ------------------------------------------------------------

    for i, line in enumerate(lines):

        if line.startswith("Current cartesian coordinates"):

            parts = line.split()

            n_values = int(parts[-1])

            values = []

            j = i + 1

            while len(values) < n_values:

                values.extend(
                    float(x)
                    for x in lines[j].split()
                )

                j += 1

            coordinates = np.asarray(
                values[:n_values],
                dtype=float
            ).reshape((-1, 3))

            break

    if coordinates is None:

        raise RuntimeError(
            "Could not find Current cartesian coordinates "
            "in {}".format(fchk_file)
        )

    # ------------------------------------------------------------
    # Gaussian FCHK coordinates are in Bohr.
    # Convert to Angstrom.
    # ------------------------------------------------------------

    BOHR_TO_ANGSTROM = 0.529177210903

    coordinates *= BOHR_TO_ANGSTROM

    if len(coordinates) != len(atomic_numbers):

        raise RuntimeError(
            "Mismatch between atomic numbers and coordinates "
            "in {}".format(fchk_file)
        )

    return atomic_numbers, coordinates


def extract_final_geometry(logfile):
    """
    Extract the final geometry from a Gaussian log file.

    Gaussian coordinate rows have the format:

        Center Atomic Atomic X Y Z

    Example:

        1  6  0  -0.706035   1.820058   0.000000

    The parser does not depend on the presence of the Gaussian
    table header. It identifies coordinate rows directly.

    Returns
    -------
    atomic_numbers : np.ndarray
    coordinates : np.ndarray
    """

    lines = read_text(logfile).splitlines()

    # ------------------------------------------------------------
    # Locate all Standard orientation sections
    # ------------------------------------------------------------

    orientation_indices = []

    for i, line in enumerate(lines):

        if "Standard orientation:" in line:
            orientation_indices.append(i)

    if not orientation_indices:

        raise RuntimeError(
            "No 'Standard orientation:' found in {}".format(logfile)
        )

    # ------------------------------------------------------------
    # Parse a Standard orientation block
    # ------------------------------------------------------------

    def parse_block(start):

        atomic_numbers = []
        coordinates = []

        # Search only until the next major Gaussian section.
        #
        # A normal geometry table is only a few tens of lines long.
        max_end = min(start + 200, len(lines))

        for i in range(start + 1, max_end):

            line = lines[i].strip()

            if not line:
                continue

            # Separator
            if line.startswith("-"):
                continue

            fields = line.split()

            # A Gaussian coordinate row must contain at least:
            #
            # atom_index
            # atomic_number
            # atomic_type
            # x
            # y
            # z
            #
            if len(fields) < 6:
                continue

            try:

                atom_index = int(fields[0])
                atomic_number = int(fields[1])
                atomic_type = int(fields[2])

                x = float(fields[3])
                y = float(fields[4])
                z = float(fields[5])

            except (ValueError, TypeError):

                # Once coordinate rows have started, a
                # non-coordinate line means that the table ended.
                if atomic_numbers:
                    break

                continue

            # ----------------------------------------------------
            # Sanity checks
            # ----------------------------------------------------

            if atom_index != len(atomic_numbers) + 1:

                if atomic_numbers:
                    break

                continue

            if atomic_number <= 0:
                if atomic_numbers:
                    break
                continue

            atomic_numbers.append(atomic_number)

            coordinates.append([x, y, z])

        # Require at least one atom
        if not coordinates:
            return None

        return (
            np.asarray(atomic_numbers, dtype=int),
            np.asarray(coordinates, dtype=float)
        )

    # ------------------------------------------------------------
    # Parse every Standard orientation
    # ------------------------------------------------------------

    geometries = []

    for start in orientation_indices:

        result = parse_block(start)

        if result is not None:

            atomic_numbers, coordinates = result

            geometries.append(
                (
                    start,
                    atomic_numbers,
                    coordinates
                )
            )

    if not geometries:

        raise RuntimeError(
            "Could not parse any Standard orientation block from {}".format(
                logfile
            )
        )

    # ------------------------------------------------------------
    # Find final optimization completion
    # ------------------------------------------------------------

    completion_indices = []

    for i, line in enumerate(lines):

        if "Optimization completed." in line:
            completion_indices.append(i)

    # ------------------------------------------------------------
    # Prefer the last Standard orientation associated with the
    # final optimization.
    # ------------------------------------------------------------

    if completion_indices:

        final_completion = completion_indices[-1]

        candidates = [
            geometry
            for geometry in geometries
            if geometry[0] < final_completion
        ]

        if candidates:

            _, atomic_numbers, coordinates = candidates[-1]

            return atomic_numbers, coordinates

    # ------------------------------------------------------------
    # Otherwise use the final successfully parsed geometry.
    # ------------------------------------------------------------

    _, atomic_numbers, coordinates = geometries[-1]

    return atomic_numbers, coordinates


def extract_frequencies(logfile):
    """
    Extract Gaussian harmonic frequencies.

    Imaginary frequencies are retained as negative values.
    """

    frequencies = []

    pattern = re.compile(
        r"Frequencies\s+--\s+(.*)$"
    )

    with open(logfile, "r") as f:

        for line in f:

            match = pattern.search(line)

            if match is None:
                continue

            for value in match.group(1).split():

                try:
                    frequencies.append(float(value))

                except ValueError:
                    pass

    if len(frequencies) == 0:

        raise RuntimeError(
            "No vibrational frequencies found in {}".format(
                logfile
            )
        )

    return np.asarray(frequencies, dtype=float)


def extract_final_frequencies(logfile, n_atoms):
    """
    Extract the final complete vibrational-frequency set from a
    Gaussian log file.

    Gaussian may contain multiple frequency calculations in the same
    log file. Each frequency set contains 3N-6 vibrational frequencies
    for a nonlinear molecule.

    The parser:
        1. collects every line containing "Frequencies --"
        2. extracts the numerical frequencies
        3. groups them into complete sets of 3N-6 frequencies
        4. returns the LAST complete set

    Parameters
    ----------
    logfile : str
        Gaussian log file.
    n_atoms : int
        Number of atoms.

    Returns
    -------
    np.ndarray
        Final complete set of vibrational frequencies in cm^-1.
    """

    expected_nfreq = 3 * n_atoms - 6

    all_frequencies = []

    with open(logfile, "r") as f:

        for line in f:

            if "Frequencies --" not in line:
                continue

            # Everything after "Frequencies --"
            # should contain the frequencies.
            part = line.split("Frequencies --", 1)[1]

            values = part.split()

            for value in values:

                try:
                    all_frequencies.append(float(value))
                except ValueError:
                    pass

    # ------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------

    print("\nFrequency extraction:")
    print("  Expected frequencies per set :", expected_nfreq)
    print("  Total frequencies found      :", len(all_frequencies))

    if len(all_frequencies) == 0:

        raise RuntimeError(
            "No vibrational frequencies found in {}".format(
                logfile
            )
        )

    # ------------------------------------------------------------
    # Check that the total number is compatible with complete sets
    # ------------------------------------------------------------

    n_complete_sets = len(all_frequencies) // expected_nfreq

    remainder = len(all_frequencies) % expected_nfreq

    print("  Complete frequency sets      :", n_complete_sets)
    print("  Remaining frequencies        :", remainder)

    if n_complete_sets == 0:

        raise RuntimeError(
            "Could not find a complete set of {} vibrational "
            "frequencies in {}".format(
                expected_nfreq,
                logfile
            )
        )

    # ------------------------------------------------------------
    # Extract complete sets
    # ------------------------------------------------------------

    frequency_sets = []

    for i in range(n_complete_sets):

        start = i * expected_nfreq
        end = start + expected_nfreq

        block = np.asarray(
            all_frequencies[start:end],
            dtype=float
        )

        frequency_sets.append(block)

    # ------------------------------------------------------------
    # Select LAST complete set
    # ------------------------------------------------------------

    final_frequencies = frequency_sets[-1]

    print("  Using frequency set          :", n_complete_sets)
    print("  Frequencies used             :", len(final_frequencies))

    return final_frequencies


# def extract_final_frequencies(logfile, n_atoms):
    """
    Extract the final complete set of vibrational frequencies
    from a Gaussian log file.

    Gaussian may contain multiple frequency calculations in
    the same log. Only the LAST complete frequency set is used.

    Parameters
    ----------
    logfile : str
        Gaussian log file.
    n_atoms : int
        Number of atoms.

    Returns
    -------
    frequencies : np.ndarray
        Vibrational frequencies in cm^-1.
    """

    lines = read_text(logfile).splitlines()

    expected_nfreq = 3 * n_atoms - 6

    frequency_sets = []
    current_set = []

    for line in lines:

        if "Frequencies --" in line:

            try:

                values = line.split("Frequencies --", 1)[1].split()


                print(values)

                frequencies = [
                    float(x)
                    for x in values
                ]

                current_set.extend(frequencies)

            except ValueError:
                continue

        else:

            # A non-frequency line terminates the current
            # contiguous frequency section.
            if current_set:

                if len(current_set) >= expected_nfreq:

                    frequency_sets.append(
                        np.asarray(
                            current_set[:expected_nfreq],
                            dtype=float
                        )
                    )

                current_set = []

    print(current_set)

    # Handle a frequency section that reaches EOF
    if current_set and len(current_set) >= expected_nfreq:

        frequency_sets.append(
            np.asarray(
                current_set[:expected_nfreq],
                dtype=float
            )
        )

    if not frequency_sets:

        raise RuntimeError(
            "Could not find a complete set of {} vibrational "
            "frequencies in {}".format(
                expected_nfreq,
                logfile
            )
        )

    # ------------------------------------------------------------
    # IMPORTANT:
    # Select the LAST frequency set
    # ------------------------------------------------------------

    final_frequencies = frequency_sets[-1]

    if len(final_frequencies) != expected_nfreq:

        raise RuntimeError(
            "Final frequency set contains {} frequencies; "
            "expected {}".format(
                len(final_frequencies),
                expected_nfreq
            )
        )

    return final_frequencies

def extract_scf_energy(logfile):
    """
    Extract the last SCF energy.
    """

    energies = []

    pattern = re.compile(
        r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+]?\d+\.\d+)"
    )

    with open(logfile, "r") as f:

        for line in f:

            match = pattern.search(line)

            if match:

                energies.append(
                    float(match.group(1))
                )

    if not energies:
        return np.nan

    return energies[-1]


# ================================================================
# Geometry RMSD
# ================================================================

def kabsch_rmsd(reference, target):
    """
    Calculate RMSD after optimal Kabsch alignment.

    Atom ordering must be identical.
    """

    reference = np.asarray(reference, dtype=float)
    target = np.asarray(target, dtype=float)

    if reference.shape != target.shape:

        raise ValueError(
            "Geometry dimensions differ: {} vs {}".format(
                reference.shape,
                target.shape
            )
        )

    ref_center = reference.mean(axis=0)
    tar_center = target.mean(axis=0)

    A = reference - ref_center
    B = target - tar_center

    covariance = np.dot(A.T, B)

    U, S, Vt = np.linalg.svd(covariance)

    determinant = np.linalg.det(
        np.dot(Vt.T, U.T)
    )

    correction = np.eye(3)

    if determinant < 0.0:
        correction[2, 2] = -1.0

    rotation = np.dot(
        np.dot(Vt.T, correction),
        U.T
    )

    A_rot = np.dot(A, rotation)

    difference = A_rot - B

    return np.sqrt(
        np.mean(
            np.sum(
                difference ** 2,
                axis=1
            )
        )
    )


def calculate_geometry_rmsd(
        qm_atomic_numbers,
        qm_coordinates,
        mm_atomic_numbers,
        mm_coordinates,
        heavy_atoms=True
):
    """
    Calculate QM/MM geometry RMSD.

    By default only heavy atoms are included.
    """

    if not np.array_equal(
        qm_atomic_numbers,
        mm_atomic_numbers
    ):

        raise ValueError(
            "QM and MM atomic-number ordering is different."
        )

    if heavy_atoms:

        mask = qm_atomic_numbers > 1

        qm_coordinates = qm_coordinates[mask]
        mm_coordinates = mm_coordinates[mask]

    return kabsch_rmsd(
        qm_coordinates,
        mm_coordinates
    )


# ================================================================
# Frequency error
# ================================================================

def calculate_frequency_error(
        qm_frequencies,
        mm_frequencies
):
    """
    Calculate relative RMS frequency error.

    Only positive frequencies are compared.

    This excludes translational/rotational zero modes and
    imaginary frequencies.
    """

    qm_frequencies = np.asarray(
        qm_frequencies,
        dtype=float
    )

    mm_frequencies = np.asarray(
        mm_frequencies,
        dtype=float
    )

    qm_positive = qm_frequencies[
        qm_frequencies > 0.0
    ]

    mm_positive = mm_frequencies[
        mm_frequencies > 0.0
    ]

    n_imag_qm = np.sum(
        qm_frequencies < 0.0
    )

    n_imag_mm = np.sum(
        mm_frequencies < 0.0
    )

    if len(qm_positive) != len(mm_positive):

        raise ValueError(
            "Different number of positive vibrational "
            "frequencies: QM={} MM={}".format(
                len(qm_positive),
                len(mm_positive)
            )
        )

    if len(qm_positive) == 0:

        raise ValueError(
            "No positive QM vibrational frequencies."
        )

    denominator = np.maximum(
        np.abs(qm_positive),
        EPS
    )

    relative_error = (
        mm_positive - qm_positive
    ) / denominator

    rms_error = np.sqrt(
        np.mean(
            relative_error ** 2
        )
    )

    return {
        "frequency_error": float(rms_error),
        "n_modes": int(len(qm_positive)),
        "n_imag_qm": int(n_imag_qm),
        "n_imag_mm": int(n_imag_mm)
    }


# ================================================================
# Gaussian FCHK reader
# ================================================================

def read_fchk_array(filename, label):
    """
    Read an array from a Gaussian formatted checkpoint file.

    This generic routine is provided for FCHK data that may be
    needed by the RIC adapter.
    """

    with open(filename, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):

        if not line.startswith(label):
            continue

        match = re.search(
            r"N=\s*(\d+)",
            line
        )

        if match is None:
            continue

        n_values = int(
            match.group(1)
        )

        values = []

        for next_line in lines[i + 1:]:

            for field in next_line.split():

                try:

                    values.append(
                        float(field)
                    )

                except ValueError:
                    continue

                if len(values) == n_values:

                    return np.asarray(
                        values,
                        dtype=float
                    )

    raise KeyError(
        "FCHK array '{}' not found in {}".format(
            label,
            filename
        )
    )


# ================================================================
# RIC Hessian adapter
# ================================================================

def load_ric_hessian(fchk_file, ric_list=None):
    """
    Read the RIC Hessian from an FCHK file.

    IMPORTANT
    ---------
    This function must be connected to the existing HessFit
    RIC-Hessian reader.

    The function must return:

        H_RIC

    as a square NumPy array:

        (N_RIC, N_RIC)

    Do NOT insert an independent Cartesian -> RIC transformation
    here.

    The QM and HessFit FCHK files must be transformed using the
    exact same RIC definition already used by HessFit.
    """

    # ============================================================
    # TODO:
    #
    # Replace this block with your existing HessFit routine.
    #
    # For example, if your HessFit code contains:
    #
    #     hessRIC = some_existing_function(filename)
    #
    # then simply use:
    #
    #     return some_existing_function(filename)
    #
    # ============================================================


    text_fchk = pgau.store_any_file(fchk_file)
    # ric_list, _ = pgau.read_RicDim_Grad(text_fchk)
    # print(ric_list)
    hessRIC = pgau.read_HessRIC_2(text_fchk, ric_list)

    return hessRIC
    
    raise NotImplementedError(
        "\n"
        "RIC Hessian reader is not yet connected.\n\n"
        "Connect load_ric_hessian() to the existing HessFit "
        "FCHK -> RIC Hessian routine.\n\n"
        "This is intentional: we must use exactly the same "
        "RIC definition as the production HessFit code."
    )


# ================================================================
# Hessian metric
# ================================================================

def calculate_hessian_error(
        qm_hessian,
        mm_hessian
):
    """
    Normalized Frobenius-norm Hessian error:

        ||H_QM - H_MM|| / ||H_QM||
    """

    qm_hessian = np.asarray(
        qm_hessian,
        dtype=float
    )

    mm_hessian = np.asarray(
        mm_hessian,
        dtype=float
    )

    if qm_hessian.shape != mm_hessian.shape:

        raise ValueError(
            "RIC Hessian dimensions differ: "
            "QM={} MM={}".format(
                qm_hessian.shape,
                mm_hessian.shape
            )
        )

    if qm_hessian.ndim != 2:

        raise ValueError(
            "RIC Hessian must be two-dimensional."
        )

    if (
        qm_hessian.shape[0]
        !=
        qm_hessian.shape[1]
    ):

        raise ValueError(
            "RIC Hessian is not square."
        )

    qm_norm = np.linalg.norm(
        qm_hessian,
        ord="fro"
    )

    if qm_norm < EPS:

        raise ValueError(
            "QM RIC Hessian norm is zero."
        )

    difference = (
        qm_hessian - mm_hessian
    )

    error = np.linalg.norm(
        difference,
        ord="fro"
    ) / qm_norm

    return float(error)


# ================================================================
# Normalization
# ================================================================

def minmax_normalize(values):
    """
    Normalize a list of errors to [0,1].
    """

    values = np.asarray(
        values,
        dtype=float
    )

    minimum = np.min(values)
    maximum = np.max(values)

    if abs(maximum - minimum) < EPS:

        return np.zeros_like(values)

    return (
        values - minimum
    ) / (
        maximum - minimum
    )


# ================================================================
# File pairing
# ================================================================

def find_qm_files(
        qm_dir,
        basename
):
    """
    QM files:

        <basename>.log
        <basename>.fchk
    """

    log_file = os.path.join(
        qm_dir,
        basename + ".log"
    )

    fchk_file = os.path.join(
        qm_dir,
        basename + ".fchk"
    )

    if not os.path.isfile(log_file):

        raise FileNotFoundError(
            "QM LOG not found: {}".format(
                log_file
            )
        )

    if not os.path.isfile(fchk_file):

        raise FileNotFoundError(
            "QM FCHK not found: {}".format(
                fchk_file
            )
        )

    return log_file, fchk_file


def find_hessfit_files(
        hessfit_dir,
        basename
):
    """
    HessFit files:

        <basename>_hessfit.log
        <basename>_hessfit.fchk
    """

    prefix = basename + "_hessfit"

    log_file = os.path.join(
        hessfit_dir,
        prefix + ".log"
    )

    fchk_file = os.path.join(
        hessfit_dir,
        prefix + ".fchk"
    )

    if not os.path.isfile(log_file):

        raise FileNotFoundError(
            "HessFit LOG not found: {}".format(
                log_file
            )
        )

    if not os.path.isfile(fchk_file):

        raise FileNotFoundError(
            "HessFit FCHK not found: {}".format(
                fchk_file
            )
        )

    return log_file, fchk_file


# ================================================================
# Evaluate one conformer
# ================================================================

def evaluate_conformer(
        basename,
        qm_dir,
        hessfit_dir,
        heavy_atoms=True
):
    """
    Evaluate one QM/HessFit pair.
    """

    print()
    print("-" * 72)
    print("Evaluating {}".format(basename))
    print("-" * 72)

    qm_log, qm_fchk = find_qm_files(
        qm_dir,
        basename
    )

    hf_log, hf_fchk = find_hessfit_files(
        hessfit_dir,
        basename
    )

    print("QM LOG       : {}".format(qm_log))
    print("QM FCHK      : {}".format(qm_fchk))
    print("HessFit LOG  : {}".format(hf_log))
    print("HessFit FCHK : {}".format(hf_fchk))

    # ------------------------------------------------------------
    # Gaussian termination
    # ------------------------------------------------------------

    if not check_gaussian_termination(
        qm_log
    ):

        raise RuntimeError(
            "QM calculation did not terminate normally: {}".format(
                qm_log
            )
        )

    if not check_gaussian_termination(
        hf_log
    ):

        raise RuntimeError(
            "HessFit calculation did not terminate normally: {}".format(
                hf_log
            )
        )

    # ------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------

    # qm_Z, qm_geometry = extract_final_geometry(
    #     qm_log
    # )

    qm_Z, qm_geometry = extract_fchk_geometry(
        qm_fchk
    )

    mm_Z, mm_geometry = extract_fchk_geometry(
        hf_fchk
    )

    print("\nGeometry comparison:")
    print("  QM atoms       :", len(qm_Z))
    print("  HessFit atoms  :", len(mm_Z))
    print("  QM coordinates :", qm_geometry.shape)
    print("  HF coordinates :", mm_geometry.shape)

    if len(qm_Z) != len(mm_Z):

        raise RuntimeError(
            "QM/HessFit atom count mismatch: {} vs {}".format(
                len(qm_Z),
                len(mm_Z)
            )
        )

    if not np.array_equal(qm_Z, mm_Z):

        raise RuntimeError(
            "QM and HessFit atomic-number ordering differ."
        )

    geometry_rmsd = calculate_geometry_rmsd(
        qm_Z,
        qm_geometry,
        mm_Z,
        mm_geometry,
        heavy_atoms=heavy_atoms
    )

    # ------------------------------------------------------------
    # Frequencies
    # ------------------------------------------------------------

    qm_frequencies = extract_final_frequencies(
        qm_log, len(qm_Z)
    )

    mm_frequencies = extract_frequencies(
        hf_log
    )

    print("\nVibrational frequencies:")
    print("  Expected modes :", 3 * len(qm_Z) - 6)
    print("  QM modes       :", len(qm_frequencies))
    print("  HessFit modes  :", len(mm_frequencies))

    frequency_data = calculate_frequency_error(
        qm_frequencies,
        mm_frequencies
    )

    # ------------------------------------------------------------
    # QM energy
    # ------------------------------------------------------------

    qm_energy = extract_scf_energy(
        qm_log
    )

    # ------------------------------------------------------------
    # RIC Hessians
    # ------------------------------------------------------------

    print("Reading QM RIC Hessian...")

    qm_hessian = load_ric_hessian(
        qm_fchk
    )

    print(
        "  QM RIC Hessian shape: {}".format(
            qm_hessian.shape
        )
    )

    print("Reading HessFit RIC Hessian...")

    mm_hessian = load_ric_hessian(
        hf_fchk
    )

    print(
        "  MM RIC Hessian shape: {}".format(
            mm_hessian.shape
        )
    )

    # ------------------------------------------------------------
    # Hessian error
    # ------------------------------------------------------------

    hessian_error = calculate_hessian_error(
        qm_hessian,
        mm_hessian
    )

    # ------------------------------------------------------------
    # Results
    # ------------------------------------------------------------

    return {
        "conformer": basename,
        "qm_energy": qm_energy,
        "hessian_error": hessian_error,
        "rmsd": geometry_rmsd,
        "frequency_error":
            frequency_data["frequency_error"],
        "n_modes":
            frequency_data["n_modes"],
        "n_imag_qm":
            frequency_data["n_imag_qm"],
        "n_imag_mm":
            frequency_data["n_imag_mm"],
        "n_ric":
            qm_hessian.shape[0]
    }


# ================================================================
# Main
# ================================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Rank selected conformers using "
            "RIC Hessian, geometry and frequency errors."
        )
    )

    parser.add_argument(
        "--conformers",
        nargs="+",
        required=True,
        help=(
            "QM basenames, e.g. "
            "02_conf_000 02_conf_010"
        )
    )

    parser.add_argument(
        "--qm-dir",
        default=".",
        help=(
            "Directory containing <basename>.log/.fchk"
        )
    )

    parser.add_argument(
        "--hessfit-dir",
        default=".",
        help=(
            "Directory containing "
            "<basename>_hessfit.log/.fchk"
        )
    )

    parser.add_argument(
        "--weights",
        nargs=3,
        type=float,
        default=[0.50, 0.25, 0.25],
        metavar=("WH", "WR", "WF"),
        help=(
            "Weights for Hessian, geometry and frequency "
            "errors. Default: 0.50 0.25 0.25"
        )
    )

    parser.add_argument(
        "--all-atoms",
        action="store_true",
        help=(
            "Use all atoms for geometry RMSD. "
            "Default: heavy atoms only."
        )
    )

    parser.add_argument(
        "--output",
        default="conformer_merit.csv",
        help=(
            "CSV output file."
        )
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Validate weights
    # ------------------------------------------------------------

    weights = np.asarray(
        args.weights,
        dtype=float
    )

    if np.any(weights < 0.0):

        parser.error(
            "Weights cannot be negative."
        )

    weight_sum = np.sum(weights)

    if weight_sum <= EPS:

        parser.error(
            "At least one weight must be non-zero."
        )

    weights /= weight_sum

    wH, wR, wF = weights

    # ------------------------------------------------------------
    # Header
    # ------------------------------------------------------------

    print()
    print("=" * 72)
    print("HESSFIT CONFORMER MERIT EVALUATION")
    print("=" * 72)

    print()
    print(
        "Hessian weight    : {:.3f}".format(wH)
    )

    print(
        "Geometry weight   : {:.3f}".format(wR)
    )

    print(
        "Frequency weight  : {:.3f}".format(wF)
    )

    print()
    print(
        "Selected conformers:"
    )

    for basename in args.conformers:

        print(
            "  {}".format(basename)
        )

    # ------------------------------------------------------------
    # Evaluate conformers
    # ------------------------------------------------------------

    results = []

    for basename in args.conformers:

        result = evaluate_conformer(
            basename=basename,
            qm_dir=args.qm_dir,
            hessfit_dir=args.hessfit_dir,
            heavy_atoms=not args.all_atoms
        )

        results.append(result)

    # ------------------------------------------------------------
    # Relative QM energies
    # ------------------------------------------------------------

    energies = np.asarray([
        result["qm_energy"]
        for result in results
    ])

    finite = np.isfinite(
        energies
    )

    if np.any(finite):

        minimum_energy = np.min(
            energies[finite]
        )

        for result in results:

            if np.isfinite(
                result["qm_energy"]
            ):

                result["delta_energy"] = (
                    result["qm_energy"]
                    -
                    minimum_energy
                ) * HARTREE_TO_KCAL

            else:

                result["delta_energy"] = np.nan

    else:

        for result in results:

            result["delta_energy"] = np.nan

    # ------------------------------------------------------------
    # Raw errors
    # ------------------------------------------------------------

    hessian_errors = np.asarray([
        result["hessian_error"]
        for result in results
    ])

    rmsd_errors = np.asarray([
        result["rmsd"]
        for result in results
    ])

    frequency_errors = np.asarray([
        result["frequency_error"]
        for result in results
    ])

    # ------------------------------------------------------------
    # Normalize
    # ------------------------------------------------------------

    hessian_normalized = minmax_normalize(
        hessian_errors
    )

    rmsd_normalized = minmax_normalize(
        rmsd_errors
    )

    frequency_normalized = minmax_normalize(
        frequency_errors
    )

    # ------------------------------------------------------------
    # Merit
    # ------------------------------------------------------------

    for i, result in enumerate(results):

        result["hessian_normalized"] = (
            hessian_normalized[i]
        )

        result["rmsd_normalized"] = (
            rmsd_normalized[i]
        )

        result["frequency_normalized"] = (
            frequency_normalized[i]
        )

        result["merit"] = (
            wH * hessian_normalized[i]
            +
            wR * rmsd_normalized[i]
            +
            wF * frequency_normalized[i]
        )

    # ------------------------------------------------------------
    # Rank
    # ------------------------------------------------------------

    results.sort(
        key=lambda x: x["merit"]
    )

    # ------------------------------------------------------------
    # Print table
    # ------------------------------------------------------------

    print()
    print("=" * 100)
    print("CONFORMER RANKING")
    print("=" * 100)

    print()

    print(
        "{:<18s} {:>9s} {:>12s} {:>10s} "
        "{:>12s} {:>10s}".format(
            "Conformer",
            "ΔE",
            "Hessian",
            "RMSD",
            "Frequency",
            "Merit"
        )
    )

    print("-" * 100)

    for result in results:

        print(
            "{:<18s} {:>9.3f} {:>12.6f} "
            "{:>10.5f} {:>12.6f} {:>10.6f}".format(
                result["conformer"],
                result["delta_energy"],
                result["hessian_error"],
                result["rmsd"],
                result["frequency_error"],
                result["merit"]
            )
        )

    # ------------------------------------------------------------
    # Components
    # ------------------------------------------------------------

    print()
    print("Normalized merit components")
    print("-" * 72)

    print(
        "{:<18s} {:>12s} {:>12s} {:>12s} {:>12s}".format(
            "Conformer",
            "Hessian",
            "Geometry",
            "Frequency",
            "Merit"
        )
    )

    print("-" * 72)

    for result in results:

        print(
            "{:<18s} {:>12.5f} {:>12.5f} "
            "{:>12.5f} {:>12.5f}".format(
                result["conformer"],
                result["hessian_normalized"],
                result["rmsd_normalized"],
                result["frequency_normalized"],
                result["merit"]
            )
        )

    # ------------------------------------------------------------
    # Frequency diagnostics
    # ------------------------------------------------------------

    print()
    print("Frequency diagnostics")
    print("-" * 72)

    for result in results:

        print(
            "{} : modes={}  "
            "imag(QM)={}  "
            "imag(HessFit)={}".format(
                result["conformer"],
                result["n_modes"],
                result["n_imag_qm"],
                result["n_imag_mm"]
            )
        )

    # ------------------------------------------------------------
    # Winner
    # ------------------------------------------------------------

    winner = results[0]

    print()
    print("=" * 72)
    print(
        "OPTIMAL CONFORMER: {}".format(
            winner["conformer"]
        )
    )

    print(
        "Merit score: {:.6f}".format(
            winner["merit"]
        )
    )

    print("=" * 72)

    # ------------------------------------------------------------
    # CSV
    # ------------------------------------------------------------

    fieldnames = [
        "rank",
        "conformer",
        "qm_energy",
        "delta_energy",
        "hessian_error",
        "rmsd",
        "frequency_error",
        "hessian_normalized",
        "rmsd_normalized",
        "frequency_normalized",
        "merit",
        "n_ric",
        "n_modes",
        "n_imag_qm",
        "n_imag_mm"
    ]

    with open(
        args.output,
        "w",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames
        )

        writer.writeheader()

        for rank, result in enumerate(
            results,
            start=1
        ):

            row = {
                "rank": rank
            }

            for field in fieldnames:

                if field == "rank":
                    continue

                row[field] = result.get(
                    field,
                    ""
                )

            writer.writerow(row)

    print()
    print(
        "Results written to {}".format(
            args.output
        )
    )


if __name__ == "__main__":
    main()

# qm_Z, qm_geometry = extract_final_geometry(
#     "../examples/conformers/qm_confs/01_conf_000.log"
# )

# print("Number of atoms:", len(qm_Z))
# print("Geometry shape:", qm_geometry.shape)

# print("\nAtomic numbers:")
# print(qm_Z)

# print("\nFirst 5 coordinates:")
# print(qm_geometry[:5])

# print("\nLast 5 coordinates:")
# print(qm_geometry[-5:])