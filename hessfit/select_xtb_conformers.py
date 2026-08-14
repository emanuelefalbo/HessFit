#!/usr/bin/env python3

"""
select_xtb_conformers.py

Select representative conformers after xTB optimization.

Workflow
--------
1. Read xTB results from xtb_results.json
2. Remove failed calculations
3. Calculate relative xTB energies
4. Apply an energy cutoff
5. Calculate heavy-atom RMSD using Kabsch alignment
6. Cluster conformers using an RMSD threshold
7. Keep the lowest-energy member of each cluster
8. Select up to N representative conformers
9. Copy selected optimized XYZ files to selected_conformers/

The selected structures are intended to be used as starting
geometries for subsequent high-level Gaussian Opt + Freq jobs.

Important
---------
- RMSD is calculated using heavy atoms only.
- Atom ordering is assumed to be identical between conformers.
- The original Gaussian RIC definition is NOT modified.
"""

import os
import json
import shutil
import argparse

import numpy as np


HARTREE_TO_KCAL = 627.509474


# ================================================================
# XYZ READER
# ================================================================

def read_xyz(filename):
    """
    Read an XYZ file.

    Returns
    -------
    symbols : list[str]
    coordinates : np.ndarray, shape (N, 3)
    """

    with open(filename) as f:
        lines = f.readlines()

    if len(lines) < 2:
        raise ValueError(
            f"Invalid XYZ file: {filename}"
        )

    natoms = int(lines[0].strip())

    symbols = []
    coordinates = []

    for line in lines[2:2 + natoms]:

        fields = line.split()

        if len(fields) < 4:
            raise ValueError(
                f"Malformed XYZ line in {filename}: {line}"
            )

        symbols.append(fields[0])

        coordinates.append([
            float(fields[1]),
            float(fields[2]),
            float(fields[3])
        ])

    if len(symbols) != natoms:
        raise ValueError(
            f"Expected {natoms} atoms but found "
            f"{len(symbols)} in {filename}"
        )

    return (
        symbols,
        np.asarray(coordinates, dtype=float)
    )


# ================================================================
# HEAVY ATOM MASK
# ================================================================

def heavy_atom_mask(symbols):
    """
    Return a boolean mask selecting non-hydrogen atoms.
    """

    return np.array(
        [
            symbol.upper() != "H"
            for symbol in symbols
        ],
        dtype=bool
    )


# ================================================================
# KABSCH ALIGNMENT
# ================================================================

def kabsch_rmsd(coords_a, coords_b):
    """
    Calculate RMSD after optimal rotational alignment.

    Parameters
    ----------
    coords_a, coords_b
        Arrays of shape (N, 3).

    Returns
    -------
    rmsd : float
        RMSD in Angstrom.
    """

    if coords_a.shape != coords_b.shape:
        raise ValueError(
            "Coordinate arrays have different shapes."
        )

    # Remove translations
    centroid_a = np.mean(
        coords_a,
        axis=0
    )

    centroid_b = np.mean(
        coords_b,
        axis=0
    )

    a = coords_a - centroid_a
    b = coords_b - centroid_b

    # Covariance matrix
    covariance = np.dot(
        a.T,
        b
    )

    # Singular value decomposition
    v, s, wt = np.linalg.svd(
        covariance
    )

    # Prevent reflection
    determinant = (
        np.linalg.det(v)
        * np.linalg.det(wt)
    )

    correction = np.eye(3)

    if determinant < 0.0:
        correction[-1, -1] = -1.0

    rotation = np.dot(
        v,
        np.dot(
            correction,
            wt
        )
    )

    # Rotate A
    a_rot = np.dot(
        a,
        rotation
    )

    difference = a_rot - b

    rmsd = np.sqrt(
        np.mean(
            np.sum(
                difference ** 2,
                axis=1
            )
        )
    )

    return float(rmsd)


# ================================================================
# CONFORMER OBJECT
# ================================================================

class Conformer:

    def __init__(
        self,
        conformer_id,
        energy_hartree,
        optimized_xyz,
        directory,
    ):

        self.id = conformer_id
        self.energy_hartree = energy_hartree
        self.optimized_xyz = optimized_xyz
        self.directory = directory

        self.relative_energy = None

        self.symbols = None
        self.coordinates = None
        self.heavy_coordinates = None

        self.cluster = None

        self.rmsd_to_representative = None

    def load_geometry(self):

        self.symbols, self.coordinates = read_xyz(
            self.optimized_xyz
        )

        mask = heavy_atom_mask(
            self.symbols
        )

        self.heavy_coordinates = (
            self.coordinates[mask]
        )


# ================================================================
# LOAD RESULTS
# ================================================================

def load_xtb_results(filename):

    with open(filename) as f:
        data = json.load(f)

    conformers = []

    for item in data:

        if not item.get("success", False):
            print(
                f"Skipping failed calculation: "
                f"{item['conformer']}"
            )
            continue

        xyz = item.get(
            "optimized_xyz"
        )

        energy = item.get(
            "energy_hartree"
        )

        if xyz is None:
            print(
                f"Skipping {item['conformer']}: "
                "no optimized XYZ."
            )
            continue

        if energy is None:
            print(
                f"Skipping {item['conformer']}: "
                "no xTB energy."
            )
            continue

        if not os.path.isfile(xyz):
            print(
                f"Skipping {item['conformer']}: "
                f"XYZ not found: {xyz}"
            )
            continue

        conformer = Conformer(
            conformer_id=item["conformer"],
            energy_hartree=float(energy),
            optimized_xyz=xyz,
            directory=item.get("directory")
        )

        conformer.load_geometry()

        conformers.append(conformer)

    return conformers


# ================================================================
# RELATIVE ENERGY
# ================================================================

def calculate_relative_energies(conformers):

    if not conformers:
        return

    minimum = min(
        c.energy_hartree
        for c in conformers
    )

    for c in conformers:

        c.relative_energy = (
            c.energy_hartree - minimum
        ) * HARTREE_TO_KCAL


# ================================================================
# ENERGY FILTER
# ================================================================

def energy_filter(
    conformers,
    energy_cutoff
):

    selected = [
        c
        for c in conformers
        if c.relative_energy <= energy_cutoff
    ]

    return sorted(
        selected,
        key=lambda c: c.relative_energy
    )


# ================================================================
# RMSD CLUSTERING
# ================================================================

def cluster_by_rmsd(
    conformers,
    rmsd_cutoff
):
    """
    Cluster conformers by heavy-atom RMSD.

    Conformers are processed in ascending energy order.

    The lowest-energy conformer becomes the representative
    of each new cluster.

    Returns
    -------
    clusters : list[list[Conformer]]
    """

    conformers = sorted(
        conformers,
        key=lambda c: c.relative_energy
    )

    clusters = []

    for conformer in conformers:

        assigned = False

        for cluster_index, cluster in enumerate(clusters):

            representative = cluster[0]

            rmsd = kabsch_rmsd(
                conformer.heavy_coordinates,
                representative.heavy_coordinates
            )

            if rmsd < rmsd_cutoff:

                cluster.append(
                    conformer
                )

                conformer.cluster = (
                    cluster_index + 1
                )

                conformer.rmsd_to_representative = rmsd

                assigned = True

                break

        if not assigned:

            cluster_index = len(clusters)

            conformer.cluster = (
                cluster_index + 1
            )

            conformer.rmsd_to_representative = 0.0

            clusters.append(
                [conformer]
            )

    return clusters


# ================================================================
# SELECT REPRESENTATIVES
# ================================================================

def select_representatives(
    clusters,
    max_conformers
):
    """
    Select the lowest-energy conformer from each cluster.

    Clusters are ranked according to the energy of their
    representative.
    """

    representatives = []

    for cluster in clusters:

        representative = min(
            cluster,
            key=lambda c: c.relative_energy
        )

        representatives.append(
            representative
        )

    representatives.sort(
        key=lambda c: c.relative_energy
    )

    return representatives[:max_conformers]


# ================================================================
# COPY SELECTED STRUCTURES
# ================================================================

def copy_selected(
    selected,
    output_dir
):

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    for rank, conformer in enumerate(
        selected,
        start=1
    ):

        destination = os.path.join(
            output_dir,
            f"{rank:02d}_{conformer.id}.xyz"
        )

        shutil.copy2(
            conformer.optimized_xyz,
            destination
        )


# ================================================================
# WRITE JSON
# ================================================================

def write_selection_json(
    selected,
    filename
):

    output = []

    for rank, c in enumerate(
        selected,
        start=1
    ):

        output.append({
            "rank": rank,
            "conformer": c.id,
            "cluster": c.cluster,
            "energy_hartree": c.energy_hartree,
            "relative_energy_kcal_mol": c.relative_energy,
            "rmsd_to_cluster_representative": (
                c.rmsd_to_representative
            ),
            "optimized_xyz": c.optimized_xyz
        })

    with open(filename, "w") as f:

        json.dump(
            output,
            f,
            indent=2
        )


# ================================================================
# PRINT SUMMARY
# ================================================================

def print_summary(
    all_conformers,
    energy_selected,
    clusters,
    selected
):

    print("\n" + "=" * 72)
    print("xTB CONFORMER SCREENING SUMMARY")
    print("=" * 72)

    print(
        f"Successful xTB calculations : "
        f"{len(all_conformers)}"
    )

    print(
        f"Within energy cutoff         : "
        f"{len(energy_selected)}"
    )

    print(
        f"Unique RMSD clusters         : "
        f"{len(clusters)}"
    )

    print(
        f"Final conformers selected    : "
        f"{len(selected)}"
    )

    print("\nAll energy-filtered conformers:")
    print(
        f"{'ID':<12}"
        f"{'Energy (Eh)':>18}"
        f"{'ΔE (kcal/mol)':>18}"
        f"{'Cluster':>10}"
        f"{'RMSD (Å)':>12}"
    )

    print("-" * 72)

    for c in sorted(
        energy_selected,
        key=lambda x: x.relative_energy
    ):

        print(
            f"{c.id:<12}"
            f"{c.energy_hartree:>18.10f}"
            f"{c.relative_energy:>18.3f}"
            f"{c.cluster:>10}"
            f"{c.rmsd_to_representative:>12.3f}"
        )

    print("\nSelected conformers for high-level QM:")
    print(
        f"{'Rank':<8}"
        f"{'ID':<12}"
        f"{'ΔE (kcal/mol)':>18}"
        f"{'Cluster':>10}"
    )

    print("-" * 52)

    for rank, c in enumerate(
        selected,
        start=1
    ):

        print(
            f"{rank:<8}"
            f"{c.id:<12}"
            f"{c.relative_energy:>18.3f}"
            f"{c.cluster:>10}"
        )

    print("=" * 72)


# ================================================================
# MAIN
# ================================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Filter and select representative "
            "xTB-optimized conformers."
        )
    )

    parser.add_argument(
        "--results",
        default="xtb_results.json",
        help="xTB results JSON file."
    )

    parser.add_argument(
        "--energy-cutoff",
        type=float,
        default=5.0,
        help=(
            "Maximum relative xTB energy in "
            "kcal/mol. Default: 5.0"
        )
    )

    parser.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=0.5,
        help=(
            "Heavy-atom RMSD clustering cutoff "
            "in Angstrom. Default: 0.5"
        )
    )

    parser.add_argument(
        "--max-conformers",
        type=int,
        default=5,
        help=(
            "Maximum number of final conformers. "
            "Default: 5"
        )
    )

    parser.add_argument(
        "--output-dir",
        default="selected_conformers",
        help=(
            "Directory for selected XYZ structures."
        )
    )

    parser.add_argument(
        "--json-output",
        default="selected_xtb_conformers.json",
        help=(
            "JSON file containing final selection."
        )
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Load
    # ------------------------------------------------------------

    conformers = load_xtb_results(
        args.results
    )

    if not conformers:

        raise RuntimeError(
            "No successful xTB conformers found."
        )

    # ------------------------------------------------------------
    # Relative energies
    # ------------------------------------------------------------

    calculate_relative_energies(
        conformers
    )

    # ------------------------------------------------------------
    # Energy filter
    # ------------------------------------------------------------

    energy_selected = energy_filter(
        conformers,
        args.energy_cutoff
    )

    if not energy_selected:

        raise RuntimeError(
            "No conformers survived the energy cutoff."
        )

    # ------------------------------------------------------------
    # RMSD clustering
    # ------------------------------------------------------------

    clusters = cluster_by_rmsd(
        energy_selected,
        args.rmsd_cutoff
    )

    # ------------------------------------------------------------
    # Select representatives
    # ------------------------------------------------------------

    selected = select_representatives(
        clusters,
        args.max_conformers
    )

    # ------------------------------------------------------------
    # Copy XYZ files
    # ------------------------------------------------------------

    copy_selected(
        selected,
        args.output_dir
    )

    # ------------------------------------------------------------
    # JSON
    # ------------------------------------------------------------

    write_selection_json(
        selected,
        args.json_output
    )

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------

    print_summary(
        conformers,
        energy_selected,
        clusters,
        selected
    )


if __name__ == "__main__":
    main()
