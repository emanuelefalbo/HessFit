#!/usr/bin/env python3

import os
import glob
import json
import argparse

from xtb_screening import run_xtb


OUTPUT_DIR = "xtb_screening"

CHARGE = 0
UHF = 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default="conformers", help="Directory with conformer XYZ files (default: conformers)")
    args = parser.parse_args()
    
    INPUT_DIR = args.input_dir

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    xyz_files = sorted(
        glob.glob(
            os.path.join(INPUT_DIR, "conf_*.xyz")
        )
    )

    print(f"Found {len(xyz_files)} candidate conformers.")

    results = []

    for i, xyz_file in enumerate(xyz_files):

        conformer_id = os.path.splitext(
            os.path.basename(xyz_file)
        )[0]

        directory = os.path.join(
            OUTPUT_DIR,
            conformer_id
        )

        print(
            f"\n[{i + 1}/{len(xyz_files)}] "
            f"{conformer_id}"
        )

        result = run_xtb(
            xyz_file,
            directory,
            charge=CHARGE,
            uhf=UHF
        )

        data = {
            "conformer": conformer_id,
            "input_xyz": xyz_file,
            "optimized_xyz": result.optimized_xyz,
            "energy_hartree": result.energy,
            "success": result.success,
            "directory": directory
        }

        results.append(data)

        if result.success:
            print(
                f"  Energy = "
                f"{result.energy:.10f} Eh"
            )
        else:
            print("  xTB optimization FAILED")

    with open("xtb_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\nResults written to xtb_results.json")


if __name__ == "__main__":
    main()