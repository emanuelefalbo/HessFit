#!/usr/bin/env python3

import subprocess
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Execute Gaussian and SmartField scripts.")
    parser.add_argument("--f1", help="Path to F1 file ")
    parser.add_argument("--f2", help="Path to F2 file ")
    parser.add_argument("--json", default="opt.json", help="Path to JSON file ")
    return parser.parse_args()

def main():
    args = parse_args()

    F1 = args.f1
    F2 = args.f2
    JSON = args.json

    GPATH = os.environ.get("g09root") + "/g09"
    SM = "Smart_harmonic.py"
    BS = "build_4Smart.py"

    subprocess.run([BS, "-f1", F1, "-f2", F2, "-path", GPATH], check=True)

    for f in ["GauHarm.gjf", "GauNonBon.gjf"]:
        print(f"Executing Gaussian on {f}")
        subprocess.run([f"{GPATH}/g09", f], check=True)

    for f in ["GauHarm.chk", "GauNonBon.chk"]:
        print(f"Formatchecking {f} file")
        subprocess.run([f"{GPATH}/formchk", "-3", f, f"{os.path.splitext(f)[0]}.fchk"], check=True)

    print(f"Executing {SM}")
    subprocess.run(["SmartField_harmonic.py", JSON])
    subprocess.run([f"{GPATH}/g09", "SmartField4gau.gjf"], check=True)

if __name__ == "__main__":
    main()
