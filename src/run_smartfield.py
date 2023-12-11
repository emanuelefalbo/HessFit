#!/usr/bin/env python3

import subprocess
import os
import readin_opts as rdin


def print_init():
     print("""
 ======================================================
   Program:      SmartField
   Creator:      Emanuele Falbo, Napoli
   Language:     Python 3.v later
   Description:  The program returns force constants for
                 bonded values and non-bonded parameters
   Mail:         falbo.emanuele@gmail.com 
 =======================================================
     
    """)


def main():
    print_init()
    parser = rdin.commandline_parser3()
    opts = parser.parse_args()

    GPATH = os.environ.get("g09root") + "/g09"
    SM = "Smart_harmonic.py"
    BS = "build_4Smart.py"
    JSON = opts.optfile

    # subprocess.run([BS, JSON, "-path", GPATH], check=True)
    subprocess.run([BS, JSON], check=True)

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
