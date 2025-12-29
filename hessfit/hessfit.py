#!/usr/bin/env python3

import subprocess
import os
import readin_opts as rdin
import parser_gau as pgau
import geom2atype as g2a

def print_init():
     print("""
 ======================================================
   Program:      HessFit
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
    parser.add_argument('--version', choices=['g09', 'g16'], default='g09', help='Select Gaussian version (g09 or g16)')
    parser.add_argument('--at', choices=["gaff", "amber"], default="gaff",
                        help='Select force field for VDW parameters (gaff or amber); default = gaff')
    parser.add_argument('--test', default=False, help='Test HessFit FF with Gaussian internal routines')
    opts = parser.parse_args()
    print(opts.at)
    
    g09root = os.environ.get('g09root')
    g16root = os.environ.get('g16root')
    
    # Determine which Gaussian version to use
    if opts.version == 'g16':
        groot = g16root if g16root else None
    else:
        groot = g09root if g09root else None
    
    # Validate the selected Gaussian path
    if groot:
        GPATH = os.path.join(groot, opts.version)
        print(f"GPATH == {GPATH}")
        if not os.path.exists(GPATH):
            print(f"Warning: {GPATH} does not exist. Falling back to user-defined path.")
            groot = opts.path
            GPATH = os.path.join(groot, opts.version)
    else:
        groot = opts.path
        GPATH = os.path.join(groot, opts.version)
    
    if not os.path.exists(GPATH):
        raise FileNotFoundError(f"Error: No valid Gaussian installation found at {GPATH}.")
    

    SM = "hessfit_harmonic.py"
    BS = "build_4_hessfit.py"
    JSON = opts.optfile


    print(f"Executing Building .gjf files: {BS}")
    subprocess.run([BS, JSON, "--path", GPATH, "--at", opts.at], check=True)
    
    gaussian_exe = "g16" if opts.version == "g16" else "g09"
    for f in ["GauHarm.gjf", "GauNonBon.gjf"]:
        print(f"Executing Gaussian on {f}")
        subprocess.run([f"{GPATH}/{gaussian_exe}", f], check=True)

    for f in ["GauHarm.chk", "GauNonBon.chk"]:
        print(f"Formatchecking {f} file")
        subprocess.run([f"{GPATH}/formchk", "-3", f, f"{os.path.splitext(f)[0]}.fchk"], check=True)

    print(f"Executing Harmonic: {SM}")
    subprocess.run([SM, JSON, "--version", opts.version, "--at", opts.at], check=True)

    if opts.test:
        print(f"Executing Gaussian on hessfit4gau.gjf")
        subprocess.run([f"{GPATH}/{gaussian_exe}", "hessfit4gau.gjf"], check=True)

if __name__ == "__main__":
    main()
