# PySmartField

SmartField returns the entire force field for any molecule.
In the repository "smart_harmonic", there is the **smart_bond_angles.py** program which returns the harmonic intra-molecular force field, i.e, stretching, bending, and dihedral force constants. As such, these latter might result too stiff for specific applications. 
As alternative, least square fit of quantum-mechanically PES can be carried out with **smart_dihedral.  

# Usage 

It currently works with the Gaussian outputs, and a json file containing the necessary files, and type of results.
A numpy python library is required to properly use the program.
It can be installed by the following command:
```
pip install numpy
```
**1. smart_harmonic**

Provided that all modules and main program of **PySmartField** are in your path,
it can be simply used as:
```
smart_bonds_angles.py option_file.json**
```
This part of program works along with a json file like the following:

```json
{
    "files": {
        "log_qm_file":  "but_qm.log",
        "fchk_qm_file": "but_qm.fchk",
        "fchk_mm_file": "mm_all.fchk",
        "fchk_nb_file": "nb_ric.fchk",
        "nonbon": "nonbonded.txt"
         },
    "mode": "all",
    "opt": "ric"
}
```

|   |   |   |   |   |
|---|---|---|---|---|
|   |   |   |   |   |
|   |   |   |   |   |
|   |   |   |   |   |



**2. smart_harmonic**



