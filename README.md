# PySmartField

SmartField returns the entire force field for any molecule.
In the repository "smart_harmonic", there is the **smart_bond_angles.py** program which returns the harmonic intra-molecular force field, i.e, stretching, bending, and dihedral force constants. As such, these latter might result too stiff for specific applications. 
As alternative, least square fit of quantum-mechanically PES can be carried out with **smart_dihedral.  

# Usage 

It currently works with the Gaussian outputs, and a json file containing the necessary files, and type of results.
A numpy python library is required to properly use the program alongside other module.
Numpy can be installed by the following command:
```
pip install numpy
```
while for the modules it is sufficient to add their directory to the **PYTHONPATH** into you .bashrc (bash) file:
```
export PYTHONPATH="${PYTHONPATH}:/path/to/modules
export PATH=$PATH:/path/to/PySmartField/smart_harmonic
```
The latter line add the main program **smart_bond_angle.py to your path (see setenv for .tcsh).

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

```
 "log_qm_file":  Gaussian output log file of QM (opt+Freq) calculation
 "fchk_qm_file": Gaussian format check-point file of QM (opt+freq) calculation
 "fchk_mm_file": Gaussian format check-point file of artificial MM (freq) calculation
 "fchk_nb_file": Gaussian format check-point file of artificial MM (freq) calculation
```



**2. smart_harmonic**



