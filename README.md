# HessFit

SmartField returns the entire harmonic force fields and performs diehdral force constant fitting for any molecule. 
In the directory "src", there is the **SmartField_harmonic.py** which returns calculated harmonic intra-molecular force field, i.e, stretching, bending, and dihedral force constants, and external inter-molecular parmateres, i.e, Van der Waals and atomi charges. As such, these latter might result too stiff for specific applications. 
As alternative, least square fit of quantum-mechanically PES can be carried out with **SmartField_dihes.py**. 
The directory utils contains instead all the necessary modules.

# Install 

It currently works with the Gaussian outputs, and json files containing the necessary input information. It has been currently tested with Gaussian09. 
A numpy python library is required to properly use the program alongside other module.
Numpy can be installed by the following command:
```
pip install numpy
```
while for the modules it is sufficient to add their directory to the **PYTHONPATH** into you .bashrc (bash) file:
```
export PYTHONPATH="${PYTHONPATH}:/path/to/SmartField/src"
export PYTHONPATH="${PYTHONPATH}:/path/to/SmartField/utils"
export PYTHONPATH="${PYTHONPATH}:/path/to/SmartField/build_4Smart"
export PATH=$PATH:/path/to/SmartField/src
export PATH=$PATH:/path/to/SmartField/utils
export PATH=$PATH:/path/to/SmartField/build_4Smart
```
with path/to/SmartField being the full path to where it is located. 
The last two lines add the programs to your bash path (see setenv for .tcsh).

# Usage 

**1.1 Harmonic Force Field***

The program is thought of as dual usage, i.e., it can be executed by launching the **run_smartfield.py** script that performs all 
operations provided by the *build_4Smart.py* first, and *SmartField_harmonic.py* secondly, or these two scripts can be run independently
It works with a json file as described below:

```
run_smartfield.py  file1.json -path $g09root
```

The file1.json is composed as :

```json
{
    "files": {
        "log_qm_file":  "file.log",
        "fchk_qm_file": "file.fchk",
        "fchk_mm_file": "GauHarm.fchk",
        "fchk_nb_file": "GauNonBon.fchk",
        "atype_file": "file_type.txt"
         },
    "mode": "mean",
    "opt": "ric"
}
```

```
 "log_qm_file":  Gaussian output log file of QM (opt+Freq) calculation
 "fchk_qm_file": Gaussian format check-point file of QM (opt+freq) calculation
 "fchk_mm_file": Gaussian format check-point file of artificial MM (freq) calculation
 "fchk_nb_file": Gaussian format check-point file of artificial MM (freq) calculation
 "atype_file": column file containing element and atom types
 "mode": "mean" string averages all force contsants over same types, while "all" leaves it unchanged
 "opt": "ric" string perform the Hessian diagonalization in redundant internal coordinates, whereas "sem" string performs the Seminario method. 
```

The "atype_file" must be a two-column file with the element, and the atom type on the first and second column respectively: 

```
N-N3
C-CT
C-C 
O-O 
C-CT
C-CA
...
```


GauHarm.gjf and GauHarm.gjf are the necessary input files for *SmartField_harmonic.py* and need to be specified in the file1.json.
The output is a *SmartField4gau.gjf* and its output (SmartField4gau.log) which contains the optimization of the system with computed harmonic force field. 
A further *SmartField_frcmod.txt* is outputted and can be transfered for other Moleculaar Mechanics software such as AMBER or GROMACS.


**1.2. Dihedral Parameters**

...


