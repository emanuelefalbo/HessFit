# SmartField

SmartField returns the entire force field for any molecule.
In the directory "src", there is the **SmartField_harmonic.py** which returns calculated harmonic intra-molecular force field, i.e, stretching, bending, and dihedral force constants, and external inter-molecular parmateres, i.e, Van der Waals and atomi charges. As such, these latter might result too stiff for specific applications. 
As alternative, least square fit of quantum-mechanically PES can be carried out with **SmartField_torsion.py**. 
The directoryutils contains instead all the necessary modules.

# Install 

It currently works with the Gaussian outputs, and a json file containing the necessary files, and type of results.
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
with path/to/PySmartField being the full path to where it is located. 
The last two lines add the programs to your bash path (see setenv for .tcsh).

# Usage 


**1.0 build_4Smart**

The program in **SmartField_harmonic.py** needs some mandatory input files that can be esily created 
with program *build_4Smart.py** which is located in the directory with same name.
By calling it from command line as:
```
build_4Smart.py -f1 file.log -f2 file.fchk -m mean -path $g09root
```
it returns two files **GauHarm.gjf** and **GauHarm.gjf**, which are the 
necessary input files for **SmartField_harmonic.py**. Info about the positional argument can be obtained as follows:
```
usage: build_4Smart.py [-h] [-f1 LOG_FILE] [-f2 FCHK_FILE] [-m {all,mea

optional arguments:
  -h, --help            show this help message and exit
  -f1 LOG_FILE, --log_file LOG_FILE
                        Gaussian QM log file 
  -f2 FCHK_FILE, --fchk_file FCHK_FILE
                        Gaussain QM fchk file
  -m {all,mean}, --mode {all,mean}
                        averaging across same types
  -path PATH            path/to/amber.prm in Gaussain root directory
```

**1.1 Harmonic FF**

Provided that all modules and main program of **PySmartField** are in your path,
it can be simply used as:
```
SmartField_harmonic.py.py option_file.json
```
This part of program works along with a json file like the following:

```json
{
    "files": {
        "log_qm_file":  "but_qm.log",
        "fchk_qm_file": "but_qm.fchk",
        "fchk_mm_file": "GauHarm.fchk",
        "fchk_nb_file": "GauNonBon.fchk",
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



**1.2. Dihedral Parameters**



