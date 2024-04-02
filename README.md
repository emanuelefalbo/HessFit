# HessFit

hessfit returns the entire harmonic force fields and performs diehdral force constant fitting for any molecule. 
In the directory "src", there is the **hessfit_harmonic.py** which returns calculated harmonic intra-molecular force field, i.e, stretching, bending, and dihedral force constants, and external inter-molecular parmateres, i.e, Van der Waals and atomi charges. As such, these latter might result too stiff for specific applications. 
As alternative, least square fit of quantum-mechanically PES can be carried out with **hessfit_dihes.py**. 
The directory utils contains instead all the necessary modules.

# Install

It currently works with the Gaussian outputs, and json files containing the necessary input information. It has been currently tested with Gaussian09. 
A numpy python library is required to properly use the program alongside other module.
Numpy can be installed by the following command:
```
pip install numpy pandas scipy
```

Then, by cloning the directory locally, user should give:
```
git clone https://github.com/emanuelefalbo/HessFit
cd HessFit
python setup.py install
```

Otherwise, it is sufficient to add their directory to the **PYTHONPATH** into you .bashrc (bash) file:
```
export PYTHONPATH="${PYTHONPATH}:/path/to/hessfit/src"
export PATH=$PATH:/path/to/hessfit/src
```
with path/to/hessfit being the full path to where it is located. 
The last line add the programs to your bash path (see setenv for .tcsh) to make it visibile anywhere.

# Usage 

**1.1 Harmonic Force Field***

The program is thought of as dual usage, i.e., it can be executed by launching the **hessfit.py** script that performs all 
operations provided by the *build_4_hessfit.py* first, and *hessfit_harmonic.py* secondly, or these two scripts can be run independently.
It works with a json file as described below:

```
hessfit.py  step1.json -path $g09root
```

The step1.json is composed as :

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
    "charge":0,
    "multiplicity":1,
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
 "charge": the molecular charge of compound
 "multiplicity": the molecular multiplicity accoding to spin state
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

If atom types are absent in the second column, a full list of atom types with their serial index will be used.
GauHarm.gjf and GauHarm.gjf are the necessary input files for *hessfit_harmonic.py* and need to be specified in the step1.json, although these are created during the run.
The output is a *hessfit4gau.gjf* and its output (hessfit4gau.log) which contains the optimization+frequency of the system with computed HessFit harmonic force field. 
A further *hessfit_frcmod.txt*, which is in AMBER-like format,  is outputted and can be transfered for other Moleculaar Mechanics software such as AMBER or GROMACS.


**1.2. Dihedral Parameters**

Further accuracy on dihedral can be achieved by using a linear least-square fitting of torsional angles related to rotational bonds.
After the first step, a directory called **dihedral** is created. Users should enter in this directory, and run: 


```
hessfit_dihes.py  step2.json -path $g09root
```

The step2.json is composed as :

```json
{
    "files": {
        "atom2type": "type_charge.txt",
        "force_file": "ff_string.txt",
        "file_xyz": ".qm_opt.xyz",
        "topol": "topol.txt"
    },
    "method": "B3LYP/6-31G*",
    "nprocs": "4"
    }
```


```
 "atom2type":  Input file including atom type and atomic charges
 "force_file": Input file containg string text 
 "topol": Input file containing topolgy information 
 "file_xyz": QM-optimized geometry in xyz format 
 "method": Computational method for running the scan (e.g. B3LYP/6-31G*)
 "nprocs": Number of processor to be used. 

```

While type_charge.txt, ff_string.txt, and topol.txt are generated internally, users must specify which optimized xyz geometry to use for the torsional scans.
It must be noted that since several methods, like ab initio or semi-empirical ones, can be chosen for the torsional scan, their completion time will depend on the number of processors and method chosen. The output files *x_qm_all.csv* (x=1,2,...) contains the angle, QM, and MM energies in 1st, 2nd, and 3rd columns, respectively for each scanned dihedral. The internal subroutine *fit4dihe.py* processes these files and return the Fourier coefficients. 



