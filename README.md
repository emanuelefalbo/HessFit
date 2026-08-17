# HessFit

hessfit returns the entire harmonic force fields and performs diehdral force constant fitting for any molecule. 
In the directory "src", there is the **hessfit_harmonic.py** which returns calculated harmonic intra-molecular force field, i.e, stretching, bending, and dihedral force constants, and external inter-molecular parmateres, i.e, Van der Waals and atomi charges. As such, these latter might result too stiff for specific applications. 
As alternative, least square fit of quantum-mechanically PES can be carried out with **hessfit_dihes.py**. 
The directory utils contains instead all the necessary modules.

# Install

It currently works with the Gaussian outputs, and json files containing the necessary input information. It has been currently tested with Gaussian09/16. 
The following python libraries are srequired to properly use the program alongside other modules. These can be installed with pip:

**(Recommended)**
```
pip install git+https://github.com/emanuelefalbo/HessFit.git
```
**(Dev Option)** or  by cloning the directory locally, user should give:
```
git clone https://github.com/emanuelefalbo/HessFit.git
cd HessFit
pip install -e .
```

Alternatively, it is sufficient to add their directory to the **PYTHONPATH** into you .bashrc (bash) file to call the main executable (**hessfit.py**):
```
export PYTHONPATH="${PYTHONPATH}:/path/to/hessfit/"
export PATH=$PATH:/path/to/hessfit/
```
with path/to/hessfit being the full path to where it is located. 
The last line add the programs to your bash path (see setenv for .tcsh) to make it visibile anywhere.

# Usage 

**1.1 Harmonic Force Field**

The program is thought of as dual usage, i.e., it can be executed by launching the **hessfit.py** script that performs all 
operations provided by the *build_4_hessfit.py* first, and *hessfit_harmonic.py* secondly, or these two scripts can be run independently.
It works with a json file as described below:

```
hessfit.py  step1.json -path $g09root --at gaff --test True
```

The step1.json is composed as :

```json
{
    "files": {
        "log_qm_file":  "file.log",
        "fchk_qm_file": "file.fchk",
        "fchk_mm_file": "GauHarm.fchk",
        "fchk_nb_file": "GauNonBon.fchk"
         },
    "mode": "mean",
    "charge":0,
    "multiplicity":1,
    "opt": "ric"
}
```

where --at option indicate to use the AMBER or GAFF atom type options. According to this variable, VDW parameters will be taken from amber.prm file in Gaussian path or extracted from tabulated data for GAFF. 

```
 "log_qm_file":  Gaussian output log file of QM (opt+Freq) calculation
 "fchk_qm_file": Gaussian format check-point file of QM (opt+freq) calculation
 "fchk_mm_file": Gaussian format check-point file of artificial MM (freq) calculation
 "fchk_nb_file": Gaussian format check-point file of artificial MM (freq) calculation
 "mode": "mean" string averages all force contsants over same types, while "all" leaves it unchanged
 "charge": the molecular charge of compound
 "multiplicity": the molecular multiplicity accoding to spin state
 "opt": "ric" string perform the Hessian diagonalization in redundant internal coordinates, whereas "sem" string performs the Seminario method.

```
An example of input file for gaussian is the following:
```
%mem=1GB
%nprocshared=1
%chk=but_qm.chk
#P B3LYP/def2TZVPP Geom=Connectivity opt=(calcall,tight,maxstep=7,maxcycles=100) Freq=intmodes

title

0 1
...

```

GauHarm.gjf and GauHarm.gjf are the necessary input files for *hessfit_harmonic.py* and need to be specified in the step1.json, although these are created during the run.
A rapid test of HessFit FF can be obtained using MM Gaussian routines with --test TRUE option. 
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


# HessFit Conformer Selection Workflow

## Overview

HessFit can evaluate several candidate molecular conformers and identify the structure that provides the best agreement between the QM reference calculation and the HessFit force-field representation.

The workflow consists of:

1. Generate candidate conformers.
2. Check their RICs and geometries.
3. Screen them with xTB.
4. Select representative conformers using energy and RMSD.
5. Perform high-level Gaussian `Opt+Freq` calculations.
6. Run HessFit independently for each selected conformer.
7. Compare QM and HessFit geometries, RIC Hessians, and vibrational frequencies.
8. Calculate the merit function.
9. Rank the conformers and select the optimal final geometry.

The QM RIC coordinates and QM RIC Hessian are used directly. The Cartesian QM Hessian is not converted into a RIC Hessian for this workflow.

---

# 1. Generate initial conformers

Starting from the initial molecular structure, construct the canonical redundant internal coordinate (RIC) representation containing bonds, angles, and dihedrals.

Example:

```text
Number of atoms       : 14
Number of RICs        : 64
Number of bonds       : 13
Number of angles      : 24
Number of dihedrals   : 27
```

Generate candidate conformers:

```bash
python3 generate_conformers.py
```

Structures are written to:

```text
conformers/
```

For example:

```text
conformers/conf_000.xyz
conformers/conf_001.xyz
...
```

---

# 2. Check the generated conformers

Run:

```bash
./check_conformers.py
```

Verify that:

- all structures contain the expected number of atoms;
- all structures use the same RIC definition;
- bond lengths remain reasonable;
- the intended dihedral degrees of freedom change.

Visually inspect the XYZ structures when appropriate.

---

# 3. Screen the conformers with xTB

Perform an inexpensive QM screening before expensive Gaussian calculations.

For each conformer, perform an xTB optimization. This stage is intended to remove high-energy and structurally redundant structures, not to provide the final HessFit reference geometry.

A test calculation can be run with:

```bash
screen_xtb_conformers.py
```

Successful calculations produce optimized structures such as:

```text
xtb_screening/conf_000/xtbopt.xyz
```

---

# 4. Select representative conformers using energy and RMSD

Run:

```bash
./select_xtb_conformers.py \
    --energy-cutoff 5.0 \
    --rmsd-cutoff 0.5 \
    --max-conformers 5
```

The options mean:

- `--energy-cutoff 5.0`: retain structures within 5 kcal/mol of the lowest-energy conformer;
- `--rmsd-cutoff 0.5`: treat structures below 0.5 Å RMSD as redundant;
- `--max-conformers 5`: retain at most five representatives.

Example:

```text
Successful xTB calculations : 20
Within energy cutoff         : 20
Unique RMSD clusters         : 2
Final conformers selected    : 2
```

The selected structures should represent both energetic and structural diversity.

---

# 5. Prepare Gaussian Opt+Freq calculations

For each selected conformer, create a Gaussian input and perform geometry optimization followed by a frequency calculation.

Use the same method, basis set, charge, multiplicity, and solvent model for all conformers.

For example:

```text
%mem=96GB
%nprocshared=48
%chk=01_conf_000.chk
#p m052x/6-31g+(d,p) opt freq

01_conf_000

0 1
```

The resulting files should be named consistently:

```text
01_conf_000.log
01_conf_000.fchk

02_conf_010.log
02_conf_010.fchk
```

Each conformer is treated independently.

---

# 6. Preserve the QM RIC representation

For every QM calculation, HessFit uses the corresponding Gaussian information to obtain:

- optimized geometry;
- canonical RIC coordinates;
- QM RIC Hessian;
- vibrational frequencies.

The QM RIC Hessian must be used directly.

Do **not** replace it with a Cartesian Hessian that is subsequently transformed into RIC coordinates.

All conformers must use the same RIC definition. For the development test:

```text
QM RIC Hessian shape: (64, 64)
```

---

# 7. Run HessFit for every selected conformer

Run the normal HessFit procedure independently for each high-level QM conformer.

The output naming convention should be:

```text
01_conf_000_hessfit.log
01_conf_000_hessfit.fchk

02_conf_010_hessfit.log
02_conf_010_hessfit.fchk
```

This creates a one-to-one correspondence:

```text
01_conf_000.log
       |
       +-- 01_conf_000_hessfit.log
       |
       +-- 01_conf_000_hessfit.fchk
```

The HessFit FCHK contains the fully force-field-derived RIC Hessian used in the comparison.

---

# 8. Compare QM and HessFit geometries

The evaluator compares the optimized QM geometry with the final HessFit geometry.

Example:

```text
Geometry comparison:
  QM atoms       : 14
  HessFit atoms  : 14
  QM coordinates : (14, 3)
  HF coordinates : (14, 3)
```

The geometry error is represented by RMSD.

Example:

```text
01_conf_000    RMSD = 0.04126 Å
02_conf_010    RMSD = 0.09749 Å
```

The atom ordering must remain consistent between the QM and HessFit structures.

---

# 9. Compare vibrational frequencies

For a nonlinear molecule with `N` atoms:

\[
N_{vib}=3N-6
\]

For 14 atoms:

\[
N_{vib}=36
\]

Gaussian logs may contain more than one frequency calculation. The evaluator therefore:

1. reads all `Frequencies --` entries;
2. determines the expected number of modes;
3. groups the values into complete sets;
4. selects the final complete set.

Example QM output:

```text
Expected frequencies per set : 36
Total frequencies found      : 72
Complete frequency sets      : 2
Remaining frequencies        : 0
Using frequency set          : 2
Frequencies used             : 36
```

The HessFit output contains one complete set:

```text
Expected frequencies per set : 36
Total frequencies found      : 36
Complete frequency sets      : 1
Remaining frequencies        : 0
Using frequency set          : 1
Frequencies used             : 36
```

Imaginary frequencies are also monitored.

Example:

```text
01_conf_000 : modes=36  imag(QM)=0  imag(HessFit)=0
02_conf_010 : modes=36  imag(QM)=0  imag(HessFit)=0
```

---

# 10. Compare the RIC Hessians

For each conformer, extract:

- QM RIC Hessian from the QM FCHK;
- HessFit RIC Hessian from the HessFit FCHK.

Example:

```text
Reading QM RIC Hessian...
  QM RIC Hessian shape: (64, 64)
Reading HessFit RIC Hessian...
  MM RIC Hessian shape: (64, 64)
```

The comparison is performed directly in the common RIC space:

\[
\mathbf H^{QM}_{RIC}
\quad	ext{vs.}\quad
\mathbf H^{FF}_{RIC}
\]

---

# 11. Calculate the merit-function components

## 11.1 Hessian error

The normalized RIC Hessian error is:

\[
E_H =
rac{
\left\|
\mathbf H^{FF}_{RIC}-\mathbf H^{QM}_{RIC}

ight\|_F
}{
\left\|
\mathbf H^{QM}_{RIC}

ight\|_F
}
\]

where the Frobenius norm is used.

## 11.2 Geometry error

The geometry error is the QM/HessFit RMSD:

\[
E_G = RMSD(\mathbf R^{QM},\mathbf R^{FF})
\]

The conformer-ranking procedure normalizes the individual conformer errors before constructing the final merit score.

## 11.3 Frequency error

A normalized RMS frequency error can be defined as:

\[
E_\nu =
rac{
\sqrt{
rac{1}{N_{vib}}
\sum_i(
u_i^{FF}-
u_i^{QM})^2
}
}{
\sqrt{
rac{1}{N_{vib}}
\sum_i(
u_i^{QM})^2
}
}
\]

---

# 12. Combine the three components

The total merit function is:

\[
oxed{
M=w_HE_H+w_GE_G+w_\nu E_\nu
}
\]

with:

\[
w_H+w_G+w_\nu=1.
\]

The example uses:

```text
Hessian weight    : 0.500
Geometry weight   : 0.250
Frequency weight  : 0.250
```

so:

\[
M=0.50E_H+0.25E_G+0.25E_\nu
\]

The lowest-merit conformer is selected.

---

# 13. Run the complete merit evaluation

Once the QM and HessFit files have been generated:

```bash
./evaluate_hessfit_confromers.py \
    --conformers 01_conf_000 02_conf_010 \
    --qm-dir ../examples/conformers/qm_confs/ \
    --hessfit-dir ../examples/conformers/confs_runs/ \
    --weights 0.50 0.25 0.25
```

The evaluator reports:

- geometry comparison;
- vibrational-mode counts;
- frequency diagnostics;
- QM RIC Hessian dimensions;
- HessFit RIC Hessian dimensions;
- Hessian error;
- geometry RMSD;
- frequency error;
- normalized merit components;
- final ranking.

Results are also written to:

```text
conformer_merit.csv
```

---

# 14. Example result

For the current development test:

```text
Conformer                 ΔE      Hessian       RMSD    Frequency      Merit
----------------------------------------------------------------------------------------------------
01_conf_000            0.000     0.147109    0.04126     0.072943   0.000000
02_conf_010            0.603     0.153490    0.09749     0.088284   1.000000
```

Normalized components:

```text
Conformer               Hessian     Geometry    Frequency        Merit
------------------------------------------------------------------------
01_conf_000             0.00000      0.00000      0.00000      0.00000
02_conf_010             1.00000      1.00000      1.00000      1.00000
```

The selected geometry is:

```text
OPTIMAL CONFORMER: 01_conf_000
Merit score: 0.000000
```

---

# 15. Recommended directory organization

A convenient layout is:

```text
examples/conformers/
│
├── initial/
│   └── initial.xyz
│
├── conformers/
│   ├── conf_000.xyz
│   ├── conf_001.xyz
│   └── ...
│
├── xtb_screening/
│   ├── conf_000/
│   │   └── xtbopt.xyz
│   └── ...
│
├── qm_confs/
│   ├── 01_conf_000.log
│   ├── 01_conf_000.fchk
│   ├── 02_conf_010.log
│   └── 02_conf_010.fchk
│
└── confs_runs/
    ├── 01_conf_000_hessfit.log
    ├── 01_conf_000_hessfit.fchk
    ├── 02_conf_010_hessfit.log
    └── 02_conf_010_hessfit.fchk
```

---

# 16. Choosing the number of conformers

The purpose of the xTB stage is to avoid performing expensive high-level `Opt+Freq` calculations for every generated structure.

A practical workflow is:

```text
Initial conformers
       |
       v
20–100 structures
       |
       v
xTB screening
       |
       v
Energy filtering
       |
       v
RMSD clustering
       |
       v
3–5 representative conformers
       |
       v
High-level QM Opt+Freq
       |
       v
HessFit
       |
       v
Merit function
       |
       v
Optimal final geometry
```

For very small ligands, more conformers can be screened because the xTB stage is relatively inexpensive.

---

# 17. User-selectable conformers

The user can decide which conformers proceed to the expensive QM/HessFit stage.

For example:

```bash
./evaluate_hessfit_confromers.py \
    --conformers 01_conf_000 02_conf_010 03_conf_014 \
    --qm-dir ../examples/conformers/qm_confs/ \
    --hessfit-dir ../examples/conformers/confs_runs/ \
    --weights 0.50 0.25 0.25
```

This allows selection based on the xTB results, RMSD clusters, or chemically interesting structures.

---

# 18. Interpretation of the final result

The optimal conformer is **not necessarily the lowest-energy conformer**.

The objective is to identify the geometry for which the HessFit force-field representation gives the best overall agreement with the QM reference.

The three components address different properties:

| Component | Question |
|---|---|
| RIC Hessian | Does the force field reproduce the local QM curvature? |
| Geometry RMSD | Does HessFit preserve the QM equilibrium structure? |
| Vibrational frequencies | Does the force field reproduce the vibrational response? |

Therefore, the final geometry is selected according to the quality of the **force-field-derived model**, rather than molecular energy alone.

---

# 19. Complete workflow

For a new ligand:

```text
1. Prepare initial ligand
        |
2. Build canonical RIC
        |
3. Identify rotatable bonds
        |
4. Generate multiple conformers
        |
5. Validate RICs and geometries
        |
6. xTB screen/optimize
        |
7. Energy filtering
        |
8. RMSD clustering
        |
9. Select 3–5 representatives
        |
10. High-level Gaussian Opt+Freq
        |
11. Extract QM geometry
        |
12. Extract QM RIC Hessian
        |
13. Run HessFit for each conformer
        |
14. Extract HessFit geometry
        |
15. Extract HessFit RIC Hessian
        |
16. Extract QM/HessFit frequencies
        |
17. Calculate Hessian error
        |
18. Calculate geometry error
        |
19. Calculate frequency error
        |
20. Combine weighted merit function
        |
21. Rank conformers
        |
22. Select minimum-merit geometry
        |
23. Use selected geometry for final HessFit parameterization
```

## Important consistency requirements

- Keep the same atom ordering for all conformers.
- Keep the same molecular connectivity.
- Keep the same canonical RIC definition.
- Use the same QM method/basis/charge/multiplicity/solvent model.
- Use the same HessFit protocol.
- Compare QM and HessFit Hessians directly in the common RIC space.
- Use the corresponding QM and HessFit files for each conformer.
- Verify that the expected number of vibrational modes is obtained.
- Monitor imaginary frequencies before accepting a conformer.

