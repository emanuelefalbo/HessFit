#!/usr/bin/env python3

import parser_gau as pgau
import force_constant_mod as fc
import average_across_types as aat
# import printout_mod as pout
import os
import readin_opts as rdin
import print_top as top
import log2topol
import numpy as np
# import scipy.optimize as optimize
from scipy.sparse import rand
from scipy.optimize import lsq_linear


def get_DiagMatrix(AM):
    for i in range(len(AM)):
        for j in range(len(AM)):
            if i != j :
                 AM[i,j] = 0.0
    return AM


def solve_lsq(A, b):
    rng = np.random.default_rng()
    b = rng.standard_normal(A.shape[0])
    lb = rng.standard_normal(b.shape[0])
    ub = lb + 1
    res = lsq_linear(A, b, bounds=(lb, ub), lsmr_tol='auto', verbose=1)
    print(res)

def fit_hessian(qm_hessian, md_hessian):
    # Extract 3x3 submatrices
    N_atoms = np.int(qm_hessian.shape[0]/3)
    subs_mm, subs_qm = [], []
    for i in range(N_atoms-1):
        for j in range(i+1, N_atoms):
            sub_mm = md_hessian[3*i:3*(i+1), 3*j:3*(j+1)]
            sub_qm = qm_hessian[3*i:3*(i+1), 3*j:3*(j+1)]
            # print(f'{np.array2string(submatrix, precision=2, suppress_small=True)}')
            subs_mm.append(sub_mm)
            subs_qm.append(sub_qm)
            
    # Display extracted submatrices
    for idx, (qm, mm) in enumerate(zip(subs_qm, subs_mm), start=0):
        print(f"Submatrix {idx}:")
        # print(submatrix)
        print(np.array2string(mm, precision=2, suppress_small=True))
        # print(np.linalg.lstsq(mm.flatten(), qm.flatten(), rcond=-1)[0])
        fit = lsq_linear(qm.reshape(-1,1), mm.flatten()).x
        print(fit)


def main():
    parser = rdin.commandline_parser1()
    opts = parser.parse_args()
    json_opts = rdin.read_optfile(opts.optfile)
    f_qm_log = json_opts['files']['log_qm_file']
    f_qm_fchk = json_opts['files']['fchk_qm_file']
    f_atype = json_opts['files']['atype_file']
    f_mm_fchk = json_opts['files']['fchk_mm_file']
    f_nb_fchk = json_opts['files']['fchk_nb_file']
    formal_chg = json_opts['charge']
    multi = json_opts['multiplicity']
    
    # Store all fiels in texts
    text_qm_log = pgau.store_any_file(f_qm_log)
    text_qm_fchk = pgau.store_any_file(f_qm_fchk)
    text_atype = pgau.store_any_file(f_atype)
    text_mm_fchk = pgau.store_any_file(f_mm_fchk)
    text_nb_fchk = pgau.store_any_file(f_nb_fchk)
    
    # Opening texts contents : XYZ, Grad, Hess, Topology & etc
    N_atoms, qm_XYZ = pgau.read_XYZ(text_qm_fchk)                 
    ele_list, atype_list = pgau.read_NamesTypes(text_atype)
    ric_list, force_1D = pgau.read_RicDim_Grad(text_qm_fchk)
    No_ric = ric_list[0]
    No_bonds = ric_list[1]
    No_angles = ric_list[2]
    No_dihes = ric_list[3]

    # Reading in Topology in RIC from log file
    charge = pgau.read_CM5(text_qm_log, N_atoms)
    bond_list, angle_list, tors_list = pgau.read_Top(text_qm_log, ric_list)

    if json_opts['opt'] == 'modsem':
        # Reading XYZ Hessians
        hessXYZ_qm = pgau.read_HessXYZ(text_qm_fchk, N_atoms)
        hessXYZ_nb = pgau.read_HessXYZ(text_nb_fchk, N_atoms)
        # hessXYZ_mm = pgau.read_HessXYZ(text_mm_fchk, N_atoms)
        hess_eff = hessXYZ_qm - hessXYZ_nb 
        k_bonds = np.empty(len(bond_list))
        k_angles = np.empty(len(angle_list))
        k_tors = np.empty(len(tors_list))
        diag_QM = np.diagonal(hessXYZ_qm)  
    else:
        # Reading RIC Hessians 
        hessRIC_qm = pgau.read_HessRIC(text_qm_fchk, ric_list)
        hessRIC_mm = pgau.read_HessRIC(text_mm_fchk, ric_list)
        hessRIC_nb = pgau.read_HessRIC(text_nb_fchk, ric_list)
        hess_eff = hessRIC_qm - hessRIC_nb
        diag_QM = np.diagonal(hess_eff)                       # Take diagonal items of H_QM
        MM_diag = get_DiagMatrix(hessRIC_mm)                  # Make sure H_MM is diagonal
        coeffs = np.linalg.solve(MM_diag, diag_QM)            # Solve Linear System for Bond and Angles only H_MM*K = H_QM ; ignoring Torsion 
        k_bonds = coeffs[0 : No_bonds]                        #  * ((627.509391)/(0.529117*0.529117))                       
        k_angles = coeffs[No_bonds : No_bonds + No_angles ]   #* (627.509391)
        k_tors = coeffs[No_ric - No_dihes : No_ric]           #* 627.509391  # Torsional Gradient; kcal/mol rad
    
    # If opt == sem
    mdin = json_opts['opt']
    mode = json_opts['mode']
    bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, hess_eff, atype_list, bond_list, k_bonds, mdin, mode)
    angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, hess_eff, atype_list, angle_list, k_angles, mdin, mode)
    tors_type_list, v1, v2, v3, tors_arr, phase, periodic_list = fc.set_torsion(qm_XYZ, atype_list, tors_list, k_tors, force_1D, mode)
    
    # Take out mirrored atom types of bonds & angles
    bonds_unique, k_bonds_unique = aat.make_list_unique(bond_type_list, k_bond_arr)
    angles_unique, k_angles_unique = aat.make_list_unique(angle_type_list, k_angle_arr)
    # Same for torsion...
    tors_unique, val_unique = aat.make_list_unique(tors_type_list, v1)

    # Print all into Gaussian Input
    top.print_GauInp(ele_list, atype_list, qm_XYZ, \
                 bonds_unique, k_bonds_unique, bond_arr, \
                 angles_unique, k_angles_unique, angle_arr, \
                 tors_unique, v1, v2, v3, phase, periodic_list, charge, \
                 formal_chg, multi)

    top.print_AmbFrcmod(atype_list, \
                 bonds_unique, k_bonds_unique, bond_arr, \
                 angles_unique, k_angles_unique, angle_arr, \
                 tors_type_list, v1, v2, v3, phase, periodic_list)
    
   # Make dihedral directory for subsquent torsion fitting
    fname = 'SmartField4gau.gjf'
    top.build_dihe_folder(fname, atype_list, charge)
    
    # Build topolo.txt file 
    bond_list, angle_list, tors_list = log2topol.read_Top(text_qm_log)
    os.chdir('dihedrals')
    log2topol.print_topol(bond_list, angle_list, tors_list)
    
    

if __name__ == "__main__":
    main()
