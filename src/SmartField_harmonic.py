#!/usr/bin/env python3

import parser_gau as pgau
import force_constant_mod as fc
import average_across_types as aat
# import printout_mod as pout
import argparse
import os
import json
import print_top as top
import numpy as np
# import scipy.optimize as optimize


def print_init():
     print("""
 ======================================================
   Program:      SmartField
   Creator:      Emanuele Falbo, Napoli
   Date:         October 2023
   Language:     Python 3.v later
   Description:  The program returns force constants for
                 bonded values and non-bonded parameters
   Mail:         falbo.emanuele@gmail.com 
 =======================================================
     
    """)

def commandline_parser():
    parser = argparse.ArgumentParser(prog='smart_bonds_angles.py', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('optfile', nargs='?', help='option file in json')
    parser.add_argument('-m', '--mode', choices=['modsem', 'ric'],
                        default='ric', help='method to compute harmonic factors')
    return parser


def read_optfile(fname):
    """
    Read the option file. expected in json style
    """
    if not os.path.exists(str(fname)):
        raise FileNotFoundError('Missing JSON file')
    with open(fname, 'r')  as fopen:
        data = json.load(fopen)
    for i in ["log_qm_file", "fchk_qm_file", "fchk_mm_file", "fchk_nb_file"]:
        if not os.path.exists(data['files'][i]):
            raise FileNotFoundError(f'Missing {i} file')
    return data


def get_DiagMatrix(AM):
    for i in range(len(AM)):
        for j in range(len(AM)):
            if i != j :
                 AM[i,j] = 0.0
    return AM


def main():
    print_init()
    parser = commandline_parser()
    opts = parser.parse_args()
    json_opts = read_optfile(opts.optfile)
    f_qm_log = json_opts['files']['log_qm_file']
    f_qm_fchk = json_opts['files']['fchk_qm_file']
    f_mm_fchk = json_opts['files']['fchk_mm_file']
    f_nb_fchk = json_opts['files']['fchk_nb_file']
    # Store all fiels in texts
    text_qm_log = pgau.store_any_file(f_qm_log)
    text_qm_fchk = pgau.store_any_file(f_qm_fchk)
    text_mm_fchk = pgau.store_any_file(f_mm_fchk)
    text_nb_fchk = pgau.store_any_file(f_nb_fchk)
    # Opening texts contents : XYZ, Grad, Hess, Topology & etc
    N_atoms, qm_XYZ = pgau.read_XYZ(text_qm_fchk)                 
    ele_list, type_list = pgau.read_NamesTypes(text_qm_log, N_atoms)
    ric_list, force_1D = pgau.read_RicDim_Grad(text_qm_fchk)
    No_ric = ric_list[0]
    No_bonds = ric_list[1]
    No_angles = ric_list[2]
    No_dihes = ric_list[3]

    # Reading in Topology in RIC from log file
    chg = pgau.read_CM5(text_qm_log, N_atoms)
    bond_list, angle_list, tors_list = pgau.read_Top(text_qm_log, ric_list)

    if json_opts['opt'] == 'modsem':
        # Reading XYZ Hessians
        hessXYZ_qm = pgau.read_HessXYZ(text_qm_fchk, N_atoms)
        hessXYZ_nb = pgau.read_HessXYZ(text_nb_fchk, N_atoms)
        hess_eff = hessXYZ_qm - hessXYZ_nb 
        k_bonds = np.empty(len(bond_list))
        k_angles = np.empty(len(angle_list))
        k_tors = np.empty(len(tors_list))
        diag_QM = np.diagonal(hessXYZ_qm)  
        # coeffs = np.linalg.lstsq(hessXYZ_mm, diag_QM, rcond=-1)[0]
    else:
        # Reading RIC Hessains 
        hessRIC_qm = pgau.read_HessRIC(text_qm_fchk, ric_list)
        hessRIC_mm = pgau.read_HessRIC(text_mm_fchk, ric_list)
        hessRIC_nb = pgau.read_HessRIC(text_nb_fchk, ric_list)
        hess_eff = hessRIC_qm - hessRIC_nb
        diag_QM = np.diagonal(hess_eff)                 # Take diagonal items of H_QM
        MM_diag = get_DiagMatrix(hessRIC_mm)            # Make sure H_MM is diagonal
        # print(f'{np.array2string(MM_diag, precision=2)}')
        coeffs = np.linalg.solve(MM_diag, diag_QM)      # Solve Linear System for Bond and Angles only 
        # fit = np.linalg.solve(hessRIC_mm, diag_QM) #*  627.509391   # Solve Linear System for Bond and Angles only 
                                                        # H_MM * K = H_QM ; ignoring Torsion 
           
        k_bonds = coeffs[0 : No_bonds]  #  * ((627.509391)/(0.529117*0.529117))                       
        k_angles = coeffs[No_bonds : No_bonds + No_angles ] #* (627.509391)
        k_tors = coeffs[No_ric - No_dihes : No_ric] #* 627.509391  # Torsional Gradient; kcal/mol rad
    
    print(k_tors)
    
    mdin = json_opts['opt']
    mode = json_opts['mode']
    
    bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, hess_eff, type_list, bond_list, k_bonds, mdin, mode)
    angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, hess_eff, type_list, angle_list, k_angles, mdin, mode)
    tors_type_list, v1, v2, v3, tors_arr, phase, periodic_list = fc.set_torsion(qm_XYZ, type_list, tors_list, k_tors, force_1D, mode)
    
    
    # Take out mirrored atom types of bonds & angles
    bonds_unique, k_bonds_unique = aat.make_list_unique(bond_type_list, k_bond_arr)
    angles_unique, k_angles_unique = aat.make_list_unique(angle_type_list, k_angle_arr)
    # Same for torsion...
    tors_unique, val_unique = aat.make_list_unique(tors_type_list, v1)
    # print(val_unique)
    # [print(i) for i in zip(tors_type_list)]
    # print("")
    # [print(j) for j in zip(tors_unique)]
    
    
    
    # Print all into Gaussian Input
    top.print_GauInp(ele_list, type_list, qm_XYZ, \
                 bonds_unique, k_bonds_unique, bond_arr, \
                 angles_unique, k_angles_unique, angle_arr, \
                 tors_unique, v1, v2, v3, phase, periodic_list, chg)

    top.print_AmbFrcmod(type_list, \
                 bonds_unique, k_bonds_unique, bond_arr, \
                 angles_unique, k_angles_unique, angle_arr, \
                 tors_type_list, v1, v2, v3, phase, periodic_list)


if __name__ == "__main__":
    main()
