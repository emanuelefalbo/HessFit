#!/usr/bin/env python3

import parser_gau as pgau
import force_constant_mod as fc
import argparse
import os
import json
import numpy as np


def print_init():
     print("""
 ======================================================
   Program:      SmartField
   Creator:      Emanuele Falbo
                 Scuola Normale Superiore, Pisa
   Date:         July 2022
   Language:     Python 3.v later
   Description:  The program returns force constants for
                 bonds and angles.
   Mail:         falbo.emanuele@gmail.com 
 =======================================================
     
    """)

def commandline_parser():
    parser = argparse.ArgumentParser(prog='smart_bonds_angles.py', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('optfile', nargs='?', help='option file in json')
    # parser.add_argument('-m', '--mode', choices=['seminario', 'ric'],
                        # default='ric', help='method to compute harmonic factors')
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

def print_GauInp(ele_ls, tp_ls, coord, \
                 bond_tp_ls, k_bond_ls, bond_ln_ls,
                 angle_tp_ls, k_angle_ls, angle_ln_ls):
    """
    Write a gaussian input file:
    ele_ls: element list 
    tp_ls: atom type list
    coord: coordinates XYZ
    bond_tp_ls:  bond type list
    bond_ln_ls:  bond length list 
    k_bond_ls:  bond constant list
    angle_tp_ls: angle type list 
    k_angle_ls: angle constant list 
    angle_ln_ls : deg angle list

    """
    header_gjf ="""%mem=4GB
%nprocshared=1
%chk=smartfield4gau.chk
#p Amber=(SoftFirst,Print) nosymm geom=nocrowd opt
 
Title

0 1
"""

    master_func = """
! Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2

"""
    fname = 'smartfield4gau.gjf'
    with open(fname, 'w') as fout:
        fout.write(header_gjf)
        for m, p, l in zip(ele_ls, tp_ls, coord):
                    #  s1 = '  '.join(str(x) for x in p)
                     s2 = '  '.join((f'{x:.6f}') for x in l)
                     fout.write(f'{m}-{p}  {s2} \n')
        fout.write(f'\n')
        fout.write(master_func)
        fout.write(f'! SMARTFIELD FF\n')
        fout.write(f'!Bonds\n')
        for k in range(len(bond_tp_ls)):
            msg = (
                  f'HrmStr1 {bond_tp_ls[k]}  {k_bond_ls[k]:.3f} ' 
                  f' {bond_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'!Angles\n')
        for k in range(len(angle_tp_ls)):
            msg = (
                  f'HrmBnd1 {angle_tp_ls[k]}  {k_angle_ls[k]:.3f} ' 
                  f' {angle_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'\n')
    

def get_DiagMatrix(AM):
    for i in range(len(AM)):
        for j in range(len(AM)):
            if i != j :
                 AM[i,j] = 0.0
    return AM


def main():
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
    No_bonds = ric_list[1]
    No_angles = ric_list[2]
    No_dihes = ric_list[3]

    # Reading in Topology in RIC from log file
    bond_list, angle_list, dihe_list = pgau.read_Top(text_qm_log, ric_list)
    hess_qm = pgau.read_Hess(text_qm_fchk, ric_list)
    hess_mm = pgau.read_Hess(text_mm_fchk, ric_list)
    hess_nb = pgau.read_Hess(text_nb_fchk, ric_list)

    diag_QM = np.diagonal(hess_qm)          # Take diagonla items of H_QM
    MM_diag = get_DiagMatrix(hess_mm)       # Make sure H_MM is diagonal
    
    coeffs = np.linalg.solve(MM_diag, diag_QM)              # Solve Linear System for Bond and Angles only 
    k_bonds = coeffs[0 : No_bonds]                          # H_MM * K = H_QM ; ignoring Torsion 
    k_angles = coeffs[No_bonds : No_bonds + No_angles ]

    if json_opts['mode'] == 'mean':
        bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, ele_list, type_list, \
                      bond_list, k_bonds, 'mean')
        angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, ele_list, type_list, \
                      angle_list, k_angles, 'mean')              
    else:
        bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, ele_list, type_list, \
                      bond_list, k_bonds, 'all')
        angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, ele_list, type_list, \
                      angle_list, k_angles, 'all')              

    # Print all into Gaussian Input
    print_GauInp(ele_list, type_list, qm_XYZ, \
                 bond_type_list, k_bond_arr, bond_arr, \
                 angle_type_list, k_angle_arr, angle_arr )
   


if __name__ == "__main__":
    main()