#!/usr/bin/env python3

import parser_gau as pgau
import set_module as smod
import sys
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
    parser.add_argument('-m', '--mode', choices=['seminario', 'ric'],
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
    # Open texts contents : XYZ, Grad, Hess, Topology & etc
    N_atoms, qm_XYZ = pgau.read_XYZ(text_qm_fchk)                 # Reading in QM CC XYZ
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
    
    # Solve Linear System for Bond and Angles only 
    # H_MM * K = H_QM ; ignoring Torsion 
    coeffs = np.linalg.solve(MM_diag, diag_QM)
    k_bonds = coeffs[0 : No_bonds]
    k_angles = coeffs[No_bonds : No_bonds + No_angles ]
    # print(coeffs)
    # print("")
    # [print(i,j) for i,j in enumerate(k_bonds)]
    # print("")
    # [print(i,j) for i,j in enumerate(k_angles)]

    mdout = 'mean'
    bond_type_list, bond_arr, k_bond_arr = smod.set_bonds(qm_XYZ, ele_list, type_list, \
                   bond_list, k_bonds, mdout)
    # print(bond_type_list)

    for k in range(len(bond_type_list)):
            msg = (
                  f'HrmStr1 {bond_type_list[k]}  {k_bond_arr[k]:.3f} ' 
                  f' {bond_arr[k]:.3f} '
                  )
            print(msg)



if __name__ == "__main__":
    main()