#!/usr/bin/env python3

import argparse
import parser_gau as pgau
import force_constant_mod as fc
import numpy as np
import average_across_types as aat
import printout_mod as pout
# import pathlib
import os

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def commandline_parser():
    parser = argparse.ArgumentParser(prog='build_4Smart.py', formatter_class=argparse.RawTextHelpFormatter)
    requiredNamed = parser.add_argument_group('mandatory arguments')
    requiredNamed.add_argument('-f1','--log_file', help='Gaussian QM log file ')
    requiredNamed.add_argument('-f2','--fchk_file', help='Gaussain QM fchk file')
    parser.add_argument('-m', '--mode', choices=['all', 'mean'],
                        default='mean', help='averaging across same types; default = mean')
    txt = """path/to/amber.prm in Gaussain root directory
default = current directory """
    parser.add_argument('-path', nargs='?', help=txt, default=os.getcwd(), type=dir_path)
    
    return parser

def reduce_bond_list(bond_list):
    splitted = list()
    [ splitted.append(i.split()) for i in bond_list]
    data = { (tuple(sorted(item))) for item in splitted}
    new = list(map(list, data))
    reduced = list(' '.join(item) for item in new)
    return reduced

def reduce_angle_list(angle_list):
    splitted = [ i.split() for i in angle_list]
    res = splitted.copy()
    for i in range(len(splitted[:])-1):
        for j in range(i+1,len(splitted[:])):
            if (splitted[i] == splitted[j][::-1]) and (i !=j) :
                res.pop(i)
                break
            
    reduced = list(' '.join(item) for item in res)
    return reduced

def main():
    parser = commandline_parser()
    opts = parser.parse_args()
    f_qm_log = opts.log_file
    f_qm_fchk = opts.fchk_file
    text_qm_log = pgau.store_any_file(f_qm_log)
    text_qm_fchk = pgau.store_any_file(f_qm_fchk)
    
    N_atoms, qm_XYZ = pgau.read_XYZ(text_qm_fchk)                 
    ele_list, type_list = pgau.read_NamesTypes(text_qm_log, N_atoms)
    ric_list, force_1D = pgau.read_RicDim_Grad(text_qm_fchk)
    No_ric = ric_list[0]
    No_bonds = ric_list[1]
    No_angles = ric_list[2]
    No_dihes = ric_list[3]

    # # Reading in Topology in RIC from log file
    chg = pgau.read_CM5(text_qm_log, N_atoms)
    bond_list, angle_list, tors_list = pgau.read_Top(text_qm_log, ric_list)
    k_bonds = np.ones(No_bonds)                          
    k_angles = np.ones(No_angles)


    if opts.mode == 'mean':
        bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, None, type_list, \
                      bond_list, k_bonds, 'ric', 'mean')
        angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, None, type_list, \
                      angle_list, k_angles, 'ric', 'mean')
    elif opts.mode == 'all':
        bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, None, type_list, \
                      bond_list, k_bonds, 'ric', 'all')
        angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, None, type_list, \
                      angle_list, k_angles, 'ric', 'all') 


    # Take out Mirror atom types of bonds & angles
    bond_reduced = reduce_bond_list(bond_type_list)
    angle_reduced, k_angles_unique = aat.reduce_angle_list(angle_type_list, k_angle_arr)

    # Print all into Gaussian Input
    pout.print_GauHarm(ele_list, type_list, qm_XYZ, \
                 bond_reduced, k_bond_arr, bond_arr, \
                 angle_reduced, k_angles_unique, angle_arr, \
                 chg)

    VDW_list = pgau.read_AmberParm(opts.path, type_list)
    
    pout.print_GauNonBon(ele_list, type_list, qm_XYZ, \
                 bond_reduced, bond_arr, \
                 angle_reduced, k_angles_unique, \
                 chg, VDW_list )




if __name__ == '__main__':
    main()