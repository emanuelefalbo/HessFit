#!/usr/bin/env python3

import argparse
from turtle import clear
from unittest import result
import parser_gau as pgau
import force_constant_mod as fc
import numpy as np
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


def print_GauHarm(*arg):
    ele_ls = arg[0]
    tp_ls= arg[1]
    coord = arg[2]
    bond_tp_ls = arg[3]
    k_bond_ls = arg[4]
    bond_ln_ls = arg[5]
    angle_tp_ls = arg[6]
    k_angle_ls = arg[7]
    angle_ln_ls = arg[8]
    chg = arg[9]

    header_gjf ="""%mem=1GB
%nprocshared=1
%chk=GauHarm.chk
#p Amber=(SoftOnly,Print) nosymm  Freq=intmodes
 
Title

0 1
"""

    master_func = """
! Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
"""
    fname = 'GauHarm.gjf'
    with open(fname, 'w') as fout:
        fout.write(header_gjf)
        for m, p, l in zip(ele_ls, tp_ls, coord):
                    #  s1 = '  '.join(str(x) for x in p)
                     s2 = '  '.join((f'{x:.6f}') for x in l)
                     fout.write(f'{m}-{p}-0.00  {s2}\n')
        # fout.write(f'\n')
        fout.write(master_func)
        fout.write(f'! SMARTFIELD FF\n')
        fout.write(f'!Bonds\n')
        for k in range(len(bond_tp_ls)):
            msg = (
                  f'HrmStr1 {bond_tp_ls[k]}  1.0 ' 
                  f' {bond_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'!Angles\n')
        for k in range(len(angle_tp_ls)):
            msg = (
                  f'HrmBnd1 {angle_tp_ls[k]}  1.0 ' 
                  f' {angle_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'\n')


def print_GauNonBon(*arg):
    ele_ls = arg[0]
    tp_ls= arg[1]
    coord = arg[2]
    bond_tp_ls = arg[3]
    bond_ln_ls = arg[4]
    angle_tp_ls = arg[5]
    angle_ln_ls = arg[6]
    chg = arg[7]
    vdw = arg[8]

    header_gjf ="""%mem=1GB
%nprocshared=1
%chk=GauNonBon.chk
#p Amber=(SoftOnly,Print) nosymm Freq=intmodes
 
Title

0 1
"""

    master_func = """
! Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
"""
    fname = 'GauNonBon.gjf'
    with open(fname, 'w') as fout:
        fout.write(header_gjf)
        for m, p, l, q in zip(ele_ls, tp_ls, coord, chg):
                #   s1 = '  '.join(str(x) for x in p)
                s2 = '  '.join((f'{x:.6f}') for x in l)
                fout.write(f'{m}-{p}-{q}  {s2}\n')
        # fout.write(f'\n')
        fout.write(master_func)
        fout.write(f'! SMARTFIELD FF\n')
        fout.write(f'!Bonds\n')
        for k in range(len(bond_tp_ls)):
            msg = (
                  f'HrmStr1 {bond_tp_ls[k]}  0.0 ' 
                  f' {bond_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'!Angles\n')
        for k in range(len(angle_tp_ls)):
            msg = (
                  f'HrmBnd1 {angle_tp_ls[k]}  0.0 ' 
                  f' {angle_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'!VDW\n')
        for k in vdw:
            s = ' '.join(k)
            fout.write(f'{s} \n')
        fout.write(f'\n')

def reduce_bond_list(bond_list):
    tmp = list()
    [ tmp.append(i.split()) for i in bond_list]
    data = { (tuple(sorted(item))) for item in tmp}
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
    angle_reduced = reduce_angle_list(angle_type_list)

    # Print all into Gaussian Input
    print_GauHarm(ele_list, type_list, qm_XYZ, \
                 bond_reduced, k_bond_arr, bond_arr, \
                 angle_reduced, k_angle_arr, angle_arr, \
                 chg)

    VDW_list = pgau.read_AmberParm(opts.path, type_list)
    
    print_GauNonBon(ele_list, type_list, qm_XYZ, \
                 bond_reduced, bond_arr, \
                 angle_type_list, angle_arr, \
                 chg, VDW_list )




if __name__ == '__main__':
    main()