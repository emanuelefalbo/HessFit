#!/usr/bin/env python3

from cmath import phase
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
                 bonded and non-bonded values
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

def print_GauInp(*arg):
    ele_ls = arg[0]
    tp_ls= arg[1]
    coord = arg[2]
    bond_tp_ls = arg[3]
    k_bond_ls = arg[4]
    bond_ln_ls = arg[5]
    angle_tp_ls = arg[6]
    k_angle_ls = arg[7]
    angle_ln_ls = arg[8]
    tors_tp_ls = arg[9]
    v1 = arg[10]
    v2 = arg[11]
    v3 = arg[12]
    phase = arg[13]
    hybrid_list = arg[14]
    chg = arg[15]
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
    header_gjf ="""%mem=1GB
%nprocshared=1
%chk=SmartField4gau.chk
#p Amber=(SoftFirst,Print) nosymm geom=nocrowd opt Freq
 
Title

0 1
"""

    master_func = """
! Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
"""
    fname = 'SmartField4gau.gjf'
    with open(fname, 'w') as fout:
        fout.write(header_gjf)
        for m, p, l, q in zip(ele_ls, tp_ls, coord, chg):
                    #  s1 = '  '.join(str(x) for x in p)
                     s2 = '  '.join((f'{x:.6f}') for x in l)
                     fout.write(f'{m}-{p}-{q}  {s2}\n')
        # fout.write(f'\n')
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
        fout.write(f'! Torsions\n')
        for k in range(len(tors_tp_ls)):
            s = '  '.join((f'{x}') for x in phase[k,:])
            msg = (
                  f'AmbTrs {tors_tp_ls[k]} {s} '
                  f' {v1[k]:.2f} {v2[k]:.2f} {v3[k]:.2f} 0. ' 
                  f' {float(hybrid_list[k])}\n'
                  )
            fout.write(msg)
        fout.write(f'\n')


def print_AmbFrcmod(*arg):
    tp_ls= arg[0]
    bond_tp_ls = arg[1]
    k_bond_ls = arg[2]
    bond_ln_ls = arg[3]
    angle_tp_ls = arg[4]
    k_angle_ls = arg[5]
    angle_ln_ls = arg[6]
    tors_tp_ls = arg[7]
    v1 = arg[8]
    v2 = arg[9]
    v3 = arg[10]
    phase=arg[11]
    hybrid_list = arg[12]

    """
    Write a AMBER-like frcmod file:
    tp_ls: atom type list
    bond_tp_ls:  bond type list
    bond_ln_ls:  bond length list 
    k_bond_ls:  bond constant list
    angle_tp_ls: angle type list 
    k_angle_ls: angle constant list 
    angle_ln_ls : deg angle list

    """
    tp_ls_unique = set(tp_ls)
    fname = 'SmartField_frcmod.txt'
    with open(fname, 'w') as fout:
        fout.write('MASS\n')
        for m in tp_ls_unique:
            fout.write(f'{m}\n')
        fout.write(f'BOND\n')
        for k in range(len(bond_tp_ls)):
            msg = (
                  f'{bond_tp_ls[k]} {k_bond_ls[k]:.3f} ' 
                  f' {bond_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'ANGLE\n')
        for k in range(len(angle_tp_ls)):
            msg = (
                  f'{angle_tp_ls[k]}  {k_angle_ls[k]:.3f} ' 
                  f' {angle_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'DIHE\n')
        for k in range(len(tors_tp_ls)):
            s = '  '.join((f'{x}') for x in phase[k,:])
            msg = (
                  f'{tors_tp_ls[k]} {s}'
                  f' {v1[k]:.2f} {v2[k]:.2f} {v3[k]:.2f} 0. ' 
                  f' {hybrid_list[k]}\n'
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
        diag_tors = np.empty(len(tors_list))
    else:
        # Reading RIC Hessains 
        hessRIC_qm = pgau.read_HessRIC(text_qm_fchk, ric_list)
        hessRIC_mm = pgau.read_HessRIC(text_mm_fchk, ric_list)
        hessRIC_nb = pgau.read_HessRIC(text_nb_fchk, ric_list)
        hess_eff = hessRIC_qm - hessRIC_nb
        diag_QM = np.diagonal(hess_eff)                 # Take diagonal items of H_QM
        MM_diag = get_DiagMatrix(hessRIC_mm)            # Make sure H_MM is diagonal
        coeffs = np.linalg.solve(MM_diag, diag_QM)      # Solve Linear System for Bond and Angles only 
                                                        # H_MM * K = H_QM ; ignoring Torsion 
        k_bonds = coeffs[0 : No_bonds]                          
        k_angles = coeffs[No_bonds : No_bonds + No_angles ]
        diag_tors = diag_QM[No_ric - No_dihes : No_ric] * 627.509391  # Torsional Gradient; kcal/mol rad
        # [print(i) for i in diag_tors]

    mdin = json_opts['opt']
    print(f'Option mode =  {mdin}' )
    if json_opts['mode'] == 'mean':
        bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, hess_eff, type_list, \
                      bond_list, k_bonds, mdin, 'mean')
        angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, hess_eff, type_list, \
                      angle_list, k_angles, mdin, 'mean')
        tors_type_list, v1, v2, v3, tors_arr, phase, periodic_list = fc.set_torsion(qm_XYZ, type_list, \
                      tors_list, diag_tors, force_1D, 'mean')               
    elif json_opts['mode'] == 'all':
        bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, hess_eff, type_list, \
                      bond_list, k_bonds, mdin, 'all')
        angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, hess_eff, type_list, \
                      angle_list, k_angles, mdin, 'all') 
        tors_type_list, v1, v2, v3, tors_arr, phase, periodic_list = fc.set_torsion(qm_XYZ, type_list, \
                      tors_list, diag_tors, force_1D, 'all')


    # Print all into Gaussian Input
    print_GauInp(ele_list, type_list, qm_XYZ, \
                 bond_type_list, k_bond_arr, bond_arr, \
                 angle_type_list, k_angle_arr, angle_arr, \
                 tors_type_list, v1, v2, v3, phase, periodic_list, chg)

    print_AmbFrcmod(type_list, \
                 bond_type_list, k_bond_arr, bond_arr, \
                 angle_type_list, k_angle_arr, angle_arr, \
                 tors_type_list, v1, v2, v3, phase, periodic_list)


if __name__ == "__main__":
    main()
