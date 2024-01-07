#!/usr/bin/env python3

import parser_gau as pgau
import force_constant_mod as fc
import numpy as np
import average_across_types as aat
import readin_opts as rdin
import os

def print_GauHarm(*args):
    ele_ls, tp_ls, coord, bond_tp_ls, k_bond_ls, bond_ln_ls, angle_tp_ls, k_angle_ls, \
    angle_ln_ls, torsion_type_list, hybrid_list, chg = args

    header_gjf = """%mem=1GB
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
            s2 = '  '.join(f'{x:.6f}' for x in l)
            fout.write(f'{m}-{p}-0.00  {s2}\n')

        fout.write(master_func)
        fout.write('! SMARTFIELD FF\n')
        fout.write('!Bonds\n')
        for bond_tp, k_bond, bond_ln in zip(bond_tp_ls, k_bond_ls, bond_ln_ls):
            fout.write(f'HrmStr1 {bond_tp}  1.0 {bond_ln:.3f}\n')

        fout.write('!Angles\n')
        for angle_tp, k_angle, angle_ln in zip(angle_tp_ls, k_angle_ls, angle_ln_ls):
            fout.write(f'HrmBnd1 {angle_tp}  1.0 {angle_ln:.3f}\n')
            
        fout.write('!Torsions\n')
        for i, (torsion_type, hybrid) in enumerate(zip(torsion_type_list, hybrid_list)):
            if hybrid == 4.0 or hybrid == 2.0:
               formatted_phase = '0 180 0 0 '
               fout.write(f'AmbTrs {torsion_type} {formatted_phase} 0. {1.0} 0. 0. {float(hybrid)}\n')
            else:
               formatted_phase = '0 0 0 0 '
               fout.write(f'AmbTrs {torsion_type} {formatted_phase} 0. 0. {1.0} 0. {float(hybrid)}\n')
            
            
        fout.write('\n')
        
        
        
def print_GauNonBon(elements, types, coordinates, bond_types, bond_lengths, angle_types, \
                    angle_lengths, torsion_type_list, hybrid_list,  charges, vdw):
    header_gjf = """%mem=1GB
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

        for element, type_, coords, charge in zip(elements, types, coordinates, charges):
            coord_str = '  '.join(f'{x:.6f}' for x in coords)
            fout.write(f'{element}-{type_}-{charge}  {coord_str}\n')

        fout.write(master_func)
        fout.write('! SMARTFIELD FF\n')
        fout.write('!Bonds\n')
        for bond_type, bond_length in zip(bond_types, bond_lengths):
            fout.write(f'HrmStr1 {bond_type}  0.0 {bond_length:.3f}\n')

        fout.write('!Angles\n')
        for angle_type, angle_length in zip(angle_types, angle_lengths):
            fout.write(f'HrmBnd1 {angle_type}  0.0 {angle_length:.3f}\n')
            
        fout.write('!Torsions\n')
        for i, (torsion_type, hybrid) in enumerate(zip(torsion_type_list, hybrid_list)):
            if hybrid == 4.0 or hybrid == 2.0:
               formatted_phase = '0 180 0 0 '
               fout.write(f'AmbTrs {torsion_type} {formatted_phase} 0. 0. 0. 0. {float(hybrid)}\n')
            else:
               formatted_phase = '0 0 0 0 '
               fout.write(f'AmbTrs {torsion_type} {formatted_phase} 0. 0. 0. 0. {float(hybrid)}\n')
            
        fout.write('!VDW\n')
        for vdw_params in vdw:
            s = ' '.join(map(str, vdw_params))
            fout.write(f'{s}\n')

        fout.write('\n')
        

def main():
    parser = rdin.commandline_parser3()
    opts = parser.parse_args()
    json_opts = rdin.read_optfile_3(opts.optfile)
    # opts = parser.parse_args()
    f_qm_log = json_opts['files']['log_qm_file']
    f_qm_fchk = json_opts['files']['fchk_qm_file']
    f_atype = json_opts['files']['atype_file']
    text_qm_log = pgau.store_any_file(f_qm_log)
    text_qm_fchk = pgau.store_any_file(f_qm_fchk)
    text_atype = pgau.store_any_file(f_atype)
    
    N_atoms, qm_XYZ = pgau.read_XYZ(text_qm_fchk)
    ele_list, atype_list = pgau.read_NamesTypes(text_atype)             
    # ele_list, type_list = pgau.read_NamesTypes(text_qm_log, N_atoms)
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
    k_tors = np.ones(No_dihes)
    # print(No_dihes)
    
    mdin = json_opts['opt']
    mode = json_opts['mode']
    bond_type_list, bond_arr, k_bond_arr = fc.set_bonds(qm_XYZ, None, atype_list, \
                      bond_list, k_bonds, mdin, mode)
    angle_type_list, angle_arr, k_angle_arr = fc.set_angles(qm_XYZ, None, atype_list, \
                      angle_list, k_angles, mdin, mode)
    tors_type_list, v1, _, _, _, _, periodic_list = fc.set_torsion(qm_XYZ, atype_list, tors_list, k_tors, force_1D, mode)
    
    # Take out Mirror atom types of bonds & angles
    bond_reduced, _ = aat.make_list_unique(bond_type_list, k_bond_arr)
    angle_reduced, _ = aat.make_list_unique(angle_type_list, k_angle_arr)
    tors_reduced, hybrid_unique = aat.make_list_unique(tors_type_list, periodic_list)
    
    # Print all into Gaussian Input
    print_GauHarm(ele_list, atype_list, qm_XYZ, \
                 bond_reduced, k_bond_arr, bond_arr, \
                 angle_reduced, k_angle_arr, angle_arr, \
                 tors_reduced, hybrid_unique, \
                 chg)
    
    path = os.environ.get("g09root") + "/g09"
    VDW_list = pgau.read_AmberParm(path, atype_list)
    
    print_GauNonBon(ele_list, atype_list, qm_XYZ, \
                 bond_reduced, bond_arr, \
                 angle_reduced, k_angle_arr, \
                 tors_reduced, hybrid_unique, \
                 chg, VDW_list)




if __name__ == '__main__':
    main()