#!/usr/bin/env python3

def print_GauInp(*args):
    ele_list, type_list, coord, bond_type_list, k_bond_list, bond_length_list, angle_type_list, k_angle_list, angle_length_list, torsion_type_list, v1_list, v2_list, v3_list, phase_list, hybrid_list, charges = args

    header_gjf = """%mem=1GB
%nprocshared=1
%chk=SmartField4gau.chk
#p Amber=(SoftFirst,Print) nosymm geom=nocrowd opt Freq

Title

0 1
"""

    master_function = """
! Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
"""
    filename = 'SmartField4gau.gjf'
    with open(filename, 'w') as file_out:
        file_out.write(header_gjf)
        for element, type, coordinates, charge in zip(ele_list, type_list, coord, charges):
            formatted_coords = '  '.join(f'{x:.6f}' for x in coordinates)
            file_out.write(f'{element}-{type}-{charge}  {formatted_coords}\n')

        file_out.write(master_function)
        file_out.write('! SMARTFIELD FF\n')
        write_bonds(file_out, bond_type_list, k_bond_list, bond_length_list)
        write_angles(file_out, angle_type_list, k_angle_list, angle_length_list)
        write_torsions(file_out, torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list)
        file_out.write('\n')


def write_bonds(file, bond_type_list, k_bond_list, bond_length_list):
    file.write('!Bonds\n')
    for bond_type, k_bond, bond_length in zip(bond_type_list, k_bond_list, bond_length_list):
        file.write(f'HrmStr1 {bond_type}  {k_bond:.3f} {bond_length:.3f}\n')


def write_angles(file, angle_type_list, k_angle_list, angle_length_list):
    file.write('!Angles\n')
    for angle_type, k_angle, angle_length in zip(angle_type_list, k_angle_list, angle_length_list):
        file.write(f'HrmBnd1 {angle_type}  {k_angle:.3f} {angle_length:.3f}\n')


def write_torsions(file, torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list):
    file.write('! Torsions\n')
    for torsion_type, phase, v1, v2, v3, hybrid in zip(torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list):
        formatted_phase = '  '.join(f'{x}' for x in phase)
        file.write(f'AmbTrs {torsion_type} {formatted_phase} {v1:.2f} {v2:.2f} {v3:.2f} 0. {float(hybrid)}\n')


def print_AmbFrcmod(*args):
    type_list, bond_type_list, k_bond_list, bond_length_list, angle_type_list, k_angle_list, angle_length_list, torsion_type_list, v1_list, v2_list, v3_list, phase_list, hybrid_list = args

    """
    Write an AMBER-like frcmod file:
    type_list: atom type list
    bond_type_list: bond type list
    bond_length_list: bond length list 
    k_bond_list: bond constant list
    angle_type_list: angle type list 
    k_angle_list: angle constant list 
    angle_length_list: degree angle list
    """

    type_list_unique = set(type_list)
    filename = 'SmartField_frcmod.txt'
    with open(filename, 'w') as file_out:
        file_out.write('MASS\n')
        for type_entry in type_list_unique:
            file_out.write(f'{type_entry}\n')

        write_bonds_amber(file_out, bond_type_list, k_bond_list, bond_length_list)
        write_angles_amber(file_out, angle_type_list, k_angle_list, angle_length_list)
        write_torsions_amber(file_out, torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list)


def write_bonds_amber(file, bond_type_list, k_bond_list, bond_length_list):
    file.write('BOND\n')
    for bond_type, k_bond, bond_length in zip(bond_type_list, k_bond_list, bond_length_list):
        file.write(f'{bond_type} {k_bond:.3f} {bond_length:.3f}\n')


def write_angles_amber(file, angle_type_list, k_angle_list, angle_length_list):
    file.write('ANGLE\n')
    for angle_type, k_angle, angle_length in zip(angle_type_list, k_angle_list, angle_length_list):
        file.write(f'{angle_type} {k_angle:.3f} {angle_length:.3f}\n')


def write_torsions_amber(file, torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list):
    file.write('DIHE\n')
    for torsion_type, phase, v1, v2, v3, hybrid in zip(torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list):
        formatted_phase = '  '.join(f'{x}' for x in phase)
        file.write(f'{torsion_type} {formatted_phase} {v1:.2f} {v2:.2f} {v3:.2f} 0. {hybrid}\n')


# def print_GauInp(*arg):
#     ele_ls = arg[0]
#     tp_ls= arg[1]
#     coord = arg[2]
#     bond_tp_ls = arg[3]
#     k_bond_ls = arg[4]
#     bond_ln_ls = arg[5]
#     angle_tp_ls = arg[6]
#     k_angle_ls = arg[7]
#     angle_ln_ls = arg[8]
#     tors_tp_ls = arg[9]
#     v1 = arg[10]
#     v2 = arg[11]
#     v3 = arg[12]
#     phase = arg[13]
#     hybrid_list = arg[14]
#     chg = arg[15]
#     """
#     Write a gaussian input file:
#     ele_ls: element list 
#     tp_ls: atom type list
#     coord: coordinates XYZ
#     bond_tp_ls:  bond type list
#     bond_ln_ls:  bond length list 
#     k_bond_ls:  bond constant list
#     angle_tp_ls: angle type list 
#     k_angle_ls: angle constant list 
#     angle_ln_ls : deg angle list

#     """
#     header_gjf ="""%mem=1GB
# %nprocshared=1
# %chk=SmartField4gau.chk
# #p Amber=(SoftFirst,Print) nosymm geom=nocrowd opt Freq
 
# Title

# 0 1
# """

#     master_func = """
# ! Master function
# NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
# """
#     fname = 'SmartField4gau.gjf'
#     with open(fname, 'w') as fout:
#         fout.write(header_gjf)
#         for m, p, l, q in zip(ele_ls, tp_ls, coord, chg):
#                     #  s1 = '  '.join(str(x) for x in p)
#                      s2 = '  '.join((f'{x:.6f}') for x in l)
#                      fout.write(f'{m}-{p}-{q}  {s2}\n')
#         # fout.write(f'\n')
#         fout.write(master_func)
#         fout.write(f'! SMARTFIELD FF\n')
#         fout.write(f'!Bonds\n')
#         for k in range(len(bond_tp_ls)):
#             msg = (
#                   f'HrmStr1 {bond_tp_ls[k]}  {k_bond_ls[k]:.3f} ' 
#                   f' {bond_ln_ls[k]:.3f}\n'
#                   )
#             fout.write(msg)
#         fout.write(f'!Angles\n')
#         for k in range(len(angle_tp_ls)):
#             msg = (
#                   f'HrmBnd1 {angle_tp_ls[k]}  {k_angle_ls[k]:.3f} ' 
#                   f' {angle_ln_ls[k]:.3f}\n'
#                   )
#             fout.write(msg)
#         fout.write(f'! Torsions\n')
#         for k in range(len(tors_tp_ls)):
#             s = '  '.join((f'{x}') for x in phase[k,:])
#             msg = (
#                   f'AmbTrs {tors_tp_ls[k]} {s} '
#                   f' {v1[k]:.2f} {v2[k]:.2f} {v3[k]:.2f} 0. ' 
#                   f' {float(hybrid_list[k])}\n'
#                   )
#             fout.write(msg)
#         fout.write(f'\n')


# def print_AmbFrcmod(*arg):
#     tp_ls= arg[0]
#     bond_tp_ls = arg[1]
#     k_bond_ls = arg[2]
#     bond_ln_ls = arg[3]
#     angle_tp_ls = arg[4]
#     k_angle_ls = arg[5]
#     angle_ln_ls = arg[6]
#     tors_tp_ls = arg[7]
#     v1 = arg[8]
#     v2 = arg[9]
#     v3 = arg[10]
#     phase=arg[11]
#     hybrid_list = arg[12]

#     """
#     Write a AMBER-like frcmod file:
#     tp_ls: atom type list
#     bond_tp_ls:  bond type list
#     bond_ln_ls:  bond length list 
#     k_bond_ls:  bond constant list
#     angle_tp_ls: angle type list 
#     k_angle_ls: angle constant list 
#     angle_ln_ls : deg angle list

#     """
#     tp_ls_unique = set(tp_ls)
#     fname = 'SmartField_frcmod.txt'
#     with open(fname, 'w') as fout:
#         fout.write('MASS\n')
#         for m in tp_ls_unique:
#             fout.write(f'{m}\n')
#         fout.write(f'BOND\n')
#         for k in range(len(bond_tp_ls)):
#             msg = (
#                   f'{bond_tp_ls[k]} {k_bond_ls[k]:.3f} ' 
#                   f' {bond_ln_ls[k]:.3f}\n'
#                   )
#             fout.write(msg)
#         fout.write(f'ANGLE\n')
#         for k in range(len(angle_tp_ls)):
#             msg = (
#                   f'{angle_tp_ls[k]}  {k_angle_ls[k]:.3f} ' 
#                   f' {angle_ln_ls[k]:.3f}\n'
#                   )
#             fout.write(msg)
#         fout.write(f'DIHE\n')
#         for k in range(len(tors_tp_ls)):
#             s = '  '.join((f'{x}') for x in phase[k,:])
#             msg = (
#                   f'{tors_tp_ls[k]} {s}'
#                   f' {v1[k]:.2f} {v2[k]:.2f} {v3[k]:.2f} 0. ' 
#                   f' {hybrid_list[k]}\n'
#                   )
#             fout.write(msg)
#         fout.write(f'\n')
    