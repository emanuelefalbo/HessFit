#!/usr/bin/env python3

import os

def build_dihe_folder(fname, type_list, charges):
    # Function to store force field strings, and
    # atomtype-charges for the GauScan2log function
    folder = 'dihedrals'
    if not os.path.exists(folder):
       os.makedirs(folder)
       
    with open(fname , 'r') as file:
        file_content = file.read()
        split_content = file_content.split('!Master function')
        if len(split_content) > 1:
            # Storing the text after encountering the '!Master function'
            text_after_mf = split_content[1].strip()
        f1_path = os.path.join(folder, 'ff_string.txt')
        with open(f1_path,'w') as f1:
            f1.write(text_after_mf)
    
    f2_path = os.path.join(folder, 'type_charge.txt')
    with open(f2_path,'w') as f2:
        for type, charge in zip(type_list, charges):
                f2.write(f'{type}-{charge}\n')
                
    
                
    

def print_GauInp(*args):
    ele_list, type_list, coord, bond_type_list, \
    k_bond_list, bond_length_list, angle_type_list, k_angle_list, \
    angle_length_list, torsion_type_list, v1_list, v2_list, \
    v3_list, phase_list, hybrid_list, charges, formal_chg, multi = args

    header_gjf = """%mem=1GB
%nprocshared=1
%chk=SmartField4gau.chk
#p Amber=(SoftFirst,Print) nosymm geom=nocrowd opt(MaxMicroiterations=2000) Freq

Title

{f_chg} {mult}
"""

    master_function = """
!Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
"""
    filename = 'SmartField4gau.gjf'
    with open(filename, 'w') as file_out:
        file_out.write(header_gjf.format(f_chg = formal_chg, mult=multi))
        for element, type, coordinates, charge in zip(ele_list, type_list, coord, charges):
            formatted_coords = '  '.join(f'{x:.6f}' for x in coordinates)
            file_out.write(f'{element}-{type}-{charge}  {formatted_coords}\n')

        file_out.write(master_function)
        file_out.write('!SMARTFIELD FF\n')
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
    file.write('!Torsions\n')
    for torsion_type, phase, v1, v2, v3, hybrid in zip(torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list):
        formatted_phase = '  '.join(f'{x}' for x in phase)
        file.write(f'AmbTrs {torsion_type} {formatted_phase} {v1:.2f} {v2:.2f} {v3:.2f} 0. {float(hybrid)}\n')


def print_AmbFrcmod(*args):
    type_list, bond_type_list, k_bond_list, bond_length_list, \
    angle_type_list, k_angle_list, angle_length_list, torsion_type_list, \
    v1_list, v2_list, v3_list, phase_list, hybrid_list = args

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
        file_out.write('\n')

        write_bonds_amber(file_out, bond_type_list, k_bond_list, bond_length_list)
        write_angles_amber(file_out, angle_type_list, k_angle_list, angle_length_list)
        write_torsions_amber(file_out, torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list)


def write_bonds_amber(file, bond_type_list, k_bond_list, bond_length_list):
    file.write('BOND\n')
    for bond_type, k_bond, bond_length in zip(bond_type_list, k_bond_list, bond_length_list):
        file.write(f'{bond_type} {k_bond:.3f} {bond_length:.3f}\n')
    file.write('\n')


def write_angles_amber(file, angle_type_list, k_angle_list, angle_length_list):
    file.write('ANGLE\n')
    for angle_type, k_angle, angle_length in zip(angle_type_list, k_angle_list, angle_length_list):
        file.write(f'{angle_type} {k_angle:.3f} {angle_length:.3f}\n')
    file.write('\n')
        


def write_torsions_amber(file, torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list):
    file.write('DIHE\n')
    for torsion_type, phase, v1, v2, v3, hybrid in zip(torsion_type_list, phase_list, v1_list, v2_list, v3_list, hybrid_list):
        formatted_phase = '  '.join(f'{x}' for x in phase)
        file.write(f'{torsion_type} {formatted_phase} {v1:.2f} {v2:.2f} {v3:.2f} 0. {hybrid}\n')
    file.write('\n')
    

