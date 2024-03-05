#!/usr/bin/env python3

import numpy as np
import os

def range2(start, end):
     return range(start, end+1)

def flat_list(lis):
    flatList = []
    # Iterate with outer list
    for element in lis:
        if type(element) is list:
            # Check if type is list than iterate through the sublist
            for item in element:
                flatList.append(item)
        else:
            flatList.append(element)
    return flatList
 
def store_any_file(fname):
    """ 
    Reading Any file 
    """
    with open(fname, 'r') as f:
        all_lines = []
        for line in f:
            all_lines.append(line.strip())

    return all_lines

def read_XYZ(all_lines):
    """ 
    Reading CC XYZ from fchk file 
    """
    # with open(fname, 'r') as f:
    #     all_lines = []
    #     for line in f:
    #         all_lines.append(line.strip())

    for s in range(len(all_lines)):                          # Get no of Atoms
        if 'Number of atoms' in all_lines[s]:
            Natom = int(all_lines[s][-10:])  
        
    No_cc = 3* Natom
    nlines = int(No_cc/5.0) + 1
    cc_xyz_list = []
    for s in range(len(all_lines)):                          #  Get CC Coordinates
        if 'Current cartesian coordinates' in all_lines[s]:                              #Reads the Input orientation information
            start = s
            for e in range2(start + 1, start + nlines):
                cc_xyz_list.append(all_lines[e].split() )
    
    cc_xyz_flat = flat_list(cc_xyz_list)

    # Take out Garbage Strings
    if len(cc_xyz_flat) != No_cc:
        diff = abs(len(cc_xyz_flat) - No_cc)
        cc_xyz_flat_mod = cc_xyz_flat[:-diff]
        cc_xyz_1D = np.array(cc_xyz_flat_mod, float)
    else:
        cc_xyz_1D = np.array(cc_xyz_flat, float)

    cc_xyz_arr = np.reshape(cc_xyz_1D, (Natom,3))
    cc_xyz_arr = cc_xyz_arr*0.529177                          # Convert from Bohr to Ang.
    
    return Natom, cc_xyz_arr

def read_RicDim_Grad(all_lines):
    """ 
    Reading Redundant Internal Dimension &
    Gradient from fchk file 
    """

    for s in range(len(all_lines)):                          # Get no of Atoms
        if 'Redundant internal dimensions' in all_lines[s]:
            tmp_list = all_lines[s+1].split()
            ric_list = list(map(int,tmp_list))

    force_list = []
    for s in range(len(all_lines)):                          # Get no of Atoms
        if 'Internal Forces' in all_lines[s]:
            N_force = int(all_lines[s][-10:])
            nlines = int(N_force/5.0) + 1
            start = s
            for e in range2(start + 1, start + nlines):
                force_list.append(all_lines[e].split() )
    
    
    force_flat = flat_list(force_list)

    # Take out Garbage Strings
    if len(force_flat) != N_force:
        diff = abs(len(force_flat) - N_force)
        force_flat_mod = force_flat[:-diff]
        force_1D = np.array(force_flat_mod, float)
    else:
        force_1D = np.array(force_flat,float)

    return ric_list, force_1D


def symmetricize(arr1D):
    ID = np.arange(arr1D.size)
    return arr1D[np.abs(ID - ID[:,None])]


def fill_lower_diag(a):
    n = int(np.sqrt(len(a)*2))+1
    mask = np.tri(n,dtype=bool, k=-1) # or np.arange(n)[:,None] > np.arange(n)
    out = np.zeros((n,n),dtype=int)
    out[mask] = a
    return out

def read_HessRIC(all_lines, ric_list):
    """ 
    Reading Internal Hessian from fchk file:
    all_lines : text file
    ric_list = [RICs, BONDs, ANGLEs, DIHEs]  
    """

    hess_list = []
    for s in range(len(all_lines)):                          # Get no of Atoms
        if 'Internal Force Constants' in all_lines[s]:
            N_hess = int(all_lines[s][-10:])
            nlines = int(N_hess/5.0) + 1
            start = s
            for e in range2(start + 1, start + nlines):
                hess_list.append(all_lines[e].split() )
    
    hess_flat = flat_list(hess_list)

    # Take out Garbage Strings
    if len(hess_flat) != N_hess:
        diff = abs(len(hess_flat)-N_hess)
        # print('diff = ', diff)
        hess_flat_mod = hess_flat[:-diff]
        hess_1D = np.array(hess_flat_mod, float)
    else:
        hess_1D = np.array(hess_flat,float)

    
    # hess_1D_mod = np.append(hess_1D, [0])
    # print(hess_1D_mod)

    # print(len(hess_1D))
    ric_len = ric_list[0]
    # print(ric_len)
    hess_arr = np.zeros((ric_len, ric_len))
    # print(hess_arr.shape, hess_arr.size)
    for i in range(ric_len+1):
        # i_low = int( 0.5 * i * (i - 1) + -1 )           # Adding n(n+1)/2 elements in up/low triangular matrix
        i_low = int( 0.5 * i * (i - 1)  )           # Adding n(n+1)/2 elements in up/low triangular matrix
        # i_up = int( 0.5 * i *  (i + 1) + 1 )
        i_up = int( 0.5 * i *  (i + 1)   )
        # print(i, i_low, i_up, hess_1D[i_low:i_up])
        hess_arr[i-1,0:i] =  hess_1D[i_low: i_up]
        hess_arr[0:i,i-1] =  hess_1D[i_low: i_up]

    return hess_arr

def read_HessXYZ(all_lines, N_atom):
    """ 
    Reading XYZ Hessian from fchk file:
    all_lines : text file
    """

    hess_list = []
    for s in range(len(all_lines)):                          # Get no of Atoms
        if 'Cartesian Force Constants' in all_lines[s]:
            N_hess = int(all_lines[s][-10:])
            nlines = int(N_hess/5.0) + 1
            start = s
            for e in range2(start + 1, start + nlines):
                hess_list.append(all_lines[e].split() )
    
    hess_flat = flat_list(hess_list)

    # Take out Garbage Strings
    if len(hess_flat) != N_hess:
        diff = abs(len(hess_flat)-N_hess)
        # print('diff = ', diff)
        hess_flat_mod = hess_flat[:-diff]
        hess_1D = np.array(hess_flat_mod, float)
    else:
        hess_1D = np.array(hess_flat,float)

    hess_1D = hess_1D * ((627.509391)/(0.529117*0.529117))   # From Hartree/bohr to kcal/mol / angstrom
    # hess_1D_mod = np.append(hess_1D, [0])
    # print(hess_1D_mod)

    # print(len(hess_1D))
    len_hess = 3*N_atom
    hess_XYZ = np.zeros((len_hess, len_hess))
    # print(hess_arr.shape, hess_arr.size)
    for i in range(len_hess + 1):
        # i_low = int( 0.5 * i * (i - 1) + -1 )        # Adding n(n+1)/2 elements in up/low triangular matrix
        i_low = int( 0.5 * i * (i - 1)  )              # Adding n(n+1)/2 elements in up/low triangular matrix
        # i_up = int( 0.5 * i *  (i + 1) + 1 )
        i_up = int( 0.5 * i *  (i + 1)   )
        # print(i, i_low, i_up, hess_1D[i_low:i_up])
        hess_XYZ[i-1,0:i] =  hess_1D[i_low: i_up]
        hess_XYZ[0:i,i-1] =  hess_1D[i_low: i_up]

    return hess_XYZ


def read_NamesTypes(all_lines):
    """ 
    Reading Atom Names from log file:
    """
    ele_list = []
    atype_list = []
    if all_lines is not None and len(all_lines) > 0: 
        for item in all_lines:
                if "-" in item:
                    ele_list.append(item.split("-")[0])
                    atype_list.append(item.split("-")[1])
                else:
                    ele_list.append(item)
                    # print(item.split("-")[0])
    
    # Build atype_list from scratch if is null
    if not atype_list:
        for i, var in enumerate(ele_list):
        # if not atype_list:
            atype_list.append(''.join(f'{var}{i}'))
        #    atype_list[i] = ''.join(f'{ele_list[i]}{i}')
        
    ele_list = flat_list(ele_list)
    atype_list = flat_list(atype_list)

    return ele_list, atype_list
    

def read_Top(all_lines, ric_list):
    """ 
    Reading Topology from RIC:
    """
    tmp_list = []
    ric_tot = ric_list[0]
    for s in range(len(all_lines)):                          # Get no of Atoms
        if 'Name' in all_lines[s]:
            start = s
            # print(all_lines[s])
            for e in range(start+2, start + ric_tot+2):
                tmp_list.append(all_lines[e][8:25].strip())
            break
   

    n_bonds = 0
    n_angles = 0
    n_dihe = 0
    bond_list = []
    angle_list = []
    dihe_list = []
    for i in range(len(tmp_list)):
        if tmp_list[i][0] == 'R':
           n_bonds += 1 
           bond_list.append(tmp_list[i][2 : -1].replace(','," ").split())
        elif tmp_list[i][0] == 'A':
           n_angles += 1 
           angle_list.append(tmp_list[i][2 : -1].replace(','," ").split())
        elif tmp_list[i][0] == 'D':
           n_dihe += 1 
           dihe_list.append(tmp_list[i][2 : -1].replace(','," ").split())

    bond_list = [list(map(int, x)) for x in bond_list]
    angle_list = [list(map(int, x)) for x in angle_list]
    dihe_list = [list(map(int, x)) for x in dihe_list]

    return bond_list, angle_list, dihe_list


def read_CM5(all_lines, N_atoms):
         """ 
         Reading CM5 charges from log file 
         """
         match = 'Hirshfeld charges, spin densities, dipoles, and CM5 charges using IRadAn=      5:'
         match_2 = 'Mulliken charges:'
         text = []

         if match in all_lines:
             for s in range(len(all_lines)):                            
                 if match in all_lines[s]:                              
                     start = s
                     for e in range2(start + 2, start + N_atoms+1):
                         text.append(all_lines[e].split())
                     break
             chg = [ x[7] for x in text ]
         else:                                                       # If CM5 not found, using Mulliken
             for s in range(len(all_lines)):                            
                 if match_2 in all_lines[s]:                            
                     start = s
                     for e in range2(start + 2, start + N_atoms+1):
                         text.append(all_lines[e].split())
                     break
             
             chg = [ x[2] for x in text ]

         # Add '+' to positive charge
         chg_mod = []
         for x in chg:
             if x[0] != '-':
                tmp = '+' + x
             else:
                tmp = x
             chg_mod.append(tmp)

         return chg_mod

def read_AmberParm(path, type_list):
    fname = path + '/amber.prm'
    with open(fname, 'r') as f:
       all_lines = []
       for line in f:
           all_lines.append(line.strip()) 

    Amber_list = []
    for s in range(len(all_lines)):                          #  Get CC Coordinates
        if 'VDW' in all_lines[s]:                              #Reads the Input orientation information
            # start = s
            # for e in range(start + 1, len(all_lines)):
            Amber_list.append(all_lines[s].split() )

    VDW_list = []
    for i in type_list:
        for j in Amber_list:
            if i == j[1]:
                VDW_list.append(j)

    VDW_new = []
    for i in VDW_list:
        if i not in VDW_new:
            VDW_new.append(i)
        
    return VDW_new


def read_mass(all_lines, N_atom, mode):
    """ 
    Reading Masses from fchk file:
    all_lines : text file
    """
    
    # semi_last_column = []
    
    # # Iterate through lines until encountering "Temperature" string
    # for line in lines:
    #     if "Temperature" in line:
    #         break
    #     words = line.split()
    #     if len(words) > 0:
    #         semi_last_column.append(words[-2])
    
    # # Display the semi-last column strings
    # print(semi_last_column)
    
    if mode == 'qm':
        mass_list = []
        for s in range(len(all_lines)):                          # Get no of Atoms
            if 'Real atomic weights' in all_lines[s]:
                N_mass = int(all_lines[s][-10:])
                nlines = int(N_mass/5.0) + 1
                start = s
                for e in range2(start + 1, start + nlines):
                    mass_list.append(all_lines[e].split() )
        
        mass_flat = flat_list(mass_list)
        
    
        # Take out Garbage Strings
        if len(mass_flat) != N_mass:
            diff = abs(len(mass_flat)-N_mass)
            # print('diff = ', diff)
            hess_flat_mod = mass_flat[:-diff]
            mass_1D = np.array(hess_flat_mod, float)
        else:
            mass_1D = np.array(mass_flat,float)
        return mass_1D
    elif mode =='mm':
        mass_col = []
        start_capture = False
        # Iterate through lines between 'Temperature' and 'Molecular mass'
        for line in all_lines:
            if "Molecular mass" in line:
                break
            if start_capture:
                words = line.split()
                if len(words) > 0:
                    mass_col.append(words[-1])
            if "Temperature" in line:
                start_capture = True
        mass_1D = np.asarray(mass_col,float)
        return mass_1D



    

   

