#!/usr/bin/env python3

import sys
import subprocess
import numpy as np
import os
import glob
import gauScan2com as scan2mm

# import gcutil as gc
# from openbabel import openbabel as ob

def printout_start():
    print("""
    =============================================================
    Author : Emanuele Falbo
    E-Mail : falbo.emanuele@gmail.com
    =============================================================
    SmartField_Torsion : a small program to fragment molecules 
                         on the basis of their dihedral angles
                         and prepare inputs for subsequent QM and 
                         MM rigid torsional scans
    =============================================================
    """)

#    B3LYP/6-31G* EmpiricalDispersion=GD3  
def print2QM(data, filename, tors_angle, nprocs):
    header = """%mem=1GB
%nprocshared={}
%chk={}.chk
#p nosymm geom=nocrowd opt=(modredundant,maxcycle=100) 
  PM6

Title

0 1
"""
    s = " ".join(map(str, tors_angle))
    with open(filename,"w") as f:
        f.write(header.format(nprocs, filename[:-4]))
        for i in data:
            f.write(f"{i[0]}  {' '.join(map(str, i[1:]))}\n")
        f.write('\n')    
        f.write('D ' + str(s) + ' S 10 36.0\n')
        f.write('\n')    
        
    
# def print2MM(data, filename, tors_angle):
#     header = """%mem=1GB
# %nprocshared=1
# %chk={}.chk
# #p nosymm geom=nocrowd opt=(z-matrix,maxcycle=300) Amber=softfirst

# Title

# 0 1
# """
#     s = " ".join(map(str, tors_angle))
#     with open(filename,"w") as f:
#         f.write(header.format(filename[:-4]))
#         for i in data:
#             f.write(f"{i[0]}  {' '.join(map(str, i[1:]))}\n")
#         f.write('\n')    
#         f.write('D ' + str(s) + ' S 10 36.0\n')
#         f.write('\n')   
        
# Define a function to extract variables and their values
        

def read_XYZ(filename):
    with open(filename, 'r') as f:
        data = []
        N_atoms = int(f.readline().split()[0])      
        f.readline()
        for i in range(0, N_atoms):
               data.append(f.readline().split())
        # ele = [ i[0] for i in data ]
        # XYZ = np.array([ i[1:4] for i in data ], dtype=float)
        # return ele, XYZ
        return data
    

def read_topology(filename):
    with open(filename, 'r') as f:
        bond_list = [] 
        angle_list = []
        tors_list = []
        temp = []
        N_bonds = int(f.readline().split()[0])        
        for i in range(0, N_bonds):
            bond_list.append(f.readline().strip().split(','))
        
        bond_list = [ list(map(int,x)) for x in bond_list ]

        N_angles = int(f.readline().split()[0])       
        for i in range(0, N_angles):
            f.readline()
            # angle_list.append(f.readline().split())   

        N_tors = int(f.readline().split()[0])        
        for i in range(0, N_tors):
            tors_list.append(f.readline().strip().split(','))
        
        tors_list = [ list(map(int,x)) for x in tors_list ]
        # tors_list = np.array([ i[] for i in tors_list ], dtype=float)
        return tors_list   
        
# def read_dihe_strings(filename):
#     with open(filename, 'r') as file:
#          lines = file.readlines()
#     found_variables = False
#     dx_strings = []
    
#     for line in lines:
#         if 'Variables:' in line:
#             found_variables = True
#         elif found_variables and line.strip().startswith('d'):
#             # Extract the number after 'd' and check if it's in the range 4 to 10
#             dx_number = line.strip()[1:]
#             dx_strings.append(dx_number)
#     dx_values = [item.split('=')[0].strip() for item in dx_strings]
    
#     return dx_values

# def find_files_matching_pattern(N):
#     file_pattern = "_zmat_mm.gjf"
#     file_list = []

#     for x in range(N + 1):
#         filename = f"{x}{file_pattern}"
#         if os.path.isfile(filename):
#             file_list.append(filename)

#     return file_list


# def header_zmat(file_name):
#     new_header = '''%mem=1GB
# %nprocshared=1
# %chk={}.chk
# #p opt=z-matrix Amber=Softfirst nosymm geom=connectivity\n'''
#     with open(file_name, 'r') as file:
#         lines = file.readlines()
#     lines = lines[2:]
#     lines.insert(0, new_header.format(file_name[:-4]))
#     with open(file_name, 'w') as file:
#         file.writelines(lines)

    
# def change_mm_COM(filename, angle):
#     with open(filename, 'r') as f:
#         text = []
#         for line in f:
#             text.append(line.split())

#     # Break it in parts
#     first_bit = []
#     second_bit = []
#     third_bit = []
#     for i in text:
#         first_bit.append(i)
#         if len(i) == 2:
#             break
#     for i in text[8:]:
#         second_bit.append(i)
#         if i[0] == 'Variables:':
#             break

#     print(second_bit[angle[-1]-1][1])
#     second_bit[angle[-1]-1][1] = angle[2]  # Replace Bond
#     second_bit[angle[-1]-1][3] = angle[1]  #  .... Angle
#     second_bit[angle[-1]-1][5] = angle[0]  #  .. Dihedral
#     # r1 = second_bit[angle[-1]-1][2]        # Bond Paramters
#     # r2 = second_bit[angle[-1]-1][4]        # Angle ....
#     r3 = second_bit[angle[-1]-1][6]        #  Dihedral ...

#     for i in text[len(second_bit)+8:]:
#         third_bit.append(i)

#     third_bit = [x for x in third_bit if x]   # Take out empty spaces in list

#     for i, x in enumerate(third_bit):
#         if x[0] == r3:
#             third_bit[i][1] = '0.000 S 15 24.0 '


#     fout = filename[:-4] + '_mod.gjf'
#     text_new = first_bit + second_bit + third_bit
#     with open(fout, 'w') as f2:
#         for i in text_new:
#             print(*i, file=f2)
#         print(" ", file=f2)
     
def main():
    printout_start()
    file_xyz = sys.argv[1]
    file_top = sys.argv[2]
    nprocs = sys.argv[3]
    # ele, coords = read_XYZ(file_xyz)
    data = read_XYZ(file_xyz)
    ele = [ i[0] for i in data ]
    cords = np.array([ i[1:4] for i in data ], dtype=float)
    tors_list = read_topology(file_top)

    # Building torsional list by element type
    ele_list = []
    for id_1, x in enumerate(tors_list):
        tmp_list = []
        for id_2, y in enumerate(x):
            tmp_list.append(ele[y-1])
        ele_list.append(tmp_list)

    # Explicit tors_unique by element strings
    ele_unique = tuple(set([tuple(i) for i in ele_list]))
    ele_unique_list = [ list(i) for i in ele_unique]
    print(" Torsional Angles by Elements:")
    print(" ----------------------------\n")
    [print(i, ": ", *j) for i, j in enumerate(ele_unique_list)]
    print("")
    
    tmp_center = tuple(set([tuple(i[1:3]) for i in tors_list]))
    print(" Rotation Around Bonds:")
    print(" ----------------------\n") 
    [ print(" {}- {} ".format(i,j)) for i, j in tmp_center ]
    print("")

    # Explicit tors_unique 
    tors_mean_list = []
    for i, x in enumerate(ele_unique_list, start=0):
        for j, y in enumerate(ele_list, start=0):
            if x == y:
               tors_mean_list.append(tors_list[j])
               break

    print(" Atomic indices in Torsional Angles:")
    print(" ---------------------------------------------\n")
    [ print(" {}- {}- {}- {} ".format(i,j, k, l)) for i,j,k,l in tors_mean_list ]
    print("")

    # Writing QM Scan files for each torsional:
    for id, x in enumerate(tors_mean_list):
        fqm = str(id) + '_zmat_qm.gjf'
        # fmm = str(id) + '_zmat_mm.gjf'
        print2QM(data, fqm, x, nprocs)
        
    # Reading in ff_string and type_charge files
    atom_types = scan2mm.read_txt_info('type_charge.txt')
    ffs = scan2mm.read_txt_info('ff_string.txt')
    
    GPATH = os.environ.get("g09root") + "/g09"
    file_pattern = '*_zmat_qm.gjf*'
    file_qm_list = glob.glob(file_pattern)      # the order doesn't matter
    # Calling Gaussian to perfrom Scan on each QM dihedral
    for id, f in enumerate(file_qm_list[:2]):
        print(f"Executing Gaussian Scan on {f}")
        file_base = os.path.splitext(f)[0]
        log_file = file_base + '.log'  
    # Check if the log files exist 
        if os.path.exists(log_file):
            print(f"Log file {log_file} already exists.\n")
        else:
            process = subprocess.Popen([f"{GPATH}/g09", f], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()  # Wait for the process to finish
            if process.returncode != 0:
                print(f"Error executing Gaussian Scan on {f}. Details: {stderr.decode()}")
            else:
                print(f"Gaussian Scan executed successfully on {f}")
        # Creating MM input geometries for each QM file 
        ele, coords, Natoms = scan2mm.read_log(log_file)
        file_mm_list = scan2mm.print_mm(f, ele, coords, Natoms, atom_types, ffs)
         # Calling Gaussian to perform Optimization on each QM input for full scan 
        for jd, fj in enumerate(file_mm_list):
            process = subprocess.Popen([f"{GPATH}/g09", fj], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()  # Wait for the process to finish
            if process.returncode != 0:
                print(f"Error executing MM Gaussian Optimization on {fj}. Details: {stderr.decode()}")
            else:
                print(f"Gaussian Optimization executed successfully on {fj}")
    
    
        # conv = ob.OBConversion()
        # conv.SetInAndOutFormats("xyz", "gzmat")
        # mol = ob.OBMol()
        # conv.ReadFile(mol, file_xyz)
        # # conv.WriteFile(mol, fqm)
        # conv.WriteFile(mol, fmm)
        # Add the header for Scan
        # header_zmat(fmm)
        # header_zmat(fqm)
    # # Replace N with the maximum value you want to search for (e.g., 10)
    # N = len(tors_mean_list)
    # files_mm = find_files_matching_pattern(N)
    # dx_strings = read_dihe_strings(files_mm[0])
    # print(dx_strings, files_mm)
    
    # Conver to ZMAT by gcutils
    # xyzarr, atomnames = gc.readxyz(xyzfilename)
    # distmat = gc.distance_matrix(xyzarr)
    # gc.write_zmat(xyzarr, distmat, atomnames, rvar=False, avar=False, dvar=True)


    

if __name__ == '__main__':
    main()
    
    