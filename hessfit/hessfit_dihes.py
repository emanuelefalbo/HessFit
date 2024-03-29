#!/usr/bin/env python3

import subprocess
import os
import glob
import numpy as np
import gauScan2com as scan2mm
import readin_opts as rdin
import log2scan
import pandas as pd

# import gcutil as gc
# from openbabel import openbabel as ob

def printout_start():
    print("""
    =============================================================
    Author : Emanuele Falbo
    E-Mail : falbo.emanuele@gmail.com
    =============================================================
    hessfit_dihes : a program to fragment molecules 
                         on the basis of their dihedral angles
                         prepare and compute inputs for subsequent QM and 
                         MM rigid torsional scans
    =============================================================
    """)

#    B3LYP/6-31G* EmpiricalDispersion=GD3  
def print2QM(filename, data, tors_angle, nprocs, qm_method):
    header = """%mem=1GB
%nprocshared={}
%chk={}.chk
#p nosymm geom=nocrowd opt=(modredundant,maxcycle=100) 
  {}

Title

0 1
"""
    s = " ".join(map(str, tors_angle))
    with open(filename,"w") as f:
        f.write(header.format(nprocs, filename[:-4], qm_method))
        for i in data:
            f.write(f"{i[0]}  {' '.join(map(str, i[1:]))}\n")
        f.write('\n')    
        f.write('D ' + str(s) + ' S 10 36.0\n')
        f.write('\n')    
        

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
        N_bonds = int(f.readline().split()[0])        
        for i in range(0, N_bonds):
            bond_list.append(f.readline().strip().split())
        
        b_list = [ list(map(int,x)) for x in bond_list ]

        N_angles = int(f.readline().split()[0])       
        for i in range(0, N_angles):
            angle_list.append(f.readline().split())  
        a_list = [ list(map(int,x)) for x in angle_list ] 
        
        N_tors = int(f.readline().split()[0]) 
        for i in range(0, N_tors):
            tors_list.append(f.readline().strip().split())
        
        t_list = [ list(map(int,x)) for x in tors_list ]
        # tors_list = np.array([ i[] for i in tors_list ], dtype=float)
        return t_list   
    

def substitute_strings(indices, strings):
    res = []
    for idx_list in indices:
        substituted = [strings[idx - 1] for idx in idx_list]
        res.append(' '.join(substituted))
    return res

def replace_nth_tors(file_mm, tors_type):
    content =[]
    with open(file_mm, 'r') as f:
        content = [line.split() for line in f]
            
    # Find the index of the '!Torsions' entry
    torsions_index = next((index for index, sublist in enumerate(content) if '!Torsions' in sublist), -1)
    content_tors = content[torsions_index + 1:] if torsions_index != -1 else []
    
    # content_list = content.split('\n')
     # split_content = [item.split() for item in content_tors]
    # Split each item of the list into a sublist based on space delimiter
    for line_idx, line in enumerate(content_tors):
        if line[1:5] == tors_type:
           arr_tmp = np.asarray(content_tors[line_idx][9:13])
           arr_float = arr_tmp.astype(float)
           # Replace elements different from 0.0 with 1.0
           arr_float[arr_float != 0.0] = 0.0
        #    print(line[1:5], tors_type, arr_float, content[line_idx][9:13])
           str_type = arr_float.astype(str)
           content_tors[line_idx][9:13] = list(str_type)
           break
    
    # Write content back to the file
    with open(file_mm, 'w') as f:
        for line in content:
            f.write(' '.join(line) + '\n')
    

     
def main():
    printout_start()
    parser = rdin.commandline_parser3()
    opts = parser.parse_args()
    json_opts = rdin.read_optfile_2(opts.optfile)
    file_xyz = json_opts['files']['file_xyz']
    file_top = json_opts['files']['topol']
    file_atype = json_opts['files']['atom2type']
    file_ff_str = json_opts['files']['force_file']
    nprocs = json_opts['nprocs']
    qm_method = json_opts['method']
    
    data = read_XYZ(file_xyz)
    ele = [ i[0] for i in data ]
    # cords = np.array([ i[1:4] for i in data ], dtype=float)
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
    
    # Reading in ff_string and type_charge files
    atype_chg = scan2mm.read_txt_info(file_atype)
    ffs = scan2mm.read_txt_info(file_ff_str)
    atypes = [item[0].split('-')[0] for item in atype_chg]
    
    # Replacing indices with string atype
    tors_tmp = substitute_strings(tors_mean_list, atypes)
    tors_type_list = {idx:x.split(' ') for idx, x in enumerate(tors_tmp)}
    print(" Atomic Type per indices in Torsional Angles:")
    print(" ---------------------------------------------\n")
    [ print(" {}- {}- {}- {} ".format(i,j, k, l)) for i,j,k,l in tors_type_list.values() ]
    print("")
    
    # # Writing QM Scan files for each torsional:
    for id, x in enumerate(tors_mean_list):
        fqm = str(id) + '_qm.gjf'
        # fmm = str(id) + '_zmat_mm.gjf'
        print2QM(fqm, data, x, nprocs, qm_method)
    
    GPATH = os.environ.get("g09root") + "/g09"
    file_pattern = '*_qm.gjf*'
    file_qm_list = sorted(glob.glob(file_pattern))      # the order matters
    
    # Calling Gaussian to perfrom Scan on each QM dihedral
    for id, f in enumerate(file_qm_list):
        print(f"Executing Gaussian Scan on {f}")
        file_base = os.path.splitext(f)[0]
        log_file = file_base + '.log'  
        
    # # Check if the log files exist 
        if os.path.exists(log_file):
            print(f"Log file {log_file} already exists.\n")
        else:
            process = subprocess.Popen([f"{GPATH}/g09", f], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()  # Wait for the process to finish
            if process.returncode != 0:
                print(f"Error executing Gaussian Scan on {f}. Details: {stderr.decode()}")
            else:
                print(f"Gaussian Scan executed successfully on {f}")
                
        #  Creating MM input geometries for each QM file 
        ele, coords, Natoms = scan2mm.read_log(log_file)
        file_mm_list = scan2mm.print_mm(f, ele, coords, Natoms, atype_chg, ffs)
        
        # Silcening n-th dihedral by replacing with 0.0 barrier terms
        replace_item = tors_type_list[id]
        for idf, file in enumerate(file_mm_list):
            replace_nth_tors(file, replace_item)
        
        # Calling Gaussian to perform Optimization on each QM input for full scan 
        for jd, fj in enumerate(file_mm_list):
            process = subprocess.Popen([f"{GPATH}/g09", file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()  # Wait for the process to finish
            if process.returncode != 0:
                print(f"Error executing MM Gaussian Optimization on {file}. Details: {stderr.decode()}")
            else:
                print(f"Gaussian Optimization executed successfully on {file}")

    
    file_pattern = '*_qm.log*'
    qm_logs_list = sorted(glob.glob(file_pattern))      

    for file in qm_logs_list:
        fout_qm = os.path.splitext(file)[0] + "_scan_energy.csv"
        # Gathering energies from QM file
        subprocess.run(["log2scan.py", "-t", "qm", "-f", file, "-o", fout_qm], check=True)

        # Gathering energies from MM files
        fext_mm = file.split("_")[0]
        pattern = '_mm*log'
        files = sorted(glob.glob(fext_mm+pattern)) 
        fout_mm = fext_mm + "_mm_scan_energy.csv"
        subprocess.run(["get_mm_energy.py", "-t", "mm"] + files +["-o", fout_mm], check=True)
        
        df1 = pd.read_csv(fout_qm, names=['col1', 'col2'])
        df2 = pd.read_csv(fout_mm, header=None, usecols=[1], names=['col3'])
        
        merged_df = pd.concat([df1, df2], axis=1)
        f_merged = os.path.splitext(file)[0] + '_all.csv'
        merged_df.to_csv(f_merged, header=None, index=False)
        
        # LLS: Fitting of MM single-point PES to QM relaxes scan
        print(f"LLS fitting on: {file}\n")
        subprocess.run(["fit4dihe.py"] + [f_merged], check=True)
        
        



    
    # Calling Gaussian to perfrom Scan on each QM dihedral
    # for id, f in enumerate(file_qm_list):
    #     print(f"Executing Gaussian Scan on {f}")
    #     file_base = os.path.splitext(f)[0]
    #     log_file = file_base + '.log'  
        
    # # Check if the log files exist 
    #     if os.path.exists(log_file):
    #         print(f"Log file {log_file} already exists.\n")

    
    # Conver to ZMAT by gcutils
    # xyzarr, atomnames = gc.readxyz(xyzfilename)
    # distmat = gc.distance_matrix(xyzarr)
    # gc.write_zmat(xyzarr, distmat, atomnames, rvar=False, avar=False, dvar=True)
    

if __name__ == '__main__':
    main()
    
    
