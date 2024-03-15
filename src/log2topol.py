#!/usr/bin/env python3

import argparse
import os


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def commandline_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f1','--log_file', help='Gaussian QM log file ')
    
    return parser
 
def store_any_file(fname):
    """ 
    Reading Any file 
    """
    with open(fname, 'r') as f:
        all_lines = []
        for line in f:
            all_lines.append(line.strip())

    return all_lines

def read_Top(all_lines):
    """ 
    Reading Topology from RIC:
    """
    
    count = 0
    match = '-'*80
    # Get the Number of RICs
    for s in range(len(all_lines)):                          
        if all_lines[s] == match:
            start = s
            break
    for e in range(start+1, len(all_lines)):                          
        if all_lines[e] == match:
            break
        count += 1
        
    tmp_list = []
    for s in range(len(all_lines)):                          # Get no of Atoms
        if 'Name' in all_lines[s]:
            start = s
            # print(all_lines[s])
            for e in range(start+2, start + count +2):
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

def print_topol(*args):
    bl = args[0]
    al = args[1]
    tl = args[2]
    fname = 'topol.txt'
    with open(fname, 'w') as fout:
        fout.write(f'{len(bl)}\n')
        for i in bl:
            fout.write(f' {i[0]} {i[1]}\n')
        fout.write(f'{len(al)}\n')
        for i in al:
            fout.write(f' {i[0]} {i[1]} {i[2]}\n')
        fout.write(f'{len(tl)}\n')
        for i in tl:
            fout.write(f' {i[0]} {i[1]} {i[2]} {i[3]}\n')
            
            

if __name__ == '__main__':
    parser = commandline_parser()
    opts = parser.parse_args()
    f_qm_log = opts.log_file
    text_qm_log = store_any_file(f_qm_log)
    
    # # Reading in Topology in RIC from log file
    bond_list, angle_list, tors_list = read_Top(text_qm_log)
    print_topol(bond_list, angle_list, tors_list)
    
    
    


