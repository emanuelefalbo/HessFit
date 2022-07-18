#!/usr/bin/env python3 


import argparse
from ast import arg
import string
import sys
import numpy as np

# def read_gjf(file_gjf):
#     match = "Variables"
#     splitted_line = []
#     with open(file_gjf, 'r' ) as f:
#         for line in f:
#             if line[:9] == match:
#                 break
#         for line in f:
#             splitted_line.append(line.split())
                
#     f.close()
#     for i in splitted_line:
#         for j in i:
#             if j == 'S':
#                 name_ang, start_ang, N_step, size_step = i[0], float(i[1]), int(i[3]), float(i[4])
#     f.close()
#     name_ang_new = name_ang.replace("=", '')
#     return name_ang_new, start_ang, N_step, size_step

def build_parser():
    """
    Build options for parser
    """
    parser = argparse.ArgumentParser()
    txt = 'The type of file to be opened: qm or mm'
    parser.add_argument('-ftype', type=str, default='qm',
                         help=txt)
    txt = "Scan log file"
    parser.add_argument('-file', type=str, help=txt)
    return parser


def read_log(file_log, ftype):
     line_0 = []
     line_1 = []
     splitted_line = []
     match_0 = " The following ModRedundant input section has been read:"
     match_01 = "       Variables:"
     match_1 = " Summary of Optimized Potential Surface Scan"
     match_2 = "Eigenvalues"
     match_3 = " ! Name  Definition              Value          Derivative Info."
     match_4 = " ---"

     # Get Angle step size:
     # If file QM.log search for match_0
     #
     with open(file_log, 'r' ) as f:
        if ftype == 'qm':                            # If log file QM
            for line in f:
                 if line[:56] == match_0:
                     break
            line_0.append(f.readline().split())
            step_ang = float(line_0[0][7])
    
           # Get Starting Angle
            for line in f:
                if line[:64] == match_3:
                    break
            f.readline()
            for line in f:
                if line[:4] == match_4:
                    break
                line_1.append(line.split())
            for i, x in enumerate(line_1):
                if len(x) == 6:
                    angle = float(x[3])
        else:                                   # If log file MM
            for line in f:
                if line[:17] == match_01:
                     break
            for line in f:
                if line == " \n":
                    break
                line_0.append(line.split())
            for i in line_0:
                if len(i) == 5:
                   angle = float(i[1])
                   step_ang = float(i[4])
            # print(angle, step_ang)


        # Get Energies
        for line in f:
            if line[:44] == match_1:
                ipos = line.find('to')
                add = float(line[50:ipos])
                for line in f:
                    ipos = line.find(match_2)
                    if line[ipos:16] == match_2:
                        splitted_line.append(line.split())
     f.close()   

     energy_list = []
     [ energy_list.append(list(map(float,i[2:]))) for i in splitted_line ] 
     energy = [ y for x in energy_list for y in x]
     energy = np.array(energy)
     energy = energy + add
     return step_ang, angle, energy

def main():
    parser = build_parser()
    args = parser.parse_args()

    file_log = args.file
    ftype = args.ftype
    size_step, ang, energy_au = read_log(file_log, ftype)
    file_out = file_log[:-4] + '_scan.dat'
    with open(file_out, 'w') as f:
        for i,x in enumerate(energy_au):
            print("{:,.8e}  {:,.8e}".format(ang, energy_au[i]), file=f)
            ang = ang + size_step
    f.close()


if __name__ == '__main__':
    main()