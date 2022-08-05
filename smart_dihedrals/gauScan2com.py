#!/usr/bin/env python3 

import sys
import numpy as np
import json
import os


def print_mm(elements, xyz, Natm, atom_names, force):

    header = """%mem=10GB
%nprocshared=1
%chk={CHK}
#p Amber=(SoftFirst,Print) geom=nocrowd nosymm 
   
Title
  
0 1
"""
    for i in range(0, len(xyz), Natm):
            j = i + Natm
            fout = str(i) + '.gjf' 
            with open(fout, 'w') as fopen:
                 fopen.write(header.format(CHK=fout[:-3]+'chk'))
                #  fopen.write(f' {Natm}\n')
                #  fopen.write('\n')
                 for m, p, l in zip(elements[i:j], atom_names, xyz[i:j]):
                     s1 = '  '.join(str(x) for x in p)
                     s2 = '  '.join(str(x) for x in l)
                     fopen.write(f'{m}-{s1}      {s2} \n')
                 fopen.write(f'\n')
                 for x in force:
                    l = ' '.join(x)
                    fopen.write(f'{l}\n')
                 fopen.write(f'\n')



def read_optfile(fname):
    """
    Read the option file. expected in json style
    """
    if not os.path.exists(fname):
        raise FileNotFoundError('Missing JSON file')
    with open(fname, 'r')  as fopen:
        data = json.load(fopen)
    for i in ['atom2type', 'force_file', 'log_file']:
        if not os.path.exists(data['files'][i]):
            raise FileNotFoundError(f'Missing {i} file')
    return data

def read_txt_info(fname):
    """
    Parse txt files
    """
    # data = ""
    data = []
    with open(fname, 'r', encoding='UTF-8') as fopen:
        for line in fopen:
            data.append(line.split())
            # data += str(line)
    return data

def read_log(logfile):
    all_lines = []  
    with open (logfile, 'r') as file:                          
        for each_line in file:
            all_lines.append(each_line.strip())
    
    i = 0
    k = 0
    data = []
    for line in range(len(all_lines)):
        if "Optimized Parameters" in all_lines[line]:
            k = k + 1
            start = line
            # print(k, start, all_lines[line])
            for s in range(start, len(all_lines)): 
                j = 0                              
                if "Input orientation:" in all_lines[s]:
                    start2 = s
                    i = i + 1
                    # print(k, i, start2, all_lines[start2])
                    for e in range(start2 + 5, len(all_lines)):
                                  if "----" in all_lines[e]:
                                     end = e
                                     break
                                  j = j + 1
                                  data.append(all_lines[e])
                    break
    
    
    atomic_number = {
        '1':'H', '2':'He', '3':'Li', '4':'Be', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '10':'Ne',
        '11':'Na', '12':'Mg', '13':'Al', '14':'Si','15':'P','16':'S','17':'Cl','18':'Ar','19':'K',
        '20':'Ca', "21" : "Sc" , "22" : "Ti" , "23" : "V"  , "24" : "Cr" , "25" : "Mn", "26" : "Fe" ,
        "27": "Co" , "28" : "Ni" , "29" : "Cu" , "30" : "Zn", "31" : "Ga" , "32" : "Ge" , "33" : "As" ,
        "34": "Se" , "35" : "Br", "36" : "Kr" , "37" : "Rb" , "38" : "Sr" , "39" : "Y"  , "40" : "Zr",
        "41": "Nb" , "42" : "Mo" , "43" : "Tc" , "44" : "Ru" , "45" : "Rh", "46" : "Pd" , "47" : "Ag" ,
        "48": "Cd" , "49" : "In" , "50" : "Sn", "51" : "Sb" , "52" : "Te" , "53" : "I"  , "54" : "Xe" ,
        "55": "Cs", "56" : "Ba" , "57" : "La" , "58" : "Ce" , "59" : "Pr" , "60" : "Nd", "61" : "Pm" ,
        "62": "Sm", "63" : "Eu" , "64" : "Gd" , "65" : "Tb", "66" : "Dy" , "67" : "Ho" , "68" : "Er" ,
        "69": "Tm", "70" : "Yb", "71" : "Lu" , "72" : "Hf" , "73" : "Ta" , "74" : "W"  , "75" : "Re",
        "76": "Os", "77" : "Ir" , "78" : "Pt" , "79" : "Au" , "80" : "Hg", "81" : "Tl" , "82" : "Pb" ,
        "83": "Bi", "84" : "Po" , "85" : "At", "86" : "Rn" , "87" : "Fr" , "88" : "Ra" , "89" : "Ac" ,
        "90": "Th", "91" : "Pa" , "92" : "U"  , "93" : "Np" , "94" : "Pu" , "95" : "Am", "96" : "Cm" ,
        "97": "Bk", "98" : "Cf" , "99" : "Es" ,"100" : "Fm", "101": "Md" ,"102" : "No" ,"103" : "Lr" ,
        "104": "Rf","105" : "Db", "106": "Sg" ,"107" : "Bh" ,"108" : "Hs" ,"109" : "Mt" ,"110" : "Ds",
        "111": "Rg","112" : "Uub","113" : "Uut","114" : "Uuq","115" : "Uup", "116": "Uuh","117" : "Uus","118" : "Uuo"
        }
    
    # Build Ele, XYZ
    Natm = j                                                       
    elements = []
    xyz_angs =[]
    for line in data:
        words = line.split()
        elements.append(atomic_number.get(words[1]))
        xyz_angs.append(words[3:])

    return elements, xyz_angs, Natm
    

if __name__ == "__main__":
  
   fname = sys.argv[1]
   data = read_optfile(fname)
   logfile = str(data["files"]["log_file"]) 
   ele, xyz, Natm = read_log(logfile)
   ffs = read_txt_info(data["files"]["force_file"])          #  Reading in force field
   names = read_txt_info(data["files"]["atom2type"])    #  Reading in atom types
   print_mm(ele, xyz, Natm, names, ffs)



