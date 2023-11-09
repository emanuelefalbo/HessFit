#!/usr/bin/env python3

import sys

# Reading in Force Field file
#filename=input("Enter log file to read\n")
filename = sys.argv[1]
#if len(sys.argv) > 1:
#   sys.exit()
with open(filename, "r") as f_log:
#with open('test06.log', "r") as f_log:
     #Initialize values
     tmp = []
     #f_log.readline()   # This is starting of MASS
     # Reading in File
     match_start = "! Name  Definition              Value          Derivative Info.                !"
     #match_end = "Trust"
     match_end = "--------------------------------------------------------------------------------"
     for line in f_log:
         if match_start in line:                                 # Break loop if this string is found
            break
     f_log.readline()
     for line in f_log:
         if match_end in line:
             break
         string = line.strip().split()
         tmp.append(string)
f_log.close() 

# Remove Characters from Temporary list 
temp = []
bond_ij = []
angle_ijk = []
dihedral_ijkl = []
for count, value in enumerate(tmp):
#    print(value[2])
    string = value[2]
    string = string[2:len(string)-1] 
    temp.append(string.split(','))

# Build list of Bonds and Angles from temp list
for i in range(len(temp)):
   # for j in range(len(temp[i])):
    if len(temp[i]) <= 2:
       bond_ij.append(temp[i])
    elif len(temp[i]) <= 3:
       angle_ijk.append(temp[i])
    elif len(temp[i]) <= 4:
       dihedral_ijkl.append(temp[i])

# Print into a Bond_Angle_list file
with open('bond_angle_list.txt', 'w') as f_out:
     print(len(bond_ij), file=f_out)
     for count, value in enumerate(bond_ij):
         print(*bond_ij[count], sep = ' ', file = f_out)    
     print(len(angle_ijk), file=f_out)
     for count, value in enumerate(angle_ijk):
         print(*angle_ijk[count], sep = ' ', file = f_out)    

#print(dihedral_ijkl)    
f_out.close()
          
    


