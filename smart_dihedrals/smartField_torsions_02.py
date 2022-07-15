#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

sys.path.append('/home/emanuele/chemical_software/proxima/PyProxima')

import pyproxima as pxma

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

def change_mm_COM(filename, angle):
    with open(filename, 'r') as f:
        text = []
        for line in f:
            text.append(line.split())

    # Break it in parts
    first_bit = []
    second_bit = []
    third_bit = []
    for i in text:
        first_bit.append(i)
        if len(i) == 2:
            break
    for i in text[8:]:
        second_bit.append(i)
        if i[0] == 'Variables:':
            break

    print(second_bit[angle[-1]-1][1])
    second_bit[angle[-1]-1][1] = angle[2]  # Replace Bond
    second_bit[angle[-1]-1][3] = angle[1]  #  .... Angle
    second_bit[angle[-1]-1][5] = angle[0]  #  .. Dihedral

    # r1 = second_bit[angle[-1]-1][2]        # Bond Paramters
    # r2 = second_bit[angle[-1]-1][4]        # Angle ....
    r3 = second_bit[angle[-1]-1][6]        #  Dihedral ...

    for i in text[len(second_bit)+8:]:
        third_bit.append(i)

    third_bit = [x for x in third_bit if x]   # Take out empty spaces in list

    for i, x in enumerate(third_bit):
        if x[0] == r3:
            third_bit[i][1] = '0.000 S 15 24.0 '


    fout = filename[:-4] + '_mod.gjf'
    text_new = first_bit + second_bit + third_bit
    with open(fout, 'w') as f2:
        for i in text_new:
            print(*i, file=f2)
        print(" ", file=f2)

    
def change_qm_COM(filename, tors_angle):
    s = " ".join(map(str, tors_angle))

    with open(filename,"r") as f:
        d = f.read()
    f.close()
    m = d.split("\n")
    p = "\n".join(m[:-1])
    # print(m,s)
    with open(filename, 'w+') as f:
        for i in range(len(p)):
            f.write(p[i])
    # f.close
        f.write('D ' + str(s) + ' S 10 36.0')
        f.write('\n')
        f.write('\n')
        f.write('@GAUSS_EXEDIR:dftba.prm')
        f.write('\n')
        f.write('\n')


def read_XYZ(filename):
    with open(filename, 'r') as f:
        data = []
        N_atoms = int(f.readline().split()[0])      
        f.readline()
        for i in range(0, N_atoms):
               data.append(f.readline().split())

        ele = [ i[0] for i in data ]
        XYZ = np.array([ i[1:4] for i in data ], dtype=float)
        return ele, XYZ
    

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
               
def main():
    printout_start()
    file_xyz = sys.argv[1]
    file_top = sys.argv[2]
    ele, _ = read_XYZ(file_xyz)
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
    [print(" ", *i) for i in ele_unique_list]
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

    print(" Original Atomic indices in Torsional Angles:")
    print(" ---------------------------------------------\n")
    [ print(" {}- {}- {}- {} ".format(i,j, k, l)) for i,j,k,l in tors_mean_list ]
    print("")

    # # Call PyProxima 
    parser = pxma.PyParser()
    # mol = parser.readXYZ(file_xyz.encode('utf-8'),b'test')
    mol = parser.readXYZ(file_xyz.encode('utf-8'), b'test')
    mol.computeBonds_Fast(0.2)
    old_serial = []
    for j in range(0, mol.getNumAtoms()):
        old_serial.append(mol.getAtomFromIndex(j).getSerial()+1)
    
    atoms_serial = []
    atoms_indices = []
    for i, x in enumerate(tors_mean_list):
        bnd = mol.getAtom(x[1]-1).findBond(0, x[2]-1)      # Iterate over central bond
        # bnd = mol.getAtom(bond_list[0][0]-1).findBond(0, bond_list[0][1]-1)      # keep current bond
        new_mol = mol.cutAroundBond(bnd, 2, 0)             # Cut bonds

        out = pxma.PyOutputter()
        new_mol.computeTopology(0, 0.4)                    # Compute topology & Add missing Hydrogens
        new_mol = new_mol.addHydrogens(0, False)
        
        serial_tmp = []
        index_tmp = []
        for j in range(0, new_mol.getNumAtoms()):
            serial_tmp.append(new_mol.getAtomFromIndex(j).getSerial()+1)
        index_tmp = [ i for i in range(1,len(serial_tmp)+1) ]
        # print('item :', i, serial_tmp)
        # print('item :', i, index_tmp)
        # print("")

        atoms_serial.append(serial_tmp)
        atoms_indices.append(index_tmp)

        # Adjust indeces of Mean Torsions
        for i in range(len(atoms_serial)):
            if atoms_indices[i] != atoms_serial[i]:
                for idx, item in enumerate(tors_mean_list[i]):
                    if item not in atoms_indices[i]:
                       if item in atoms_serial[i]:
                           index = atoms_serial[i].index(item) 
                           tors_mean_list[i][idx] = atoms_indices[i][index]
        
        print(tors_mean_list[i][3])
    
        file_out_qm = str(i) + '_qm.gjf'
        string_qm_chk = '%chk=' + str(i) + '_qm.chk'
        file_out_mm = str(i) + '_mm.gjf'
        string_mm_chk = '%chk=' + str(i) + '_mm.chk'
        title = 'D: ' + str(x)

        # out.toMol2(new_mol, 0, file_out.encode('utf-8'))
        out.toCOM(new_mol, 0, file_out_qm.encode('utf-8'), \
                        [b'%mem=1GB', \
                         b'%nprocshared=4', \
                        #  b'%chk=file.chk', \
                         string_qm_chk.encode('utf-8'), \
                         b'#p opt=modredundant  nosymm DFTBA'], title.encode('utf-8') )

        # out.toCOM_withTorsion(new_mol, 0, file_out_mm.encode('utf-8'), \
        #                 [b'%mem=1GB', \
        #                  b'%nprocshared=1', \
        #                  string_mm_chk.encode('utf-8'), \
        #                  b'#p opt=Z-matrix  geom=nocrowd nosymm AMBER=HardFisrt'], title.encode('utf-8'), tors_mean_list[i][2])

    
    print(" Fragments Atomic indices in Torsional Angles:")
    print(" ---------------------------------------------\n")
    [ print(" {}- {}- {}- {} ".format(i,j, k, l)) for i,j,k,l in tors_mean_list ]
    print("")



    # Add specifics for running Scan
    for i, x in enumerate(tors_mean_list):
        fqm = str(i) + '_qm.gjf'
        fmm = str(i) + '_mm.gjf'
        change_qm_COM(fqm, x)
        # change_mm_COM(fmm, x)

    

if __name__ == '__main__':
    main()