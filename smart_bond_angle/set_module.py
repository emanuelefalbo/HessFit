#!/usr/bin/env python3

import numpy as np
from parser_gau import flat_list

def avg_dups(genes, values):
    """ 
    Average values on the basis of 
    their genes
    """
    folded, indices, counts = np.unique(genes, return_inverse=True, return_counts=True)

    output = np.zeros((folded.shape[0], values.shape[1]))
    np.add.at(output, indices, values)
    output /= counts[:, np.newaxis]

    return folded, output


def set_bonds(coords, ele_list, type_list,
              bond_list, k_bonds, mdout):
    """
    Order Bond Force Constant
    """

    bond_length_list = []
    bond_type_list = []
    for k in bond_list:
        i = k[0] -1 
        j = k[1] -1
        diff_AB = coords[i,:] - coords[j,:]
        r_AB = np.linalg.norm(diff_AB)
        bond_length_list.append(r_AB)
        bond_type_list.append(type_list[i] + ' ' + type_list[j]) 
    
    
    bond_type_list = flat_list(bond_type_list)

    if mdout == 'all':
        return bond_type_list, bond_length_list
        
    elif mdout == 'mean':
        tmp_arr = np.array(bond_length_list)
        bond_length_2d = np.reshape(tmp_arr, ((tmp_arr.shape[0], 1)) )
        folded, out_1 = avg_dups(bond_type_list, bond_length_2d)
        bond_length_mean = np.reshape(out_1, (out_1.shape[0]))

        k_bonds_2d = np.reshape(k_bonds, ((k_bonds.shape[0], 1)) )
        _, out_2 = avg_dups(bond_type_list, k_bonds_2d)
        k_bonds_mean = np.reshape(out_2, (out_2.shape[0]))

        return list(folded), bond_length_mean, k_bonds_mean



def set_angles(coords, ele_list, type_list,
              bond_list, k_bonds):
    """
    Order Bond Force Constant
    """

    # # print(len(type_list))
    # bond_lenght_list = []
    # bond_type_list = []
    # for k in bond_list:
    #     i = k[0] -1 
    #     j = k[1] -1
    #     diff_AB = coords[i,:] - coords[j,:]
    #     r_AB = np.linalg.norm(diff_AB)
    #     bond_lenght_list.append(r_AB)
    #     bond_type_list.append(type_list[i] + ' ' + type_list[j]) 
    
    # bond_type_list = flat_list(bond_type_list)
    # print(bond_type_list)
        
    # for k in range(len(bond_list)):
    #     msg = (
    #           f'HrmStr1 {bond_type_list[k]}  {k_bonds[k]:.3f} ' 
    #           f' {bond_lenght_list[k]:.3f} '
    #           )
    #     print(msg)


    