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
    Order Bond Force Constants
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

    # Average over values if duplicates found,
    # return all of'em
    if mdout == 'mean':
        tmp_arr = np.array(bond_length_list)
        bond_length_2d = np.reshape(tmp_arr, ((tmp_arr.shape[0], 1)) )
        folded, out_1 = avg_dups(bond_type_list, bond_length_2d)
        bond_length_mean = np.reshape(out_1, (out_1.shape[0]))

        k_bonds_2d = np.reshape(k_bonds, ((k_bonds.shape[0], 1)) )
        _, out_2 = avg_dups(bond_type_list, k_bonds_2d)
        k_bonds_mean = np.reshape(out_2, (out_2.shape[0]))

        return list(folded), bond_length_mean, k_bonds_mean
    elif mdout == 'all':
        return bond_type_list, bond_length_list, k_bonds
        
def set_angles(coords, ele_list, type_list,
              angle_list, k_angles, mdout):
    """
    Order Angles Force Constants
    """

    angle_length_list = []
    angle_type_list = []
    for p in angle_list:
        i = p[0] - 1 
        j = p[1] - 1
        k = p[2] - 1
        diff_AB = coords[i,:] - coords[j,:]
        diff_BC = coords[j,:] - coords[k,:]
        r_AB = np.linalg.norm(diff_AB)
        r_BC = np.linalg.norm(diff_BC)
        
        u_AB = diff_AB / r_AB
        u_BC = diff_BC / r_BC
        cos_theta = np.dot(u_AB, u_BC)
        theta = np.arccos(cos_theta)
        theta = 180 - theta * 180 / np.pi
        angle_length_list.append(theta)
        angle_type_list.append(type_list[i] + ' ' + type_list[j] +  \
                               ' ' + type_list[k]) 
    
    angle_type_list = flat_list(angle_type_list)

    # # Average over values if duplicates found,
    # # return all of'em
    if mdout == 'mean':
        tmp_arr = np.array(angle_length_list)
        angle_length_2d = np.reshape(tmp_arr, ((tmp_arr.shape[0], 1)) )
        folded, out_1 = avg_dups(angle_type_list, angle_length_2d)
        angle_length_mean = np.reshape(out_1, (out_1.shape[0]))

        k_angles_2d = np.reshape(k_angles, ((k_angles.shape[0], 1)) )
        _, out_2 = avg_dups(angle_type_list, k_angles_2d)
        k_bonds_mean = np.reshape(out_2, (out_2.shape[0]))

        return list(folded), angle_length_mean, k_bonds_mean
    elif mdout == 'all':
        return angle_type_list, angle_length_list, k_angles
        
        