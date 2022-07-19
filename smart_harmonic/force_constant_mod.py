#!/usr/bin/env python3

import numpy as np
from parser_gau import flat_list

def solve_2Dsys(n1, n2, x, grad, k_tors):
    AM = np.array([[n1 * 0.5* np.sin(n1*x), -n2 *0.5* np.sin(n2*x) ],
                   [n1*n1 * 0.5* np.cos(n1*x), -n2*n2 * 0.5* np.cos(n2*x) ]])
    b = np.array([-grad, k_tors])
    coeffs = np.linalg.solve(AM, b)

    return abs(coeffs)

    

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
        
        

def set_tors(coords, ele_list, type_list,
              tors_list, k_tors, grad, mdout):
    """
    Order Dihedral Force Constants
    """

    center_list = [ x[1:3] for x in tors_list]
    hybrid_list = [ center_list.count(x) for x in center_list]
    v1_eq = []
    v2_eq = []
    for m, p in enumerate(tors_list):
            i = p[0] - 1 
            j = p[1] - 1
            k = p[2] - 1
            l = p[3] - 1
            diff_AB = coords[i,:] - coords[j,:]
            diff_BC = coords[j,:] - coords[k,:]
            diff_CD = coords[k,:] - coords[l,:]
            u_ABC = np.cross(diff_AB, diff_BC)
            u_BCD = np.cross(diff_BC, diff_CD)
            u_ABC_n = np.linalg.norm(u_ABC)
            u_BCD_n = np.linalg.norm(u_BCD)
            cos_phi = np.dot(u_ABC, u_BCD) / ( u_ABC_n * u_BCD_n)
            phi = np.arccos(cos_phi)
            phi_deg = phi * 180 / np.pi
            
            if hybrid_list[m] == 9:
                del_1 = abs(phi_deg - 60.)
                del_2 = abs(phi_deg - 180.)
                eps = 2.0
                if del_1 < eps or del_2 < eps:
                     n = 3.0
                     d = 1.0
                     v1 = abs( -2* (d * k_tors[m])/(n*n* np.cos(n*phi_deg)) )
                     v1_eq.append(v1)
                     print(f' {phi_deg:.2f}  {v1:.2f} {hybrid_list[m]}')
                else:
                    v1, v2 = solve_2Dsys(2, 3, phi, grad[m], k_tors[m])
                    print(f' {phi_deg:.2f} {v1:.2f} {v2:.2f} {hybrid_list[m]}')
            elif hybrid_list[m] == 4:
                del_1 = abs(phi_deg - 0.)
                del_2 = abs(phi_deg - 180.)
                eps = 4.0
                if del_1 < eps or del_2 < eps:
                     n = 2.0
                     d = 1.0
                     v1 = abs( -2* (d * k_tors[m])/(n*n* np.cos(n*phi_deg)) )
                     v1_eq.append(v1)
                     print(f' {phi_deg:.2f}  {v1:.2f} {hybrid_list[m]}')
                else:
                    v1, v2 = solve_2Dsys(1, 2, phi, grad[m], k_tors[m])
                    print(f' {phi_deg:.2f} {v1:.2f} {v2:.2f} {hybrid_list[m]}')


                    # np.linalg.solve(MM_diag, diag_QM)
                    
            # print(f' {phi_deg:.2f}  {v2:.2f} {hybrid_list[m]}')
            
    #     r_BC = np.linalg.norm(diff_BC)
        
    #     u_AB = diff_AB / r_AB
    #     u_BC = diff_BC / r_BC
    #     cos_theta = np.dot(u_AB, u_BC)
    #     theta = np.arccos(cos_theta)
    #     theta = 180 - theta * 180 / np.pi
    #     angle_length_list.append(theta)
    #     angle_type_list.append(type_list[i] + ' ' + type_list[j] +  \
    #                            ' ' + type_list[k]) 
    
    # angle_type_list = flat_list(angle_type_list)

    # # Average over values if duplicates found,
    # # return all of'em
    # if mdout == 'mean':
    #     tmp_arr = np.array(angle_length_list)
    #     angle_length_2d = np.reshape(tmp_arr, ((tmp_arr.shape[0], 1)) )
    #     folded, out_1 = avg_dups(angle_type_list, angle_length_2d)
    #     angle_length_mean = np.reshape(out_1, (out_1.shape[0]))

    #     k_angles_2d = np.reshape(k_angles, ((k_angles.shape[0], 1)) )
    #     _, out_2 = avg_dups(angle_type_list, k_angles_2d)
    #     k_bonds_mean = np.reshape(out_2, (out_2.shape[0]))

    #     return list(folded), angle_length_mean, k_bonds_mean
    # elif mdout == 'all':
    #     return angle_type_list, angle_length_list, k_angles
            