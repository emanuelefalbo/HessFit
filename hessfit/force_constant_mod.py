#!/usr/bin/env python3

from math import floor
import numpy as np
from parser_gau import flat_list
import seminario_module as sem_mod

def solve_2Dsys(n1, n2, x, grad, k_tors):
    AM = np.array([[n1 * 0.5* np.sin(n1*x), -n2 *0.5* np.sin(n2*x) ],
                   [n1*n1 * 0.5* np.cos(n1*x), -n2*n2 * 0.5* np.cos(n2*x) ]])
    b = np.array([-grad, k_tors])
    coeffs = np.linalg.solve(AM, b)

    return abs(coeffs)

def flatList_to_2Darray(my_list):
    tmp_arr = np.array(my_list)
    array_2d = np.reshape(tmp_arr, ((tmp_arr.shape[0], 1)) )
    return array_2d
    

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


def set_bonds(coords, hess, type_list,
              bond_list, k_bonds, mdin, mdout):
    """
    Order Bond Force Constants
    """
    bond_length_list = []
    bond_type_list = []
    for m, k in enumerate(bond_list):
        i = k[0] -1 
        j = k[1] -1
        diff_AB = coords[i,:] - coords[j,:]
        r_AB = np.linalg.norm(diff_AB)
        bond_length_list.append(r_AB)
        bond_type_list.append(type_list[i] + ' ' + type_list[j]) 
        if mdin == 'sem':            # Triggers on Seminario 
            k_bonds[m] = sem_mod.get_ModSem_FcBonds(i, j, diff_AB, r_AB, hess)

    bond_type_list = flat_list(bond_type_list)

    # Average over values if duplicates found,
    # return all of'em
    if mdout == 'mean':
        tmp_arr = np.array(bond_length_list)
        bond_length_2d = np.reshape(tmp_arr, ((tmp_arr.shape[0], 1)) )
        folded, out_1 = avg_dups(bond_type_list, bond_length_2d)
        bond_type_unique = list(folded)
        bond_length_mean = np.reshape(out_1, (out_1.shape[0]))

        k_bonds_2d = np.reshape(k_bonds, ((k_bonds.shape[0], 1)) )
        _, out_2 = avg_dups(bond_type_list, k_bonds_2d)
        k_bonds_mean = np.reshape(out_2, (out_2.shape[0]))
        return bond_type_unique, bond_length_mean, k_bonds_mean
    elif mdout == 'all':
        return bond_type_list, bond_length_list, k_bonds
        
def set_angles(coords, hess, type_list,
              angle_list, k_angles, mdin, mdout):
    """
    Order Angles Force Constants
    """
    angle_length_list = []
    angle_type_list = []
    for m, p in enumerate(angle_list):
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
        if mdin == 'sem':            # Triggers on ModSeminario 
           k_angles[m] = sem_mod.get_ModSem_FcAngles(i, j, k, r_AB, r_BC, \
                                                     u_AB, u_BC, hess)

    angle_type_list = flat_list(angle_type_list)

    # # Average over values if duplicates found,
    # # return all of'em
    if mdout == 'mean':
        tmp_arr = np.array(angle_length_list)
        angle_length_2d = np.reshape(tmp_arr, ((tmp_arr.shape[0], 1)) )
        folded, out_1 = avg_dups(angle_type_list, angle_length_2d)
        angle_type_unique = list(folded)
        angle_length_mean = np.reshape(out_1, (out_1.shape[0]))

        k_angles_2d = np.reshape(k_angles, ((k_angles.shape[0], 1)) )
        _, out_2 = avg_dups(angle_type_list, k_angles_2d)
        k_bonds_mean = np.reshape(out_2, (out_2.shape[0]))

        return angle_type_unique, angle_length_mean, k_bonds_mean
    elif mdout == 'all':
        return angle_type_list, angle_length_list, k_angles
        
        

def set_torsion(coords, type_list, tors_list, \
                k_tors, grad, mdout):
    """
    Order Dihedral Force Constants
    """
    center_list = [ x[1:3] for x in tors_list]
    hybrid_list = [ center_list.count(x) for x in center_list]
    v1_eq = [0] * len(tors_list)
    v2_eq = [0] * len(tors_list)
    v3_eq = [0] * len(tors_list)
    tors_length_list = []
    tors_type_list = []
    phase = np.zeros((len(tors_list), 4))
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
            # Floor cos_phi if out of [-1,1] range
            if cos_phi > 1.0 or cos_phi < -1.0:
                cos_phi = float(floor(cos_phi))
            phi = np.arccos(cos_phi)
            phi_deg = phi * 180 / np.pi
            tors_length_list.append(phi_deg)
            my_str = type_list[i] + ' ' + type_list[j] + \
                     ' ' + type_list[k] + ' ' + type_list[l]
            tors_type_list.append(my_str) 
            
            # Compute Torsional Force Constant 
            # on the basis  of torsion periodicity: sp3,sp2,sp
            if hybrid_list[m] == 9:           
                del_1 = abs(phi_deg - 60.)
                del_2 = abs(phi_deg - 180.)
                eps = 5.0
                if del_1 < eps or del_2 < eps:
                     n, d = 3.0, 1.0
                     v3 = abs( - (d * np.abs(k_tors[m]))/(n*n* np.cos(n*phi_deg)) )
                    #  print(k_tors[m], v3, n)
                     if v3 > 5.0:                         # If force constants too stiff
                        v3_eq[m] = np.exp(-1.4/v3)*1.4     # use AMBER X-C-C-X values
                     else:
                        v3_eq[m] = v3
                    #  print(f' {phi_deg:.2f}  {v3_eq[m]:.2f} {hybrid_list[m]}')
                else:
                    v2, v3 = solve_2Dsys(2, 3, phi, grad[m], k_tors[m])
                    if v2 > 5.0 or v3 > 5.0:
                        v2_eq[m] = 0.0
                        v3_eq[m] = 1.4
                    else:
                        v2_eq[m] = v2
                        v3_eq[m] = v3
                    # print(f' {phi_deg:.2f} {v2_eq[m]:.2f} {v3_eq[m]:.2f} {hybrid_list[m]}')
            elif hybrid_list[m] == 4 or hybrid_list[m] == 2:          # sp2-sp2 bonds
                del_1 = abs(phi_deg - 0.)
                del_2 = abs(phi_deg - 180.)
                eps = 5.0
                phase[m, 1] = 180
                if del_1 < eps or del_2 < eps:
                     n, d = 2.0, 1.0
                     v2 = abs( (d * np.abs(k_tors[m]))/(n*n* np.cos(n*phi_deg)) )
                    #  print(k_tors[m], v2, n)
                    #  if v2 > 30 or 0 < v2 < 9.0:
                     if v2 > 30:
                        v2_eq[m] = np.exp(-30./v2)*v2       # Reduce to a full C=C
                        # print(v2_eq[m])
                        # print("")
                     elif 0 < v2 < 14.5:
                          v2_eq[m] = np.exp(-v2/30.)*30. - np.exp(-v2/14.5)*14.5   # Scale for full C=C & 
                     else:                                 # Increase to a partial C=C
                        v2_eq[m] = v2
                    #  print(f'{phase[m, 1]}  {v2:.2f} {hybrid_list[m]}')
                else:
                    v1, v2 = solve_2Dsys(1, 2, phi, grad[m], k_tors[m])
                    if v1 > 30 or v2 > 30:
                        v1_eq[m] = 0.0
                        v2_eq[m] = np.exp(-v2/30.)*30.
                    else:
                        v1_eq[m] = v1
                        v2_eq[m] = v2
                    # print(f' {phi_deg:.2f} {v1:.2f} {v2:.2f} {hybrid_list[m]}')
            else:
                n, d = 2.0, 1.0
                v2 = abs( -2* (d * k_tors[m])/(n*n* np.cos(n*phi_deg)) )
                if v2 > 30:
                   v2_eq[m] = 14.5
                else:
                   v2_eq[m] = v2
                # print(f' {phi_deg:.2f} {v2:.2f} {hybrid_list[m]}')
    
    tors_type_list = flat_list(tors_type_list)
    phase = phase.astype(int)
    
    # for i in range(len(tors_list)):
    #     print(f'{hybrid_list[i]} {v1_eq[i]:.1f}  {v2_eq[i]:.1f}  {v3_eq[i]:.1f}')
    # # Average over values if duplicates found,
    # # return all of'em
    if mdout == 'mean':
        tors_length_2d = flatList_to_2Darray(tors_length_list)
        folded, out_0 = avg_dups(tors_type_list, tors_length_2d)
        tors_type_unique = list(folded)
        tors_length_mean = np.reshape(out_0, (out_0.shape[0]))

        hybrid_2d = flatList_to_2Darray(hybrid_list)
        _, out_1 = avg_dups(tors_type_list, hybrid_2d)
        tors_type_unique = list(folded)
        hybrid_mean = np.reshape(out_1, (out_1.shape[0]))

        v1_eq_2d = flatList_to_2Darray(v1_eq)
        _, out_1 = avg_dups(tors_type_list, v1_eq_2d)
        v1_eq_mean = np.reshape(out_1, (out_1.shape[0]))

        v2_eq_2d = flatList_to_2Darray(v2_eq)
        _, out_2 = avg_dups(tors_type_list, v2_eq_2d)
        v2_eq_mean = np.reshape(out_2, (out_2.shape[0]))

        v3_eq_2d = flatList_to_2Darray(v3_eq)
        _, out_3 = avg_dups(tors_type_list, v3_eq_2d)
        v3_eq_mean = np.reshape(out_3, (out_3.shape[0]))

        return tors_type_unique, v1_eq_mean, \
               v2_eq_mean, v3_eq_mean,  tors_length_mean, phase, hybrid_mean
    elif mdout == 'all':
        return tors_type_list, v1_eq, \
               v2_eq, v3_eq, tors_length_list, phase, hybrid_list
        


