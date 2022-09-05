#!/usr/bin/env python3

import numpy as np


def unit_vector_N(u_BC, u_AB):
    # Calculates unit normal vector which is perpendicular to plane ABC
    cross_product = np.cross(u_BC, u_AB)
    norm_u_N = np.linalg.norm(cross_product)
    u_N = cross_product / norm_u_N
    return u_N

def angle_force_constant(r_AB, r_BC, u_AB, u_BC, \
                         WR_AB, VR_AB, WR_BC, VR_BC):

    u_N = unit_vector_N(u_BC, u_AB)
    u_PA = np.cross(u_N, u_AB)
    u_PA = u_PA/np.linalg.norm(u_PA)
    u_PC = np.cross(u_N, u_BC)
    u_PC = u_PC / np.linalg.norm(u_PC)

    k_PA = 0.0 
    k_PC = 0.0 
    for i in range(0,3):
        # dot_product = abs(dot_product(u_PA, VR_AB[:,i]))
        k_PA = k_PA +  WR_AB[i] * abs(np.dot(u_PA, VR_AB[:,i]))    
        k_PC = k_PC +  WR_BC[i] * abs(np.dot(u_PC, VR_BC[:,i]))    

    k_theta = ( 1 / ( (r_AB**2) * k_PA) ) + ( 1 / ( (r_BC**2) * k_PC) ) 
    k_theta = 1/k_theta
    k_theta = - k_theta #Change to OPLS form
    k_theta = abs(k_theta * 0.5) #Change to OPLS form

    return k_theta

def bond_force_constant(diff_AB, r_AB, WR, VR,):

    unit_vector = diff_AB / r_AB
    k_AB = 0 
    for i in range(0,3):
        k_AB = k_AB +  WR[i]* abs(np.dot(unit_vector, VR[:,i]))
        # k_AB = k_AB + ( WR[i] * dot_product )

    k_AB = - k_AB * 0.5 # Convert to OPLS form

    return k_AB



def get_ModSem_FcBonds(i, j, diff_AB, r_AB, hess):
    sub_hess = hess[(i * 3):((i + 1)*3),(j * 3):((j + 1)*3)]
    WR, VR = np.linalg.eig(sub_hess)
    k_AB = bond_force_constant(diff_AB, r_AB, WR, VR)
    sub_hess = hess[(j * 3):((j + 1)*3),(i * 3):((i + 1)*3)]
    WR, VR = np.linalg.eig(sub_hess)
    k_BA = bond_force_constant(diff_AB, r_AB, WR, VR)
    kk = (k_AB + k_BA)/2.0

    return kk

def get_ModSem_FcAngles(i, j, k, r_AB, r_BC, u_AB, u_BC, hess):
    sub_hess_AB = hess[(i * 3):((i + 1)*3),(j * 3):((j + 1)*3)]
    WR_AB, VR_AB = np.linalg.eig(sub_hess_AB)
    sub_hess_BC = hess[(j * 3):((j + 1)*3),(k * 3):((k + 1)*3)]
    WR_BC, VR_BC = np.linalg.eig(sub_hess_BC)
    k_angle = angle_force_constant(r_AB, r_BC, u_AB, u_BC, \
                                   WR_AB, VR_AB, WR_BC, VR_BC)

    
    


    
