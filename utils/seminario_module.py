#!/usr/bin/env python3

import numpy as np


def get_force_constant(i, j, diff_AB, r_AB, WR, VR,):

    unit_vector = diff_AB / r_AB
    k_AB = 0 

    for i in range(0,3):
        dot_product = abs(np.dot(unit_vector, VR[:,i]))
        k_AB = k_AB + ( WR[i] * dot_product )

    k_AB = - k_AB * 0.5 # Convert to OPLS form

    return k_AB

def get_ModSem_FcBonds(i, j, diff_AB, r_AB, hess):

    k_bond = 0
    sub_hess = hess[(i * 3):((i + 1)*3),(j * 3):((j + 1)*3)]
    WR, VR = np.linalg.eig(sub_hess)
    k_AB = get_force_constant(i, j, diff_AB, r_AB, WR, VR)
    sub_hess = hess[(j * 3):((j + 1)*3),(i * 3):((i + 1)*3)]
    WR, VR = np.linalg.eig(sub_hess)
    k_BA = get_force_constant(j, i, diff_AB, r_AB, WR, VR)
    kk = (k_AB + k_BA)/2.0

    return kk

def get_ModSem_FcAngles(i, j, diff_AB, r_AB, hess):

    k_bond = 0
    sub_hess = hess[(i * 3):((i + 1)*3),(j * 3):((j + 1)*3)]
    WR, VR = np.linalg.eig(sub_hess)
    k_AB = get_force_constant(i, j, diff_AB, r_AB, WR, VR)
    sub_hess = hess[(j * 3):((j + 1)*3),(i * 3):((i + 1)*3)]
    WR, VR = np.linalg.eig(sub_hess)
    k_BA = get_force_constant(j, i, diff_AB, r_AB, WR, VR)
    kk = (k_AB + k_BA)/2.0

    return kk
    
    


    
