#!/usr/bin/env python3

import numpy as np

def set_bonds(coords, ele_list, type_list,
              bond_list, k_bonds):
    """
    Order Bond Force Constant
    """

    bond_lenght_list = []
    for k in bond_list:
        i = k[0] -1 
        j = k[1] -1
        diff_AB = coords[i,:] - coords[j,:]
        r_AB = np.linalg.norm(diff_AB)
        bond_lenght_list.append(r_AB)

    print(bond_lenght_list)

    