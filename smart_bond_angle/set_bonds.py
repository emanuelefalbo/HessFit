#!/usr/bin/env python3

import numpy as np
from parser_gau import flat_list

def set_bonds(coords, ele_list, type_list,
              bond_list, k_bonds):
    """
    Order Bond Force Constant
    """

    # print(len(type_list))
    bond_lenght_list = []
    bond_type_list = []
    for k in bond_list:
        i = k[0] -1 
        j = k[1] -1
        diff_AB = coords[i,:] - coords[j,:]
        r_AB = np.linalg.norm(diff_AB)
        bond_lenght_list.append(r_AB)
        bond_type_list.append(type_list[i] + ' ' + type_list[j]) 
    
    bond_type_list = flat_list(bond_type_list)
    print(bond_type_list)
        
    for k in range(len(bond_list)):
        msg = (
              f'HrmStr1 {bond_type_list[k]}  {k_bonds[k]:.3f} ' 
              f' {bond_lenght_list[k]:.3f} '
              )
        print(msg)


    