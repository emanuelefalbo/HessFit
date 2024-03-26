#!/usr/bin/env python3

import numpy as np
np.seterr(divide='ignore', invalid='ignore')

def get_one_hot(values, indexes):
    one_hot = np.eye(np.max(indexes) + 1)[indexes]
    counts = np.sum(one_hot, axis=0)
    average = np.sum((one_hot.T * values), axis=1) / counts
    average = average[~np.isnan(average)]             # Remove NaN 
    return average

def reduce_bond_list(bond_list):
    splitted = list()
    [ splitted.append(i.split()) for i in bond_list]
    data = { (tuple(sorted(item))) for item in splitted}
    new = list(map(list, data))
    reduced = list(' '.join(item) for item in new)
    return reduced

def make_list_unique(var_list, k_values):
    splitted = [ i.split() for i in var_list]
    indexes = [i for i in range(len(splitted))]
    # Get indexes to be removed
    id2del = list()
    for i in range(len(splitted[:])-1):
        for j in range(i+1,len(splitted[:])):
            if (splitted[i] == splitted[j][::-1]) and (i !=j) :
                id2del.append(j)
                indexes[j] = indexes[i]
                break

    # pop indexes out; change res_bonds name 
    if len(id2del) != 0:
        j = id2del[0]
        res_bonds = splitted.copy()
        for n in range(len(id2del)):
            j = id2del[n] - n
            res_bonds.pop(j)
        reduced = list(' '.join(item) for item in res_bonds)
    
        indexes = np.array(indexes)
        k_values_ave = get_one_hot(k_values, indexes)
        return reduced, k_values_ave
    else:
        # print("No duplicates have been found")
        return var_list, k_values



