#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import parser_gau as pgau

def commandline_parser():
    parser = argparse.ArgumentParser(prog='2hessian.py', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-file',"--f", help='input file containig xyz coordinate')
    return parser

def calculate_wilson_b_matrix(cartesian_coordinates, internal_coordinates):
    num_atoms = len(cartesian_coordinates)
    num_internal_coords = len(internal_coordinates)
    b_matrix = np.zeros((num_internal_coords, 3 * num_atoms))
    
    # Compute the partial derivatives of internal coordinates with respect to Cartesian coordinates
    for i in range(num_internal_coords):
        for j in range(3 * num_atoms):
            # Here, you would compute the partial derivative based on the specific internal coordinate definition
            # Update the B matrix accordingly
            print("out<", i, j)
            b_matrix[i][j] = compute_partial_derivative(internal_coordinates[i], cartesian_coordinates, j)
    
    return b_matrix

def compute_partial_derivative(internal_coord, coords, j):
    if internal_coord == 'bond_length':
        atom1_index = 0  # index of the first atom in the bond
        atom2_index = 1  # index of the second atom in the bond
        print(atom1_index, atom2_index, j)
        if j == 3 * atom1_index:
            a = (coords[j] - coords[atom2_index])
            # dist = np.linalg.norm(coords[atom1_index:atom1_index + 3] - coords[atom2_index:atom2_index + 3])
            # dist = np.linalg.norm(coords[3* atom1_index:3*atom1_index + 3] - coords[3*atom2_index:3*atom2_index + 3])
            # return(a/dist) 
            return(1.0)
            print("")
            # return(a/dist)
            # return (coords[j] - coords[3 * atom2_index]) / np.linalg.norm(
            #     coords[3 * atom1_index:3 * atom1_index + 3] - coords[3 * atom2_index:3 * atom2_index + 3]
            # )
        elif j == 3 * atom2_index:
        #     a = (coords[j] - coords[3 * atom2_index])
        #     dist = np.linalg.norm(coords[3 * atom1_index:3 * atom1_index + 3] - coords[3 * atom2_index:3 * atom2_index + 3])
        #     print(a/dist) 
        #     print("")
            return 2.0
        #     print(coords[j])
        #     return (coords[j] - coords[3 * atom1_index]) / np.linalg.norm(
        #         coords[3 * atom1_index:3 * atom1_index + 3] - coords[3 * atom2_index:3 * atom2_index + 3]
        #     )
        else:
            return 0.0  # derivative is 0 for other Cartesian coordinates
    else:
        return 0.0  # derivative is 0 for other internal coordinates


def main():
    parser = commandline_parser()
    opts = parser.parse_args()
    text_qm_fchk = pgau.store_any_file(opts.f)
    N_atoms, qm_XYZ = pgau.read_XYZ(text_qm_fchk) 
    ric_list, force_1D = pgau.read_RicDim_Grad(text_qm_fchk)
    print(ric_list)
    print(qm_XYZ)
    print("")
    print(qm_XYZ[0][:])


    # Example internal coordinate definitions
    # Replace these with your specific internal coordinate definitions
    internal_coordinates = ['bond_length', 'bond_length', 'bond_length',
                            'bond_angle', 'bond_angle', 'dihedral_angle']
    # internal_coordinates = ['internal_coord_1', 'internal_coord_2', 'internal_coord_3']
    
    # Example Cartesian coordinates as a 1D array
    # Replace this with your specific Cartesian coordinates
    # Example: cartesian_coordinates = np.array([x1, y1, z1, x2, y2, z2, ...])
    # cartesian_coordinates = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])

    # Calculate the Wilson B matrix
    wilson_b_matrix = calculate_wilson_b_matrix(qm_XYZ, internal_coordinates)
    print("Wilson B matrix:")
    print(wilson_b_matrix)

if __name__ == "__main__":
    main()
