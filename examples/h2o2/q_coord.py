import numpy as np
import sys


class Myclass:
     
    def compute_internal_coordinates(self, xyz):
        """
        Compute internal coordinates from Cartesian coordinates of a molecule.
        Parameters:
        - xyz: numpy array, the Cartesian coordinates of the atoms in the molecule
        Returns:
        - internal_coords: numpy array, the internal coordinates of the molecule
        """
        num_atoms = xyz.shape[0]
        
        # Calculate bond lengths
        bond_lengths = []
        for i in range(num_atoms):
            for j in range(i+1, num_atoms):
                bond_lengths.append(np.linalg.norm(xyz[i] - xyz[j]))
                
        print(bond_lengths)
        print("")
        
        # Calculate bond angles
        bond_angles = []
        for i in range(num_atoms):
            for j in range(i+1, num_atoms):
                for k in range(j+1, num_atoms):
                    v1 = xyz[i] - xyz[j]
                    v2 = xyz[k] - xyz[j]
                    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                    bond_angles.append(np.arccos(cos_angle))
        
        print(bond_angles)
        print("")
        
        # Calculate dihedral angles
        dihedral_angles = []
        for i in range(num_atoms):
            for j in range(num_atoms):
                if i != j:
                    for k in range(num_atoms):
                        if k != i and k != j:
                            for m in range(num_atoms):
                                if m != i and m != j and m != k:
                                    v1 = xyz[i] - xyz[j]
                                    v2 = xyz[k] - xyz[j]
                                    v3 = xyz[m] - xyz[k]
                                    normal1 = np.cross(v1, v2)
                                    normal2 = np.cross(v2, v3)
                                    cos_angle = np.dot(normal1, normal2) / \
                                                (np.linalg.norm(normal1) * np.linalg.norm(normal2))
                                    dihedral_angles.append(np.arccos(cos_angle))
                                    
        print(dihedral_angles)
        print("")
        
        internal_coords = np.concatenate((bond_lengths, bond_angles, dihedral_angles))
        return internal_coords
    
    def compute_wilson_matrix(self, internal_coords):
        """
        Compute the Wilson matrix from internal coordinates.

        Parameters:
        - internal_coords: numpy array, the internal coordinates of the molecule

        Returns:
        - wilson_matrix: numpy array, the Wilson matrix
        """
        num_internal_coords = internal_coords.shape[0]
        num_atoms = 3  # Assuming 3-dimensional Cartesian coordinates
        
        wilson_matrix = np.zeros((num_internal_coords, num_atoms * num_atoms))
        index = 0
        
        # Populate the Wilson matrix with partial derivatives of internal coordinates with respect to Cartesian coordinates
        for i in range(num_internal_coords):
            perturbation = np.zeros_like(internal_coords)
            perturbation[i] = 0.001  # Small perturbation
            xyz_plus = self.cartesian_from_internal(internal_coords + perturbation)
            xyz_minus = self.cartesian_from_internal(internal_coords - perturbation)
            partial_derivative = (xyz_plus - xyz_minus) / (2 * 0.001)
            wilson_matrix[i] = partial_derivative.flatten()
        
        return wilson_matrix
    
    def compute_internal_gradient(self, xyz, grad_cartesian):
        """
        Compute the internal coordinate gradient from Cartesian coordinate gradient.

        Parameters:
        - xyz: numpy array, the Cartesian coordinates of the atoms in the molecule
        - gradq_cartesian: numpy array, the Cartesian coordinate gradient

        Returns:
        - gradq_internal: numpy array, the internal coordinate gradient
        """
        # Compute internal coordinates
        internal_coords = self.compute_internal_coordinates(xyz)
        
        # Construct the rectangular Wilson matrix using internal coordinates
        wilson_matrix = self.compute_wilson_matrix(internal_coords)
        
        # Compute the internal coordinate gradient by multiplying with the Wilson matrix
        gradq_internal = np.dot(wilson_matrix.T, grad_cartesian)
        
        print(gradq_internal)
        
        return gradq_internal
    
    
    
def read_xyz_file(file_path):
    coord = []
    with open(file_path, 'r') as file:
        num_atoms = int(file.readline())
        file.readline()  # skip comment line
        for _ in range(num_atoms):
            atom_info = file.readline().split()
            atom_symbol = atom_info[0]
            x, y, z = map(float, atom_info[1:4])
            coord.append((x, y, z))
    
    coord = np.asarray(coord)
    return coord

def store_any_file(fname):
    """ 
    Reading Any file 
    """
    with open(fname, 'r') as f:
        all_lines = []
        for line in f:
            all_lines.append(line.strip())

    return all_lines

def range2(start, end):
     return range(start, end+1)

def read_gradient_file(all_lines):
    gradient=[]
    with open(file_path, 'r') as file:
        for s in range(len(all_lines)): 
            if "Cartesian Gradient" in all_lines[s]:# Get no of Atoms
                N_force = int(all_lines[s].split()[4])
                nlines = int(N_force/5.0) + 1
                start = s
                for e in range2(start + 1, start + nlines):
                    gradient.append(all_lines[e].split() )

    return gradient

# Example usage:
file_path = sys.argv[1] # Provide the path to your XYZ file
file_2 = sys.argv[2] # Provide the path to your XYZ file
xyz = read_xyz_file(file_path)

text = store_any_file(file_2)
xyz_grad = read_gradient_file(text)
print(xyz)
 
obj = Myclass()
obj.compute_internal_gradient(xyz, xyz_grad)