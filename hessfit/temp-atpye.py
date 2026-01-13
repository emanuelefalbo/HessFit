import geom2atype as g2a
import numpy as np


elements = ["C"]*6 + ["H"]*6
coords = np.array([
    [ 1.4,  0.0, 0.0],
    [ 0.7,  1.2, 0.0],
    [-0.7,  1.2, 0.0],
    [-1.4,  0.0, 0.0],
    [-0.7, -1.2, 0.0],
    [ 0.7, -1.2, 0.0],
    [ 2.4,  0.0, 0.0],
    [ 1.2,  2.1, 0.0],
    [-1.2,  2.1, 0.0],
    [-2.4,  0.0, 0.0],
    [-1.2, -2.1, 0.0],
    [ 1.2, -2.1, 0.0]
])

print(g2a.assign_amber_atom_types(elements, coords))
print(g2a.assign_gaff_atom_types(elements, coords))