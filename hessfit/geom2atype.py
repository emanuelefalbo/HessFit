import numpy as np

# Covalent radii (Å)

COV_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71,
    "O": 0.66, "S": 1.05, "P": 1.07,
    "F": 0.57, "Cl": 1.02, "Br": 1.20, "I": 1.39
}


def build_bonds(elements, coords, scale=1.2):
    bonds = [[] for _ in elements]
    n = len(elements)

    for i in range(n):
        for j in range(i+1, n):
            ri = COV_RADII.get(elements[i], 0.7)
            rj = COV_RADII.get(elements[j], 0.7)
            cutoff = scale * (ri + rj)
            dist = np.linalg.norm(coords[i] - coords[j])
            if dist <= cutoff:
                bonds[i].append(j)
                bonds[j].append(i)
    return bonds

def infer_hybridization(atom, bonds, elements):
    nb = len(bonds[atom])

    if elements[atom] == "C":
        if nb == 4:
            return "sp3"
        elif nb == 3:
            return "sp2"
        elif nb == 2:
            return "sp"
    if elements[atom] == "N":
        if nb == 3:
            return "sp3"
        elif nb == 2:
            return "sp2"
    return "unknown"

def is_aromatic(atom, bonds, elements):
    if elements[atom] not in ("C", "N"):
        return False
    if len(bonds[atom]) != 3:
        return False
    return True  # crude but effective for benzene-like rings

def assign_amber_type(i, elements, bonds):
    elem = elements[i]
    nb = len(bonds[i])
    neigh_elems = [elements[j] for j in bonds[i]]

    if elem == "C":
        if is_aromatic(i, bonds, elements):
            return "CA"
        if nb == 4:
            return "CT"
        if nb == 3:
            if "O" in neigh_elems:
                return "C"   # carbonyl
            return "CG"     # alkene
        if nb == 2:
            return "C1"

    if elem == "H":
        parent = bonds[i][0]
        if elements[parent] == "C":
            if is_aromatic(parent, bonds, elements):
                return "HA"
            else:
                return "HC"
        return "H"

    if elem == "O":
        if nb == 1:
            parent = bonds[i][0]
            if elements[parent] == "C":
                return "O"
        if nb == 2:
            return "OS"

    if elem == "N":
        if is_aromatic(i, bonds, elements):
            return "NA"
        if nb == 3:
            return "N3"
        if nb == 2:
            return "N"

    return elem  # fallback


def assign_gaff_type(i, elements, bonds, aromatic_atoms):
    elem = elements[i]
    nb = len(bonds[i])
    neigh = bonds[i]
    neigh_elems = [elements[j] for j in neigh]
    hyb = infer_hybridization(i, elements, bonds)

    # ----- CARBON -----
    if elem == "C":
        if i in aromatic_atoms:
            return "ca"
        if hyb == "sp3":
            return "c3"
        if hyb == "sp2":
            if "O" in neigh_elems:
                return "c"
            return "c2"
        if hyb == "sp":
            return "c1"

    # ----- HYDROGEN -----
    if elem == "H":
        parent = neigh[0]
        if parent in aromatic_atoms:
            return "ha"
        parent_hyb = infer_hybridization(parent, elements, bonds)
        if parent_hyb == "sp3":
            return "h3"
        if parent_hyb == "sp2":
            return "h2"
        return "h1"

    # ----- OXYGEN -----
    if elem == "O":
        if nb == 1:
            return "o"
        if nb == 2:
            return "os"

    # ----- NITROGEN -----
    if elem == "N":
        if i in aromatic_atoms:
            return "na"
        if hyb == "sp3":
            return "n3"
        if hyb == "sp2":
            if "C" in neigh_elems:
                return "n"
            return "n2"

    # ----- FALLBACK -----
    return elem.lower()

def find_rings(bonds, max_len=7):
    rings = set()
    n = len(bonds)

    def dfs(path, start):
        current = path[-1]
        for nb in bonds[current]:
            if nb == start and len(path) > 2:
                rings.add(tuple(sorted(path)))
            elif nb not in path and len(path) < max_len:
                dfs(path + [nb], start)

    for i in range(n):
        dfs([i], i)

    return [list(r) for r in rings]


def detect_aromatic_atoms(elements, bonds):
    rings = find_rings(bonds)
    aromatic = set()

    for ring in rings:
        if len(ring) not in (5, 6):
            continue

        ok = True
        for i in ring:
            hyb = infer_hybridization(i, elements, bonds)
            if hyb != "sp2":
                ok = False
                break

        if ok:
            aromatic.update(ring)

    return aromatic


def assign_amber_atom_types(elements, coords):
    bonds = build_bonds(elements, coords)
    types = []
    for i in range(len(elements)):
        t = assign_amber_type(i, elements, bonds)
        types.append(t)
    return types

def assign_gaff_atom_types(elements, coords):
    bonds = build_bonds(elements, coords)
    aromatic_atoms = detect_aromatic_atoms(elements, bonds)

    types = []
    for i in range(len(elements)):
        t = assign_gaff_type(i, elements, bonds, aromatic_atoms)
        types.append(t)

    return types
