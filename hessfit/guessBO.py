import numpy as np

def guess_bond_order(i, j, elements, coords):
    """
    Guess bond order for mol2 export.
    Returned values are strings: '1', '2', '3'
    """
    ei = elements[i]
    ej = elements[j]
    d = np.linalg.norm(coords[i] - coords[j])

    # ---- C–C ----
    if ei == "C" and ej == "C":
        if d < 1.22:
            return "3"   # C≡C
        elif d < 1.34:
            return "2"   # C=C
        else:
            return "1"

    # ---- C–O / O–C ----
    if (ei, ej) in [("C", "O"), ("O", "C")]:
        if d < 1.30:
            return "2"   # C=O
        else:
            return "1"

    # ---- C–N / N–C ----
    if (ei, ej) in [("C", "N"), ("N", "C")]:
        if d < 1.30:
            return "2"
        else:
            return "1"

    # ---- Default ----
    return "1"
