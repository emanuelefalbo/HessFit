
def gaff_vdw_parameters():
    """
    Returns GAFF/AMBER van der Waals parameters:
    Radius = Rmin/2 (Å)
    Well-depth = epsilon (kcal/mol)
    """

    return {
        # ---- CARBON ----
        "c":  (1.9080, 0.0860),
        "c2": (1.9080, 0.0860),
        "c3": (1.9080, 0.1094),
        "ca": (1.9080, 0.0860),
        "c1": (1.9080, 0.0860),

        # ---- HYDROGEN ----
        "h1": (1.3870, 0.0157),
        "h2": (1.3870, 0.0157),
        "h3": (1.4870, 0.0157),
        "ha": (1.4590, 0.0150),

        # ---- OXYGEN ----
        "o":  (1.6612, 0.2100),
        "os": (1.6837, 0.1700),
        "oh": (1.7210, 0.2104),

        # ---- NITROGEN ----
        "n":  (1.8240, 0.1700),
        "n2": (1.8240, 0.1700),
        "n3": (1.8240, 0.1700),
        "na": (1.8240, 0.1700),

        # ---- SULFUR ----
        "s":  (2.0000, 0.2500),  # generic divalent sulfur
        "sh": (2.0000, 0.2500),  # thiol sulfur
        "ss": (2.0000, 0.2500),  # disulfide sulfur
        "s2": (2.0000, 0.2500),  # sp2 sulfur
        "s4": (2.0000, 0.2500),  # hypervalent sulfur (sulfone/sulfoxide)
        "sx": (2.0000, 0.2500),  # aromatic sulfur
        "sy": (2.0000, 0.2500),  # aromatic sulfur variant
    }

def compute_vdw_table(atom_types):
    """
    Returns unique VDW parameters for each atom type used
    """
    vdw_db = gaff_vdw_parameters()
    vdw_table = {}

    for at in atom_types:
        if at not in vdw_db:
            raise KeyError(f"Missing VDW parameters for atom type '{at}'")

        vdw_table[at] = vdw_db[at]

    return vdw_table

def write_gaussian_vdw(vdw_table):
    lines = []
    for at, (radius, epsilon) in vdw_table.items():
        lines.append(["VDW", at, f"{radius:.4f}", f"{epsilon:.4f}"])
    return lines

def write_nbdir(
    v_type=3,  # AMBER arithmetic
    c_type=1,
    v_cutoff=12.0,
    c_cutoff=12.0
):
    return (
        "NBDir\n"
        f"{v_type:d} {c_type:d} {v_cutoff:6.2f} {c_cutoff:6.2f}"
    )

def write_nbterm(
    v_type=3,
    c_type=1,
    v_cutoff=12.0,
    c_cutoff=12.0,
    v_scale=1.0,
    c_scale=1.0
):
    return (
        "NBTerm\n"
        f"{v_type:d} {c_type:d} {v_cutoff:6.2f} {c_cutoff:6.2f} "
        f"{v_scale:6.3f} {c_scale:6.3f}"
    )

def write_nonbon():
    return (
        "NonBon\n"
        "3 1 12.0 12.0 "
        "0.0 0.0 0.5 "
        "0.0 0.0 0.8333"
    )

# #example usage
# atype_list = ['c', 'h1', 'o', 'n3', 'ha', 'ca']
# vdw_table = compute_vdw_table(atype_list)

# print(write_gaussian_vdw(vdw_table))
# print(write_nbdir())
# print(write_nbterm())
# print(write_nonbon())
