#!/usr/bin/env python3

import sys

sys.path.insert(0, ".")

from ric_reference import (
    build_canonical_ric,
    validate_canonical_ric
)


log_file = "but_qm.log"
fchk_file = "but_qm.fchk"


ric = build_canonical_ric(
    log_file,
    fchk_file
)

validate_canonical_ric(ric)


print("\nCanonical RIC successfully constructed")
print("---------------------------------------")

print("Number of atoms       :", ric.natoms)
print("Number of RICs        :", ric.n_ric)
print("Number of bonds       :", ric.n_bonds)
print("Number of angles      :", ric.n_angles)
print("Number of dihedrals   :", ric.n_dihedrals)

print("\nRIC Hessian shape:")
print(ric.hessian_ric.shape)

print("\nFirst 5 bonds:")
for x in ric.bonds[:5]:
    print(x)

print("\nFirst 5 angles:")
for x in ric.angles[:5]:
    print(x)

print("\nFirst 5 dihedrals:")
for x in ric.dihedrals[:5]:
    print(x)

print("\nRIC Hessian symmetry check:")
print(
    abs(
        ric.hessian_ric -
        ric.hessian_ric.T
    ).max()
)