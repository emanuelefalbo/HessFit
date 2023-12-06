
# import openbabel
from openbabel import openbabel
import sys

# def read_xyz(file_path):
#     obConversion = openbabel.OBConversion()
#     obConversion.SetInFormat("xyz")

#     mol = openbabel.OBMol()
#     obConversion.ReadFile(mol, file_path)
#     print(mol)
    
#     force_field = openbabel.OBForceField.FindForceField("GAFF")
#     if force_field.Setup(mol):
#         print("Force field set to GAFF")
#     else:
#         print("Failed to set force field to GAFF")
#         return
    
#     # Get GAFF atom types
#     for atom in openbabel.OBMolAtomIter(mol):
#         force_field.GetAtomTypes(atom)
#         atom_types = atom.GetData()
#         print(atom_types)
        # if atom_types is not None:
        #     print(f"Atom {atom.GetIdx()}: GAFF Type - {atom_types}")
        # print(atom.GetIdx())
        # gaff_type = force_field.GetAtomType(atom)
        # print(f"Atom {atom.GetIdx()}: GAFF Type - {gaff_type}")

    # # Get atomic types
    # for atom in openbabel.OBMolAtomIter(mol):
    #     atomic_num = atom.GetAtomicNum()
    #     atom_type = openbabel.OBElementTable.GetSymbol(atomic_num)
    #     print(f"Atom {atom.GetIdx()}: {atom_type}")

from rdkit import Chem
from rdkit.Chem import AllChem

# Load molecule from XYZ file
def read_xyz(file):
    with open(file, 'r') as f:
        content = f.readlines()
        num_atoms = int(content[0])
        mol_block = ''.join(content[2:])  # Skipping the first two lines in XYZ format
        print(mol_block)
    
    # # Create an RDKit molecule object
    mol = Chem.MolFromMolBlock(mol_block)
    
    # # Generate 3D coordinates if not present
    # AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    
    # Assign MMFF94 atom types
    # AllChem.MMFFSanitizeMolecule(mol)
    # AllChem.MMFFAddHs(mol)
    # AllChem.MMFFOptimizeMolecule(mol)
    
    # Print atom types for each atom
    # for atom in mol.GetAtoms():
    #     atom_index = atom.GetIdx()
    #     atom_type = atom.GetProp('_TriposAtomType')
    #     print(f"Atom {atom_index + 1}: Type - {atom_type}")

if __name__ == "__main__":
    file_path = sys.argv[1]
    read_xyz(file_path)
