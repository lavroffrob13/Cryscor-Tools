from ase import Atoms
from ase.io import read, write
import numpy as np
import sys
def displace_along_bond(atoms, atom_index1, atom_index2, displacement):
    pos1 = atoms.positions[atom_index1]
    pos2 = atoms.positions[atom_index2]

    #calculate the bond vector and normalize it
    bond_vector = pos2 - pos1
    bond_length = np.linalg.norm(bond_vector)
    bond_direction = bond_vector / bond_length

    atoms.positions[atom_index1] += displacement * bond_direction

def main():
    file_path = sys.argv[1]
    
    atoms = read(file_path)
    
    atom_index1 = 1  #index of the atom to displace
    atom_index2 = 2  #index of the atom defining the bond direction
    displacement = -0.4  #displacement distance along the bond

    displace_along_bond(atoms, atom_index1, atom_index2, displacement)
    
    write('displaced_structure.fort.34', atoms) #crystal
    write('displaced_structure.vasp', atoms, direct=True) #vasp POSCAR, easier to read by a human
    print("Displaced structure saved as 'displaced_structure.fort.34 and .vasp'")

if __name__ == "__main__":
    main()
