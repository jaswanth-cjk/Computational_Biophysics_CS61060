import math  # math module for mathematical functions
from Bio import PDB  # PDB module from Biopython for working with PDB files
# specific PDB-related functions from Biopython
from Bio.PDB import PDBParser, PPBuilder, calc_dihedral
import numpy as np  # numpy module for numerical operations
import warnings  # warnings module for handling warnings
from Bio import BiopythonWarning  # Biopython warning related to modules
warnings.simplefilter('ignore', BiopythonWarning)  # ignoring Biopython warning


# Open PDB file
with open("insulin.pdb", "r") as f:
    lines = f.readlines()


# Load PDB file
'''
Here we are loading the PDB file of any protein from the PDB source and save it as protein.pdb in the same directory
For example we are using the heamoglobin(1A3N) file here
'''
parser = PDB.PDBParser()
structure = parser.get_structure("Hemoglobin", "insulin.pdb")


# Extract backbone atoms
'''
This code snippet creates an empty list backbone_atoms, then loops through each line of the lines list, 
which was obtained by reading a PDB file. If the line starts with "ATOM" and the atom name (obtained from positions 12-15 of the line, stripped of whitespace)
is one of "N", "CA", "C", or "O", then the residue ID is extracted (obtained from positions 22-25 of the line and converted to an integer), 
 as well as the x, y, and z coordinates of the atom (obtained from positions 30-37, 38-45, and 46-53 of the line, respectively, and converted to floats). 
 Finally, a tuple of the residue ID, atom name, and coordinates is appended to the backbone_atoms list. 
 This code effectively extracts the coordinates of the backbone atoms (N, CA, C, and O) from the PDB file and stores them in a list.
'''
backbone_atoms = []
for line in lines:
    if line.startswith("ATOM") and line[12:16].strip() in ["N", "CA", "C", "O"]:
        residue_id = int(line[22:26])
        atom_name = line[12:16].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        backbone_atoms.append((residue_id, atom_name, x, y, z))


# Print coordinates of backbone atoms
'''
This code snippet is taking the backbone_atoms list created in the previous code block and using it to group atoms by residue ID and print them out in a formatted way.

The set(atom[0] for atom in backbone_atoms) creates a set of unique residue IDs by extracting the first element (residue ID) from each tuple in backbone_atoms.

For each residue ID, the code creates a new list residue_atoms that contains all the atoms in backbone_atoms with that residue ID. It then extracts the residue name from the first atom in residue_atoms and prints it along with the residue ID.

Finally, the code loops through each atom in residue_atoms and prints its name and coordinates.
'''
for residue_id in set(atom[0] for atom in backbone_atoms):
    residue_atoms = [atom for atom in backbone_atoms if atom[0] == residue_id]
    residue_name = residue_atoms[0][1][0:3]
    print(f"{residue_name} {residue_id}:")
    for atom in residue_atoms:
        print(f"{atom[1]}: {atom[2]:.3f} {atom[3]:.3f} {atom[4]:.3f}")


# Count number of residues in each chain and in total
'''
This code snippet iterates over the models, chains, and residues in a protein PDB file to count the number of amino acid residues in each chain and the total number of residues in the protein.

The variable total_res_count is initialized to 0. This variable will be used to keep track of the total number of amino acid residues in the protein.

The code iterates over the models in the structure using a for loop. If there is only one model in the structure, this loop will run only once.

Within each model, the code iterates over the chains using another for loop.

For each chain, the code initializes a variable res_count to 0. This variable will be used to count the number of amino acid residues in the chain.

The code then iterates over the residues in the chain using another for loop.

For each residue, the code checks whether it is a standard amino acid residue using the PDB.is_aa() function. If the residue is a standard amino acid, the res_count variable is incremented.

After counting the number of amino acid residues in the chain, the res_count value is added to the total_res_count variable.

The code then prints out the chain ID and the number of amino acid residues in the chain.

Finally, after iterating over all the chains in all the models, the code prints out the total number of amino acid residues in the protein.
'''
total_res_count = 0
for model in structure:
    for chain in model:
        res_count = 0
        for residue in chain:
            if PDB.is_aa(residue.get_resname(), standard=True):
                res_count += 1
        total_res_count += res_count
        print(f"Chain {chain.id}: {res_count} residues")
print(
    f"-----------Total number of residues in the Protein PDB File are: {total_res_count}-----------------")

'''
Here we are assigning the secondary structures(alpha helicx,beta strand,beta turn, coil) 
to the residues in PDB file following the Stride algorithm after calculating the dihedral angles of each residue using their backbone atoms.
Based on the Ramachandran plot the values are also matching 
'''


def assign_secondary_structure(phi, psi):

    # check if residue is in alpha helix region

    if (((-91 <= phi <= -40) & (-50 < psi < -5)) | ((20 < phi < 60) & (0 < psi < 95))):
        return "alpha helix"

    # check if residue is in beta strand region
    # elif (phi, psi) in beta_region:
    #     return "beta strand"
    elif ((-145 <= phi <= -119) & (117 < psi < 150)):
        return "beta strand"

    elif (((-100 <= phi <= -30) & (0 < psi < 30)) | ((120 < psi < 160) & (-70 < phi < 30))):

        return "beta turn"

    # otherwise, assign to coil region
    else:
        return "coil"


# Set up the PDB parser and parse the structure
parser = PDBParser()
structure = parser.get_structure('protein', 'insulin.pdb')

'''
This code snippet is used for calculating the phi and psi angles of every residue in every polypeptide chain of a given protein structure, 
using the PPBuilder and calc_dihedral functions from Biopython. 
It then assigns a secondary structure to each residue based on its phi and psi angles, 
using the assign_secondary_structure function .

The code then iterates over each polypeptide chain, and for each residue in the chain 
(excluding the first and last residues), it retrieves the atoms needed to calculate the phi and psi angles 
(C, N, CA, C' and N' atoms). It then uses the calc_dihedral function to calculate the dihedral angles, 
and converts them to degrees using the math.degrees function.

After calculating the phi and psi angles for a residue, the code assigns a secondary structure to that residue using the assign_secondary_structure function. 
The results (residue name, residue number, phi angle, psi angle, and secondary structure) are then printed to the console.
'''

# Get polypeptide chains from the structure
ppb = PPBuilder()
polypeptides = ppb.build_peptides(structure)

# Iterate over every residue in every polypeptide chain
for polypeptide in polypeptides:
    for i in range(-1, len(polypeptide)-1):
        # Get the atoms needed for the dihedral angle calculation
        C_i = polypeptide[i-1]['C']
        N_i = polypeptide[i]['N']
        CA_i = polypeptide[i]['CA']
        C_i1 = polypeptide[i]['C']
        N_i1 = polypeptide[i+1]['N']

        # Calculate the phi and psi angles
        phi = math.degrees(calc_dihedral(
            C_i.get_vector(), N_i.get_vector(), CA_i.get_vector(), C_i1.get_vector()))
        psi = math.degrees(calc_dihedral(
            N_i.get_vector(), CA_i.get_vector(), C_i1.get_vector(), N_i1.get_vector()))

        struct = assign_secondary_structure(phi, psi)
        # Print out the results
        print(
            f"Residue {polypeptide[i].get_resname()} {polypeptide[i].get_full_id()[3][1]}")
        print(f"Phi: {phi:.2f} degrees")
        print(f"Psi: {psi:.2f} degrees")
        print(struct)
