from mod import *
from rdkit.Chem import MolFromMolFile
from rdkit.Chem import rdDepictor
from rdkit.Chem import MolToMolBlock
from rdkit.Chem.AllChem import EmbedMolecule

g = smiles("CCO")

ids = []
chemical_elements = []
atom_ids = []
charges = []
degrees = []
for v in g.vertices:
    ids.append(v.id)
    chemical_elements.append(v.stringLabel)
    atom_ids.append(int(v.atomId))
    charges.append(int(v.charge))
    degrees.append(v.degree)

print(ids)
print(chemical_elements)
print(atom_ids)
print("charge!")
print(charges)
print(degrees)

source_ids = []
target_ids = []
string_labels = []

for e in g.edges:
    source_ids.append(e.source.id)
    target_ids.append(e.target.id)
    string_labels.append(e.stringLabel)

print(source_ids)
print(target_ids)
print(string_labels)
mol_title = "  mol file\n"
mol_program = "  ModToMol\n"
mol_comment = "\n"
mol_properties_block = ""
# number of atoms, number of bonds, number of atom list, Chiral flag, number of stext entries, number of lines of additional properties, mol version
mol_counts = str(len(ids)) + "  " + str(len(source_ids)) + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0999" + "  " + "V2000\n"

mol_atom_block = ""
for vertex in g.vertices:
    charge_value = "0"
    if vertex.charge != 0:
        charge_value = "5"
        mol_properties_block += "M  CHG  1  " + str(vertex.id + 1) + "  " + str(vertex.charge) + "\n"
    atom_string = "    0.0000    0.0000    0.0000 " + str(vertex.stringLabel).replace("-", "").replace("+", "") + "  " + "0" + "  " + charge_value + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0" "  " + "0" + "  " + "0" "  " + "0" + "  " + "0\n"
    mol_atom_block += atom_string

mol_bond_block = ""
for edge in g.edges:
    if edge.stringLabel == "-":
        bond_type = "1"
    elif edge.stringLabel == "=":
        bond_type = "2"
    elif edge.stringLabel == "#":
        bond_type = "3"
    elif edge.stringLabel == ":":
        bond_type = "4"
    else:
        bond_type = "5"

    bond_string = "  " + str(edge.source.id + 1) + "  " + str(edge.target.id + 1) + "  " + bond_type + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0\n"
    mol_bond_block += bond_string

mol_properties_block += "M  END\n"


f = open("molfile.mol", "w")
complete_mol_string = mol_title + mol_program + mol_comment + mol_counts + mol_atom_block + mol_bond_block + mol_properties_block
f.write(complete_mol_string)
f.close()

mol = MolFromMolFile('molfile.mol')       # read mol from file
m = rdDepictor.Compute2DCoords(mol)     # generate 2d coordinates
EmbedMolecule(mol, randomSeed=0xf00d)   # generate 3d coordinates

# save to file
f = open("molfile3D.mol", "a")
f.write(MolToMolBlock(mol))
f.close()
