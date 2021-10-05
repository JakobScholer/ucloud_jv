from mod import *
from rdkit.Chem import MolFromMolFile
from rdkit.Chem import rdDepictor
from rdkit.Chem import MolFromMolBlock
from rdkit.Chem.AllChem import EmbedMolecule
from rdkit.Chem.rdmolfiles import MolToXYZFile


def mod_to_xyz(g, toFile=True):
    # Count vertices
    vertex_counter = 0
    for v in g.vertices:
        vertex_counter += 1
    # count edges
    edge_counter = 0
    for e in g.edges:
        edge_counter += 1

    mol_string = ""
    mol_properties_block = ""
    # Title line (can be blank but line must exist)
    mol_string += "  mod->mol->xyz\n"
    # Program / file timestamp line
    # (Name of source program and a file timestamp)
    mol_string += "  ModToMol\n"
    # Comment line (can be blank but line must exist)
    mol_string += "\n"
    # Counts line
    # number of atoms, number of bonds, number of atom list, Chiral flag, number of stext entries, properties, mol version
    mol_string += str(vertex_counter) + "  " + str(edge_counter) + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0999" + "  " + "V2000\n"

    # Atom block
    # (1 line for each atom): x, y, z (in angstroms), element, etc.
    for vertex in g.vertices:
        charge_value = "0"
        if vertex.charge != 0:
            charge_value = "5"
            mol_properties_block += "M  CHG  1  " + str(vertex.id + 1) + "  " + str(vertex.charge) + "\n"
        atom_string = "    0.0000    0.0000    0.0000 " + str(vertex.stringLabel).replace("-", "").replace("+", "") + "  " + "0" + "  " + charge_value + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0" "  " + "0" + "  " + "0" "  " + "0" + "  " + "0\n"
        mol_string += atom_string

    # Bond block
    # (1 line for each bond): 1st atom, 2nd atom, type, etc.
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
        # each bond position is only allowed 3 spaces meaning the number "10" only has 1 space to the left for padding while "1" has two spaces of padding
        bond_string = str(edge.source.id + 1).rjust(3) + str(edge.target.id + 1).rjust(3) + "  " + bond_type + "  " + "0" + "  " + "0" + "  " + "0" + "  " + "0\n"
        mol_string += bond_string

    # Properties block
    mol_properties_block += "M  END\n"
    mol_string += mol_properties_block

    mol = MolFromMolBlock(mol_string, removeHs=False)       # Convert to rdkit mol format
    rdDepictor.Compute2DCoords(mol)     # generate 2d coordinates
    EmbedMolecule(mol, randomSeed=0xf00d)   # generate 3d coordinates

    # save to file
    if toFile:
        MolToXYZFile(mol, "mod_coordinates.xyz")
    else:
        return mol
