import openbabel.pybel as pybel
from openbabel import openbabel
from src.stringfile_helper_functions import build_bond_map
from src.cut_molecule import cut_molecule_main, make_cut_molecule, find_all_cuts, make_cut
from src.stringfile_to_rdkit import stringfile_to_rdkit
from rdkit.Chem import Draw
from rdkit.Chem.AllChem import Compute2DCoords

#from stringfile_helper_functions import build_bond_map
#from cut_molecule import cut_molecule_main, make_cut_molecule, find_all_cuts, make_cut
#from stringfile_to_rdkit import stringfile_to_rdkit

def get_educt(strfile):
    with open(strfile) as f:
        content = f.readlines()
    num_atoms = int(content[0])
    xyz_str_educt: str = "".join(content[:(num_atoms + 2)])

    educt = pybel.readstring("xyz", xyz_str_educt)
    return educt

def get_product(strfile):
    # read xyz data as string
    with open(strfile) as f:
        content = f.readlines()
    num_atoms = int(content[0])
    xyz_str_product: str = "".join(content[len(content) - (num_atoms + 2):])    # string representing product

    product = pybel.readstring("xyz", xyz_str_product)
    return product

def get_removed_atoms(cuts, molecule, lookup_dict):
    ban_list = []
    for c in cuts:
        for child in molecule[lookup_dict.get(c)].children:
            #ban_list.append(child)
            ban_list.append(child+1)
        ban_list += get_removed_atoms(molecule[lookup_dict.get(c)].children, molecule, lookup_dict)
    return ban_list

def check_product(original_strfile, modified_strfile, cuts, ordering, molecule, lookup_dict):
    if modified_strfile == "NO REACTION":
        return False
    # read the product of both files
    original_product = get_product(original_strfile)
    modified_product = get_product(modified_strfile)
    # get the bond maps of both products
    original_bmap = build_bond_map(original_product)
    modified_bmap = build_bond_map(modified_product)

    banned_atoms = get_removed_atoms(cuts, molecule, lookup_dict)

    #print("ORIGINAL!")
    original_bonds = set()
    for bond in original_bmap.keys():
        if bond[0] not in banned_atoms and bond[1] not in banned_atoms:
            original_bonds.add((ordering.get(bond[0]), ordering.get(bond[1]), original_bmap.get(bond)))
            #original_bonds.add((int(ordering.get(str(bond[0]))),int(ordering.get(str(bond[1])))))

    #print("MODIFIED!")
    modified_bonds = set()
    for bond in modified_bmap:
        #print(bond)
        modified_bonds.add((bond[0], bond[1], modified_bmap.get(bond)))
        #modified_bonds.add((bond[0], bond[1]))

    if original_bonds.difference(modified_bonds) == set() and modified_bonds.difference(original_bonds) == set():
        return True
    else:
        return False

def check_educt_to_product(stringfile):
    # get xyz data for eduxt and product
    educt = get_educt(stringfile)
    product = get_product(stringfile)

    # get bmap of both
    educt_bmap = build_bond_map(educt)
    product_bmap = build_bond_map(product)

    educt_bonds = set()
    for bond in educt_bmap:
        #modified_bonds.add((bond[0], bond[1], modified_bmap.get(bond)))
        educt_bonds.add((bond[0], bond[1]))

    product_bonds = set()
    for bond in product_bmap:
        #modified_bonds.add((bond[0], bond[1], modified_bmap.get(bond)))
        product_bonds.add((bond[0], bond[1]))

    if educt_bonds.difference(product_bonds) == set() and product_bonds.difference(educt_bonds) == set(): # no reaction happened
        return False
    else:
        return True

def stringfile_tester_main():
    #original_strfile = "xyz_test_files/reaction0001/stringfile.xyz0001"
    modified_strfile = "../xyz_test_files/reaction0182/7/stringfile.xyz0000"

    original_strfile = "../xyz_test_files/reaction0182/stringfile.xyz0182"

    cuts = {7}

    print("starting!")
    #stringfile_to_rdkit(modified_strfile, False)
    rdk_mol, atom_core, energy_curve = stringfile_to_rdkit(original_strfile, False)
    #Compute2DCoords(rdk_mol)
    #Draw.MolToFile(rdk_mol,'derp.png')

    #print(atom_core)
    molecule, lookup_dict = make_cut_molecule(rdk_mol, atom_core)
    xyz_file, ordering = make_cut(rdk_mol, cuts, molecule, lookup_dict)
    #print(xyz_file)

    print(check_product(original_strfile, modified_strfile, cuts, ordering, molecule, lookup_dict))

if __name__ == '__main__':
    stringfile_tester_main()
