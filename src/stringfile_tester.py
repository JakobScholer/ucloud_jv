import openbabel.pybel as pybel
from openbabel import openbabel
from src.stringfile_to_rdkit import build_bond_map

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
    return ban_list

def check_product(original_strfile, modified_strfile, cuts, ordering, molecule, lookup_dict):
    if modified_strfile == "NO REACTION":
        return False
    print(modified_strfile)
    # read the product of both files
    original_product = get_product(original_strfile)
    modified_product = get_product(modified_strfile)
    # get the bond maps of both products
    original_bmap = build_bond_map(original_product)
    modified_bmap = build_bond_map(modified_product)

    banned_atoms = get_removed_atoms(cuts, molecule, lookup_dict)
    print(banned_atoms)
    print(" ")
    print(cuts)
    print(" ")
    print(original_bmap)
    print(" ")
    print(ordering)
    print(" ")
    original_bonds = set()
    for bond in original_bmap.keys():
        print(bond)
        if bond[0] not in banned_atoms and bond[1] not in banned_atoms:
            #original_bonds.add((int(ordering.get(str(bond[0]))),int(ordering.get(str(bond[1]))),original_bmap.get(bond)))
            original_bonds.add((int(ordering.get(str(bond[0]))),int(ordering.get(str(bond[1])))))
    print("### original_bonds ###")
    print(original_bonds)

    modified_bonds = set()
    for bond in modified_bmap:
        print(bond)
        #modified_bonds.add((bond[0], bond[1], modified_bmap.get(bond)))
        modified_bonds.add((bond[0], bond[1]))
    print("### modified_bonds ###")
    print(modified_bonds)

    if original_bonds.difference(modified_bonds) == set():
        print("true")
        return True
    else:
        print("False")
        return False

def stringfile_tester_main():
    print("Hello")

if __name__ == '__main__':
    stringfile = "../xyz_test_files/GCD_test_files/stringfile.xyz0009"
    #original_product = get_product(stringfile)
    print((1,2,3))
