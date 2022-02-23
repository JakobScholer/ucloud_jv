import openbabel.pybel as pybel
from src.stringfile_helper_functions import build_bond_map

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
        if len(molecule[lookup_dict.get(c)].id) > 1: # if the node is a big node
            for id in molecule[lookup_dict.get(c)].id:
                if id != c:
                    print(id)
                    print(c)
                    ban_list.append(id+1)
        for child in molecule[lookup_dict.get(c)].children: # for all childs add the to ban list
            #ban_list.append(child)
            ban_list.append(child+1)
        ban_list += get_removed_atoms(molecule[lookup_dict.get(c)].children, molecule, lookup_dict) # reapeat until leaf nodes are reached
    return ban_list

def get_removed_atoms(cuts, molecule, lookup_dict, rdk_mol):
    #rdk_mol = RWMol(rdk_mol) # typecast mol as mol object
    ban_list = set()
    keep_list = set()
    for c in cuts: # for each cut add them to ban list and their childs
        ban_list.add(c)
        for child in molecule[lookup_dict.get(c)].children:
            ban_list.add(child)
    for cut in cuts: # for each cut , find if it has neighbors which is not banned. Hence they are not banned
        for neighbor_atom in rdk_mol.GetAtomWithIdx(cut).GetNeighbors():
            if neighbor_atom.GetIdx() not in ban_list:
                keep_list.add(cut)
    return [x+1 for x in ban_list if x not in keep_list]


def check_product(original_strfile, modified_strfile, cuts, ordering, molecule, lookup_dict, rdk_mol):
    # read the product of both files
    original_product = get_product(original_strfile)
    modified_product = get_product(modified_strfile)
    # get the bond maps of both products
    original_bmap = build_bond_map(original_product)
    modified_bmap = build_bond_map(modified_product)

    banned_atoms = get_removed_atoms(cuts, molecule, lookup_dict, rdk_mol)

    original_bonds = set()
    for bond in original_bmap.keys():
        if bond[0] not in banned_atoms and bond[1] not in banned_atoms:
            original_bonds.add((ordering.get(int(bond[0]),int(bond[0])), ordering.get(int(bond[1]),int(bond[1])), original_bmap.get(bond)))

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
