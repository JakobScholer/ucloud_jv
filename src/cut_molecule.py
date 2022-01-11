from rdkit.Chem.rdchem import AtomValenceException
from rdkit.Chem import rdmolops, GetSymmSSSR, Atom
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdDepictor import Compute2DCoords

class MoleculeNode:
    def __init__(self, molecule_id, node_type):
        self.id = molecule_id  #set of atom ID's
        self.children = set()  # A list of ints representing the list placement of the children
        self.root = node_type  # 1 = root, 0 = not root

# core = set(atom ID's) a set of ints
# rdk_mol is the molecule object from RDKit
def make_cut_molecule(rdk_mol, core):
    # for handling of bignodes
    def big_node_update(b_nodes, look_up, atom_list, data):
        intersection = data.intersection(atom_list) # all atoms that are already in a big node
        if not intersection == set(): # If atoms are shared, merge all possbile groups and update everything
            # for each intersected atom. find all the big nodes
            big_nodes_placement = set()
            for atom in intersection:
                big_nodes_placement.add(look_up.get(atom))
            # combine the full set
            new_big_node = data
            for p in big_nodes_placement:
                new_big_node = new_big_node.union(b_nodes[p])
                b_nodes.pop(p)    # remove all big nodes, pointed to by the intersect
            # for each atom, update and add the lookup dict for big nodes
            pos = len(b_nodes)
            for atom in new_big_node:
                look_up[atom] = pos
            # add the new big node
            b_nodes.append(new_big_node)
        else: # complete new big node
            # for each atom, add lockup dict for big nodes.
            pos = len(b_nodes)
            for atom in data:
                look_up[atom] = pos
            # insert the ring as a set in big nodes
            b_nodes.append(data)

        for atom in data: # update the atoms_in_big_nodes
            atom_list.add(atom)

    def core_ring(cut_molecule, look_up): # check for a ring, where core connects duo to the chemical reaction. Using recursive depth first search
        # use a dict to keep track on parents
        parent_dict = {}
        # Make a list of all atoms already visited and to look at
        visited_nodes = [0] # root has already been visited
        current_nodes = [0] # starting with root
        # Go over each atom starting from the core with breath first search
        while len(current_nodes) > 0:
            node_id = current_nodes[0]
            children = cut_molecule[node_id].children # get the kids of the molecule
            for child in children: # for each child, tjek if they have already been visted
                if child in visited_nodes: # WE HAVE FUND A RING!
                    ring_list = [0] # add all placements of the ring, core already added
                    # child and current node is the last edge in the ring. follow them up through the parent dict to the core
                    ring_id = node_id # for the first connecter in the ring
                    while ring_id != 0:
                        ring_list.append(ring_id) # add to ring list
                        ring_id = parent_dict.get(ring_id)
                    ring_id = child # for the second connecter in the ring
                    while ring_id != 0:
                        ring_list.append(ring_id) # add to ring list
                        ring_id = parent_dict.get(ring_id)
                    return ring_list
                else: # no ring add to visted and parent dict
                    parent_dict[child] = node_id
                    visited_nodes.append(child)
                    current_nodes.append(child) # add the child to the current node list
            current_nodes.pop(0)
        return None

    lookup = {} # look up dict
    cut_molecule = [] # The "atom" list for the cut molecule
    ### Find all double bonds, triple adn rings and insert the atoms together ###
    big_nodes = [core.copy()] # all the nodes with more than one id
    big_nodes_look_up = {} # atom maps to which big node its in
    atoms_in_big_nodes = core.copy() # all atoms in a set, which is in the big node

    # for each atom in core add to lookup dict
    for atom in core:
        big_nodes_look_up[atom] = 0

    # find alle rings first, and add them as big nodes.
    rings = GetSymmSSSR(rdk_mol)
    for r in rings:
        data = set()
        for atom in r:
            data.add(int(atom))
        big_node_update(big_nodes, big_nodes_look_up, atoms_in_big_nodes, data)

    # for each bond that is not a single bond. make a big node or add to already existing node
    bonds = []
    for bond in rdk_mol.GetBonds(): # find all bonds that are not single, and make a big node update
        start_atom = int(bond.GetBeginAtomIdx())
        end_atom = int(bond.GetEndAtomIdx())
        if not big_nodes_look_up.get(start_atom) == big_nodes_look_up.get(end_atom) or (big_nodes_look_up.get(end_atom) == None and big_nodes_look_up.get(start_atom) == None): # Not in the same big node or both not in any bignodes
            if not str(bond.GetBondType()) == "SINGLE": # make a big node if not a single bond
                big_node_update(big_nodes, big_nodes_look_up, atoms_in_big_nodes, {start_atom, end_atom})
            else:
                bonds.append(bond) # get the bonds that are not internal in a big node, to reduce later iteration over bonds

    # add all nodes to the cut molecule, based on big nodes and not big nodes
    node_holder = []
    for node in big_nodes:
        if not node.intersection(core) == set(): # its the core
            # add the core to the molecule and insert it on place zero
            for atom in node:
                lookup[atom] = 0
            cut_molecule.append(MoleculeNode(node, 1))
        else: # not the core
            # add the node to the placeholder for now, and add the look_up dict with a plus one. Since the core is the first node always
            for atom in node:
                lookup[atom] = len(node_holder) + 1
            node_holder.append(MoleculeNode(node, 0))
    # concat the to list of nodes to complete the cut moleucle, with the core on place zero
    cut_molecule = cut_molecule + node_holder

    # add all none big node atom to the cut molecule
    for atom in rdk_mol.GetAtoms():
        if not atom.GetIdx() in atoms_in_big_nodes: # the atom its not part of any big_node
            lookup[int(atom.GetIdx())] = len(cut_molecule)
            cut_molecule.append(MoleculeNode({int(atom.GetIdx())}, 0))

    # add all edges missing. One child layer at the time
    parent_list = cut_molecule[0].id
    while len(parent_list) > 0:
        # list of edges to work with
        new_parent = set()
        new_bonds = bonds.copy()
        # loop over every edge and parent to find matches
        for b in bonds:
            for p in parent_list:
                # check if a edge belongs to a parent
                if b.GetBeginAtomIdx() == p:
                    cut_molecule[lookup.get(p)].children.add(int(b.GetEndAtomIdx()))
                    if not big_nodes_look_up.get(int(b.GetEndAtomIdx())) == None:
                        for atom in big_nodes[big_nodes_look_up.get(int(b.GetEndAtomIdx()))]:
                            new_parent.add(atom)
                    else:
                        new_parent.add(int(b.GetEndAtomIdx()))
                    new_bonds.remove(b)
                    break
                elif b.GetEndAtomIdx() == p:
                    cut_molecule[lookup.get(p)].children.add(int(b.GetBeginAtomIdx()))
                    if not big_nodes_look_up.get(int(b.GetBeginAtomIdx())) == None:
                        for atom in big_nodes[big_nodes_look_up.get(int(b.GetBeginAtomIdx()))]:
                            new_parent.add(atom)
                    else:
                        new_parent.add(int(b.GetBeginAtomIdx()))
                    new_bonds.remove(b)
                    break
        parent_list = new_parent
        bonds = new_bonds

    #print("CORE RING: ")
    #print(core_ring(cut_molecule, lookup, [core], 0))

    return cut_molecule, lookup

# cut_molecule is the mocule to find cuts on
# Cuts is a set for all cuts already performed, this must include all atom ids even if part of a big node, hence flatten the list
# Lookup is a dict, for lookup placement with ids in the cut_molecule
# node is always 0, since its the placement of the root. This is due to the recursive nature of the function
def find_all_cuts(cut_molecule: [MoleculeNode], cuts: set, lookup: dict):
    # if node has no children return empty cuts list. This case should only happen if all atoms is the core
    all_childs_are_cut = True
    for child in cut_molecule[0].children: # check if the childs of the core is already cut
        if not child in cuts: # if any child is nok cut, break and search for a cut
            all_childs_are_cut = False
            break
    if all_childs_are_cut: # nothing to cut, return empty set
        return []
    else: # else find cuts return set with possible cuts
        return cut_search(cut_molecule, cuts, lookup, 0)

def cut_search(cut_molecule: [MoleculeNode], cuts: set, lookup: dict, node: int):
    # Check if node is a leaf based on different attributes
    def is_cut(node, none_leafs):
        cut_check = True
        for child in node.children:
            # if child has no childs or have been cut before.
            child_node = cut_molecule[lookup.get(child)]
            if len(child_node.children) > 0 and child not in cuts: # ordinary leaf
                cut_check = False
                none_leafs.append(child)
            elif len(child_node.children) == 0 and child not in cuts and len(child_node.id) > 1: # check if leaf is a big node
                cut_check = False
                none_leafs.append(child)
        return cut_check

    # check if possible cut
    none_leaf_childs = []
    new_cuts = []
    if is_cut(cut_molecule[node], none_leaf_childs):
        # if not root add cut
        if not cut_molecule[node].root:
            new_cuts.append(cut_molecule[node].id)
    # if not go over childs
    else:
        for c in none_leaf_childs:
            deeper_cuts = cut_search(cut_molecule, cuts, lookup, lookup.get(c))
            if len(deeper_cuts) > 0:
                new_cuts = new_cuts + deeper_cuts
    return new_cuts


def make_cut(mol, cuts, molecule, lookup_dict):
    """"takes an rdkit mol, the ids of atoms to cut off, the molecule object and lookup dictionary, returns modified mol and ordering"""
    # based on cuts to be performed decides which atoms are removed and which are replaced
    ban_list = set()
    replace_list = set()
    for c in cuts:
        ban_list.add(c)
        for child in molecule[lookup_dict.get(c)].children:
            ban_list.add(child)

    big_ban_list = set()
    big_ban_list.update(ban_list)
    big_ban_list.update(cuts)
    for cut in cuts:
        for neighbor_atom in mol.GetAtomWithIdx(cut).GetNeighbors():
            if neighbor_atom.GetIdx() not in big_ban_list:
                replace_list.add(cut)
    ban_list = [x for x in ban_list if x not in replace_list]

    # perform atom replacement and removal
    ordering = {}
    counter = 1
    atoms_to_remove = []
    atoms_to_compute_coordinates = []
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in ban_list:
            if atom.GetIdx() in replace_list:
                mol.ReplaceAtom(atom.GetIdx(), Atom("H"), updateLabel=True, preserveProps=False)
                atoms_to_compute_coordinates.append(atom)
            ordering[atom.GetIdx() + 1] = counter
            counter += 1
        else:
            atoms_to_remove.append(atom.GetIdx()) # adding id to the list of atoms needed to remove

    # perform bond removal
    bonds_to_remove = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetIdx() in ban_list or bond.GetEndAtom().GetIdx() in ban_list:
            bonds_to_remove.append(bond) # add to remove list
    bonds_to_remove.reverse()
    for bond in bonds_to_remove: # removing bond
        mol.RemoveBond(bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx())
    # remove atoms
    atoms_to_remove.reverse()
    for atom_id in atoms_to_remove:
        mol.RemoveAtom(atom_id)

    # recompute coordinates of replaced atoms
    for atom in atoms_to_compute_coordinates:
        try:
            rdmolops.SetTerminalAtomCoords(mol, atom.GetIdx(), atom.GetNeighbors()[0].GetIdx())
        except:
            pass
    return mol, ordering


def recompute_coordinates_of_mol(mol):
    try:
        mol.UpdatePropertyCache()
        mol.RemoveAllConformers()
        Compute2DCoords(mol)  # add coordinates with a comformer
        EmbedMolecule(mol)
    except AtomValenceException:
        print("atomvalencexception")
        pass
    return mol
