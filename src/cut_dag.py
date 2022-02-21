from src.blackbox2 import run_gsm_cuts
from src.cut_molecule import make_cut_molecule, find_all_cuts, make_cut
from src.stringfile_helper_functions import mol_to_xyz
from src.stringfile_to_rdkit import stringfile_to_rdkit
from src.stringfile_tester import check_product, check_educt_to_product


class CutDagNode:
    def __init__(self, cuts):
        self.energy = []  # List for energy leves at different notes for the reaction
        self.RMS = 1  # the root mean square base on the original molecule reaction
        self.stringfile = ""  # The string file with the reaction, used for GML rule
        self.cuts = cuts  # which cuts on the molecule was made
        self.childs = []  # Childs made from the molecule

class CutDag:
    def __init__(self, cut_molecule, lookup_dict):
        self.layers = {} # Dictionary for each layer of the cut tree. Keys are ints of cuts made for a certain tree_layers
        self.cut_molecule = cut_molecule # the molecule to perfom cuts on
        self.cut_molecule_lookup_dict = lookup_dict # dict for cut_molecule

def removeDuplicates(arr): # midlertidig methode. lav core_ring_check i cut molecule for at fikse det.
    temp = []
    for e in arr:
        if e not in temp:
            temp.append(e)
    return temp

# generate childs of a node
def make_childs(tree: CutDag, node: CutDagNode, layer: int):
    # find all cuts på moleculet
    child_cuts = removeDuplicates(find_all_cuts(tree.cut_molecule, node.cuts, tree.cut_molecule_lookup_dict)) ################# fjern duplicate function når core ring er done #######################
    # generate all childs
    child_sets = []
    for cut in child_cuts:
        child_sets.append(node.cuts.union(cut))
    # return childs List
    child_nodes = []
    # chek if new layer exist
    if not layer+1 in tree.layers.keys():
        # add all nodes to the list and as childs in parent node
        for c in child_sets:
            node.childs.append(len(child_nodes)) # no -1 is needed, since we do it before adding he child
            child_nodes.append(CutDagNode(c))
        tree.layers[layer+1] = child_nodes
    else: # generate all childs and check if the exist before adding
        layer_list = tree.layers[layer+1]
        for c in child_sets:
            child_not_done = True
            # run over the layer in the tree
            for layer_node_placement in range(len(layer_list)):
                if layer_list[layer_node_placement].cuts == c:
                    # add child node to parent node and go to next child
                    node.childs.append(layer_node_placement)
                    child_not_done = False
                    break
            if child_not_done:
                # add child to parent and add child to list
                node.childs.append(len(layer_list))
                child = CutDagNode(c)
                layer_list.append(child)
                child_nodes.append(child)
    return child_nodes

# find child cuts and return them
def make_childs_mp(stringfile, cuts, placement): # stringfile for make cut molecule, cuts for what have already been cut on the molecule, placement contaisn l and p tha are the location in the cut dag being layer and placement
    # make cut molecule
    rdk_mol, atom_core, energy_curve = stringfile_to_rdkit(stringfile, False)
    molecule, lookup_dict = make_cut_molecule(rdk_mol, atom_core)
    # find cuts on molecule
    child_cuts = find_all_cuts(molecule, cuts, lookup_dict)
    # generate all child cuts
    child_sets = []
    for cut in child_cuts:
        child_sets.append(cuts.union(cut))

    return (child_sets, placement)

# take childs cuts and insert them as nodes in the cutdag, return a list of tasks
def insert_childs_mp(stringfile, cd, child_sets, placement):
    child_nodes = []
    task_list = []
    node = cd.layers[placement[0]][placement[1]]
    # placement knoewledge for tasks
    layer_nr = placement[0]+1
    placement_nr = []

    # chek if new layer exist
    if not layer_nr in cd.layers.keys():
        # add all nodes to the list and as childs in parent node
        for c in child_sets:
            node.childs.append(len(child_nodes)) # no -1 is needed, since we do it before adding he child
            placement_nr.append(len(child_nodes)) # get placements for tasks
            child_nodes.append(CutDagNode(c))
        cd.layers[layer_nr] = child_nodes
    else: # generate all childs and check if the exist before adding
        layer_list = cd.layers[layer_nr]
        for c in child_sets:
            child_not_done = True
            # run over the layer in the tree
            for layer_node_placement in range(len(layer_list)):
                if layer_list[layer_node_placement].cuts == c:
                    # add child node to parent node and go to next child
                    node.childs.append(layer_node_placement)
                    child_not_done = False
                    break
            if child_not_done:
                # add child to parent and add child to list
                node.childs.append(len(layer_list))
                placement_nr.append(len(layer_list)) # get placements for tasks
                child = CutDagNode(c)
                layer_list.append(child)
                child_nodes.append(child)
    # making tasks for processes
    if len(child_nodes) > 0:
        tasks = []
        for i in range(len(child_nodes)):
            tasks.append((make_childs_mp, (stringfile, child_nodes[i].cuts, (layer_nr, placement_nr[i]))))
        return tasks
    return []

def run_blackbox(stringfile, overall_folder, cuts, placement, reaction_folder):
    # make cut molecules
    rdk_mol, atom_core, energy_curve = stringfile_to_rdkit(stringfile, False)
    molecule, lookup_dict = make_cut_molecule(rdk_mol, atom_core)
    # make cuts on it
    modified_mol, order = make_cut(rdk_mol, cuts, molecule, lookup_dict)
    #modified_mol = recompute_coordinates_of_mol(modified_mol)
    xyz_file = mol_to_xyz(modified_mol)

    cut_folder = "/"    # make cut folder name
    cuts = sorted(cuts)           # sort cuts to make cuts_folder naming consistent
    for cut in cuts:
        cut_folder = cut_folder + str(cut) + "_"
    cut_folder = cut_folder[0:-1] + "/"

    # call true black box
    stringfile_path = run_gsm_cuts([xyz_file], overall_folder, reaction_folder, cut_folder, order)
    #stringfile_path = run_zstruct_and_gsm([xyz_file], overall_folder, order, atom_core, reaction_folder, cut_folder, logfile=True)
    if stringfile_path != "NO REACTION" and check_product(stringfile, stringfile_path, cuts, order, molecule, lookup_dict) and check_educt_to_product(stringfile_path): # check if reaction is the same
        return [stringfile_path, placement]
    else: # return data
        return ["NO REACTION", placement]

# generate root node
# input: Stringfile from Xtb, boolean for making visuals of the cut molecute
def make_root(stringfile: str, visuals: bool):
    # make cut molecule and energy_curve
    rdk_mol, atom_core, energy_curve = stringfile_to_rdkit(stringfile, visualize=visuals)

    if len(atom_core) == 0:
        return None
    molecule, lookup_dict = make_cut_molecule(rdk_mol, atom_core)
    # generate cut_tree
    ct = CutDag(molecule, lookup_dict)

    # generate node for root, and add data
    root = CutDagNode(set())
    root.energy = energy_curve
    root.stringfile = stringfile

    # Insert root in layers
    ct.layers[0] = [root]

    return ct
