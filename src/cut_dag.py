from mod import *
from src.root_mean_square import root_mean_square
from src.cut_molecule import cut_molecule_main, make_cut_molecule, find_all_cuts, make_cut
from src.generate_tree import generate_tree_main, reaction_and_product_to_gml, read_energy_profiles
from src.mod_to_xyz import mod_to_xyz_main, mod_to_xyz

class CutDagNode:
    def __init__(self, cuts):
        self.energy = []  # List for energy leves at different notes for the reaction
        self.RMS = 0  # the root mean square base on the original molecule reaction
        self.stringfile = ""  # The string file with the reaction, used for GML rule
        self.cuts = cuts  # which cuts on the molecule was made
        self.childs = []  # Childs made from the molecule

class CutDag:
    def __init__(self, cut_molecule, lookup_dict):
        self.layers = {} # Dictionary for each layer of the cut tree. Keys are ints of cuts made for a certain tree_layers
        self.cut_molecule = cut_molecule # the molecule to perfom cuts on
        self.cut_molecule_lookup_dict = lookup_dict # dict for cut_molecule

# generate childs of a node
def make_childs(node: CutDagNode, tree: CutDag):
    # find alle cuts på moleculet
    child_cuts = find_all_cuts(tree.cut_molecule, node.cuts, tree.cut_molecule_lookup_dict, 0)
    # generate all childs
    child_sets = []
    for cut in child_cuts:
        child_sets.append(node.cuts.union({cut}))
    # return childs List
    child_nodes = []
    # chek if new layer exist. MUTEX
    if not len(node.cuts)+1 in tree.layers.keys():
        # add all nodes to the list and as childs in parent node
        for c in child_sets:
            node.childs.append(len(child_nodes)) # no -1 is needed, since we do it before adding he child
            child_nodes.append(CutDagNode(c))
        tree.layers[len(node.cuts)+1] = child_nodes
    else: # generate all childs and check if the exist before adding
        layer_list = tree.layers[len(node.cuts)+1]
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
    gml_string, atom_core, energy_curve = reaction_and_product_to_gml(stringfile, visualize=visuals)
    g = graphGMLString(gml_string)
    molecule, lookup_dict = make_cut_molecule(g, atom_core)

    # find cuts on molecule
    child_cuts = find_all_cuts(molecule, cuts, lookup_dict, 0)
    # generate all child cuts
    child_sets = []
    for cut in child_cuts:
        child_sets.append(cuts.union({cut}))

    return (child_sets, placement)

# take childs cuts and insert them as nodes in the cutdag, return a list of tasks
def insert_childs_mp(stringfile, cd, childs_sets, placement):
    child_nodes = []
    task_list = []
    node = cd.layers[placement[0]][placement[1]]
    # placement knoewledge for tasks
    layer_nr = len(node.cuts)+1
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
            tasks.append((make_childs_mp, (stringfile, child_nodes[i].cuts, layer_nr, placement_nr[i])))
    return tasks

# generate root node
# input: Stringfile from Xtb, boolean for making visuals of the cut molecute
def make_root(stringfile: str, visuals: bool):
    # make cut molecule and energy_curve
    gml_string, atom_core, energy_curve = reaction_and_product_to_gml(stringfile, visualize=visuals)
    g = graphGMLString(gml_string)
    molecule, lookup_dict = make_cut_molecule(g, atom_core)
    # generate cut_tree
    ct = CutDag(molecule, lookup_dict)

    # generate node for root, and add data
    root = CutDagNode(set())
    root.energy = energy_curve
    root.stringfile = stringfile

    # Insert root in layers
    ct.layers[0] = [root]

    return ct

def cut_dag_main():
    # for testing
    test = 'test/testfiles/stringfile_ring.xyz0000'
    graph = False
    tree = make_root(test, graph)

    # find root og check op på den
    print("first test")
    print(tree.layers[0][0].stringfile)
    print(tree.layers[0][0].energy)
    print(tree.layers[0][0].cuts)
    # Lav børn på root og chek op på dem
    print("Second test")
    make_childs(tree.layers[0][0],tree)
    for node in tree.layers[1]:
        print(node.cuts)
        make_childs(node,tree)
    print("Second test")
    for node in tree.layers[2]:
        print(node.cuts)
    # lav børn på børnene :D det skal nok blive magisk

if __name__ == "__main__":
    cut_dag_main()
