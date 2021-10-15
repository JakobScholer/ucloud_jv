from root_mean_square import root_mean_square
from runner import runner_main, make_cut_molecule, find_all_cuts, make_cut
from generate_tree import generate_tree_main, reaction_and_product_to_gml, read_energy_profiles
from mod_to_xyz import mod_to_xyz_main, mod_to_xyz

class CutTreeNode:
    def __init__(self, cuts):
        self.energy = []  # List for energy leves at different notes for the reaction
        self.RMS = 0  # the root mean square base on the original molecule reaction
        self.stringfile = ""  # The string file with the reaction, used for GML rule
        self.cuts = cuts  # which cuts on the molecule was made
        self.childs = []  # Childs made from the molecule

class CutTree:
    def __init__(self, cut_molecule, lookup_dict):
        self.tree_layers = {} # Dictionary for each layer of the cut tree. Keys are ints of cuts made for a certain tree_layers
        self.cut_molecule = cut_molecule # the molecule to perfom cuts on
        self.cut_molecule_lookup_dict = lookup_dict # dict for cut_molecule

# generate childs of a nodes
def make_childs(node: CutTreeNode, tree: CutTree)
    # find alle cuts på moleculet
    find_all_cuts(tree.cut_molecule, node.cuts, tree.cut_molecule_lookup_dict, 0)
    # opret alle børn
        # tjek hvis laget findes
        # Tjek hvis childet findes
        # indsæt børn på layer
        # indsæt childs i nodes childs
    # done?

# Check if child already exist

# generate root node
# input: Stringfile from Xtb, boolean for making visuals of the cut molecute
def make_root(stringfile: string, visuals: boolean):
    # make cut molecule and energy_curve
    gml_string, atom_core, energy_curve = reaction_and_product_to_gml(stringfile, visualize=visuals)
    g = graphGMLString(gml_string)
    molecule, lookup_dict = make_cut_molecule(g, atom_core)
    # generate cut_tree
    ct = CutTree(molecule, lookup_dict)

    # generate node for root, and add data
    root = CutTreeNode(set())
    root.energy = energy_curve
    root.stringfile = stringfile

    # Insert root in layers
    ct.tree_layers[0] = [root]

    return ct

def cut_tree_main():
    # for testing
    'src/stringfile.xyz0000'

if __name__ == "__main__":
    cut_tree_main()
