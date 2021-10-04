from mod import *
from generate_tree import reaction_and_product_to_gml


class CutTreeNode:
    def __init__(self, molecule, cuts):
        self.energy = []  # List for energy leves at different notes for the reaction
        self.RMS = 0  # the root mean square base on the original molecule reaction
        self.stringfile = ""  # The string file with the reaction, used for GML rule
        self.cuts = cuts  # which cuts on the molecule was made
        self.childs = []  # Childs made from the molecule


class MoleculeNode:
    def __init__(self, molecule_id, node_type):
        self.id = molecule_id
        self.children = []
        self.root = node_type


# add a edge from parent to child in cut molecule
def add_child(child_id, core, parent_molecule, parents_list):
    # alter target id
    if child_id < core[0][0]:
        place_id = child_id + 1
    else:
        place_id = child_id + 1 - len(core[0])
    parent_molecule.children.append(place_id)
    parents_list.append(([child_id], place_id))


# core = [[atom ID's],[edge ID's]]
def make_cut_molecule(g_mod, core):
    cut_molecule = [MoleculeNode(core[0], 1)]
    # Insert core
    # insert all nodes
    for v in g_mod.vertices:
        if v.id not in core[0]:
            cut_molecule.append(MoleculeNode([v.id], 0))
    # remove internal core edges
    edges = []
    edge_counter = 0
    for e in g_mod.edges:
        if edge_counter not in core[1]:
            edges.append(e)
        edge_counter += 1

    # add all edges missing one child layer at the time
    parent_list = [(cut_molecule[0].id, 0)]
    while len(parent_list) > 0:

        new_edges = []
        new_parent = []
        for e in edges:
            for p in parent_list:
                if e.source.id in p[0]:
                    add_child(e.target.id, core, cut_molecule[p[1]], new_parent)
                elif e.target.id in p[0]:
                    add_child(e.source.id, core, cut_molecule[p[1]], new_parent)
                else:
                    new_edges.append(e)
        parent_list = new_parent
        edges = new_edges

    return cut_molecule

    # calculate how to shift everything else
    # insert all new atoms in molecule list


def find_all_cuts(cut_molecule: [MoleculeNode], cuts: set):
    def is_leaf(node):
        for child in node.children:
            if cut_molecule[child].children or cut_molecule[child].id not in cuts:
                return False
        return True




if __name__ == "__main__":
    gml, ac, bc = reaction_and_product_to_gml('stringfile.xyz0000', visualize=False)
    g = graphGMLString(gml)

    m = make_cut_molecule(g, [ac, bc])
    for n in m:
        print("id: " + str(n.id))
        print("bond: " + str(n.children))
