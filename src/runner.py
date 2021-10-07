from mod import *
from src.generate_tree import reaction_and_product_to_gml, fig_plot


class CutTreeNode:
    def __init__(self, molecule, cuts):
        self.energy = []        # List for energy leves at different notes for the reaction
        self.RMS = 0            # the root mean square base on the original molecule reaction
        self.stringfile = ""    # The string file with the reaction, used for GML rule
        self.cuts = cuts        # which cuts on the molecule was made
        self.childs = []        # Childs made from the molecule


class MoleculeNode:
    def __init__(self, molecule_id, node_type):
        self.id = molecule_id   # List of atom ID's
        self.children = set()      # A list of ints representing the list placement of the children
        self.root = node_type   # 1 = root, 0 = not root


# core = [[atom ID's],[edge ID's]]
def make_cut_molecule(g_mod, core):
    lookup = {}
    for c in core:
        lookup[c] = 0
    cut_molecule = [MoleculeNode(core, 1)]
    # insert all nodes
    for v in g_mod.vertices:
        if v.id not in core:
            lookup[v.id] = len(cut_molecule)
            cut_molecule.append(MoleculeNode([v.id], 0))
    # remove internal core edges
    edges = []
    for e in g_mod.edges:
        # only add edge if its not in the core
        if not e.source.id in core or not e.target.id in core:
            edges.append(e)
    # add all edges missing one child layer at the time
    parent_list = cut_molecule[0].id
    while len(parent_list) > 0:
        # list of edges to work with
        new_parent = []
        new_edges = edges.copy()
        # loop over every edge and parent to find matches
        for e in edges:
            for p in parent_list:
                # check if a edge belongs to a parent
                if e.source.id == p:
                    cut_molecule[lookup.get(p)].children.add(e.target.id)
                    new_parent.append(e.target.id)
                    new_edges.remove(e)
                    break
                elif e.target.id == p:
                    cut_molecule[lookup.get(p)].children.add(e.source.id)
                    new_parent.append(e.source.id)
                    new_edges.remove(e)
                    break
        parent_list = new_parent
        edges = new_edges
    return cut_molecule, lookup


def find_all_cuts(cut_molecule: [MoleculeNode], cuts: set, lookup: dict, node: int):
    # Check if node is a leaf based on different attributes
    def is_cut(nod, none_leafs):
        cut_check = True
        for child in nod.children:
            # if child has no childs or have been cut before.
            if len(cut_molecule[lookup.get(child)].children) > 0 and child not in cuts:
                cut_check = False
                none_leafs.append(child)
        return cut_check

    # if node has no children return empty cuts list. This case should only happen if all atoms is the core
    if not cut_molecule[node].children:
        return True

    # check if possible cut
    none_leaf_childs = []
    if is_cut(cut_molecule[node], none_leaf_childs):
        # if not root add cut
        if not cut_molecule[node].root:
            cuts.add(cut_molecule[node].id[0])
    # if not go over childs
    else:
        for c in none_leaf_childs:
            find_all_cuts(cut_molecule, cuts, lookup, lookup.get(c))
    return True


def runner_main():
    #fig_plot('src/gmlstring.gml', [1, 3])
    with open('src/gmlstring.gml', 'r') as file:
        gml = file.read()
    g = graphGMLString(gml)
    m, l = make_cut_molecule(g, [13, 8, 6, 5, 4, 11, 12])
    #print(l)
    for n in m:
        print("id: " + str(n.id))
        print("bond: " + str(n.children))

    cuts = set()
    find_all_cuts(m, cuts, l, 0)
    print(cuts)
    '''
    gml, ac, ec = reaction_and_product_to_gml('src/stringfile.xyz0000', visualize=True)
    g = graphGMLString(gml)

    m = make_cut_molecule(g, ac)
    for n in m:
        print("id: " + str(n.id))
        print("bond: " + str(n.children))

    cuts = set()
    find_all_cuts(m, cuts, 0)
    print(cuts)
    '''


