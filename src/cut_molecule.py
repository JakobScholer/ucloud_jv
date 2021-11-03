from mod import *
from src.generate_tree import reaction_and_product_to_gml, fig_plot

class MoleculeNode:
    def __init__(self, molecule_id, node_type):
        self.id = molecule_id  # List of atom ID's
        self.children = set()  # A list of ints representing the list placement of the children
        self.root = node_type  # 1 = root, 0 = not root

# core = [[atom ID's],[edge ID's]]
# g_mod is the graph object from MØD
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

# cut_molecule is the mocule to find cuts on
# Cuts is a set for all cuts already performed
# Lookup is a dict, for lookup placement with ids in the cut_molecule
# node is always 0, since its the placement of the root. This is due to the recursive nature of the function
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
        return set()

    # check if possible cut
    none_leaf_childs = []
    new_cuts = set()
    if is_cut(cut_molecule[node], none_leaf_childs):
        # if not root add cut
        if not cut_molecule[node].root:
            new_cuts.add(cut_molecule[node].id[0])
    # if not go over childs
    else:
        for c in none_leaf_childs:
            deeper_cuts = find_all_cuts(cut_molecule, cuts, lookup, lookup.get(c))
            if len(deeper_cuts) > 0:
                new_cuts = new_cuts.union(deeper_cuts)
    return new_cuts


def make_cut(mod_graph, molecule_to_cut, molecules):
    ban_list = []
    replace_list = []
    for mc in molecule_to_cut:
        for c in molecules[mc + 1].children:
            ban_list.append(c)
        replace_list.append(mc)
    print("--------------før---------------")
    print(mod_graph.getGMLString())
    gml_string = "graph [\n"
    ordering = []
    for vertex in mod_graph.vertices:
        if vertex.id not in ban_list:
            if vertex.id in replace_list:
                gml_string += "    node [ id " + str(vertex.id) + " label \"" + "H" + "\" ]\n"
                ordering.append(vertex.id)
            else:
                gml_string += "    node [ id " + str(vertex.id) + " label \"" + str(vertex.stringLabel) + "\" ]\n"
                ordering.append(vertex.id)
    for edge in mod_graph.edges:
        if edge.source.id not in ban_list and edge.target.id not in ban_list:
            gml_string += "    edge [ source " + str(edge.source.id) + " target " + str(
                edge.target.id) + " label \"" + str(edge.bondType) + "\"]\n"
    print("--------------efter---------------")
    print(gml_string)
    print(ordering)
    return gml_string, ordering


def cut_molecule_main():
    gml, atom_core, ep = reaction_and_product_to_gml('test/testfiles/stringfile.xyz0143', visualize=True)
    g = graphGMLString(gml)
    m, l = make_cut_molecule(g, atom_core)
    for n in m:
        print("id: " + str(n.id))
        print("bond: " + str(n.children))
    cuts = find_all_cuts(m, set(), l, 0)
    print(cuts)
    gml, order = make_cut(g, cuts, m)
