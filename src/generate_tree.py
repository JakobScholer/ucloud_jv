from mod import *
import openbabel.pybel as pybel
from openbabel import openbabel
from mod_to_xyz import mod_to_xyz
from igraph import *
import plotly.graph_objects as go

from root_mean_square import root_mean_square


def fig_plot(gmlfile, core_atoms):
    g = Graph.Read_GML(gmlfile)
    labels = list(g.vs['label'])
    edge_label_list = list(g.es['label'])
    vertex_amount = len(labels)
    edge_list = [e.tuple for e in g.es]  # list of edges
    core_bonds = []
    edge_counter = 0
    for e in edge_list:
        if e[0] in core_atoms and e[1] in core_atoms:
            core_bonds.append(edge_counter)
        edge_counter += 1
    layt = g.layout('kk')  # kamada-kawai layout

    atom_x = [layt[k][0] for k in range(vertex_amount)]  # x coordinate of atoms
    atom_y = [layt[k][1] for k in range(vertex_amount)]  # y coordinate of atoms
    atom_core_x = []
    atom_core_y = []
    for atom in core_atoms:
        atom_core_x.append(atom_x[atom])
        atom_core_y.append(atom_y[atom])
    bond_x = []  # x coordinates of line (start and end)
    bond_y = []  # y coordinates of line (start and end)
    bond_middle_x = []  # x coordinate of middle of line
    bond_middle_y = []  # y coordinate of middle of line
    for edge in edge_list:
        bond_x += [layt[edge[0]][0], layt[edge[1]][0], None]
        bond_y += [layt[edge[0]][1], layt[edge[1]][1], None]
        bond_middle_x.append((float(layt[edge[0]][0]) + float(layt[edge[1]][0])) / 2)
        bond_middle_y.append((float(layt[edge[0]][1]) + float(layt[edge[1]][1])) / 2)
    bond_core_x = []
    bond_core_y = []
    for bond in core_bonds:
        bond_core_x.append(bond_middle_x[bond])
        bond_core_y.append(bond_middle_y[bond])

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=atom_core_x,
                             y=atom_core_y,
                             mode='markers',
                             name='core atoms',
                             marker=dict(symbol='circle-dot',
                                         size=25,
                                         color='#cf0202'
                                         ),
                             hoverinfo='skip',
                             ))

    fig.add_trace(go.Scatter(x=bond_x,
                             y=bond_y,
                             mode='lines',
                             name='bonds',
                             line=dict(color='rgb(210,210,210)', width=1),
                             text=edge_label_list,
                             hoverinfo='skip'
                             ))
    fig.add_trace(go.Scatter(x=bond_middle_x,
                             y=bond_middle_y,
                             mode='text',
                             name='bondIDs',
                             text=list(range(0, len(edge_label_list))),
                             hoverinfo='text',
                             textfont_size=1
                             ))
    fig.add_trace(go.Scatter(x=bond_core_x,
                             y=bond_core_y,
                             mode='markers',
                             name='core bonds',
                             marker=dict(symbol='circle-dot',
                                         size=25,
                                         color='#cf0202'
                                         ),
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=bond_middle_x,
                             y=bond_middle_y,
                             mode='text',
                             name='bondtypes',
                             text=edge_label_list,
                             hoverinfo='skip',
                             textfont_size=20
                             ))
    fig.add_trace(go.Scatter(x=atom_x,
                             y=atom_y,
                             mode='text',
                             name='atomIDs',
                             text=list(range(0, len(labels))),
                             hoverinfo='text',
                             textfont_size=15
                             ))
    fig.add_trace(go.Scatter(x=atom_x,
                             y=atom_y,
                             mode='markers+text',
                             name='atoms',
                             marker=dict(symbol='circle-dot',
                                         size=18,
                                         color='#61c1ab'
                                         ),
                             text=labels,
                             hoverinfo='skip',
                             ))
    fig.show()


def build_bond_map(mol):
    order_map = {
        1: "-",
        2: "=",
        3: "#",
        1.5: ":"
    }
    bmap = {}
    b: openbabel.OBBond
    for b in openbabel.OBMolBondIter(mol.OBMol):
        src, tar = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if src > tar:
            src, tar = tar, src
        bmap[(src, tar)] = order_map[b.GetBondOrder()]
    return bmap


# takes the content of a stringfile and puts all energy profiles into a list
def read_energy_profiles(string_content):
    energy_profiles = []
    for i in range(1, len(string_content), (int(string_content[0]) + 2)):
        energy_profiles.append(float(string_content[i]))
    return energy_profiles


def reaction_and_product_to_gml(filename, visualize=False):
    with open(filename) as f:
        content = f.readlines()
    num_atoms = int(content[0])
    xyz_str_reactant: str = "".join(content[:(num_atoms + 2)])
    xyz_str_product: str = "".join(content[len(content) - (num_atoms + 2):])

    # read all energy profiles
    energy_profiles = read_energy_profiles(content)

    reactant = pybel.readstring("xyz", xyz_str_reactant)
    product = pybel.readstring("xyz", xyz_str_product)

    num_atoms: int = len(reactant.atoms)
    assert (num_atoms == len(product.atoms))

    left_verts, ctx_verts, right_verts = [], [], []
    for i in range(num_atoms):
        a1, a2 = reactant.atoms[i], product.atoms[i]
        assert (a1.idx == a2.idx and a1.atomicnum == a2.atomicnum)
        symbol = openbabel.GetSymbol(a1.atomicnum)
        a1_lbl = symbol
        a2_lbl = symbol
        if a1.formalcharge != 0:
            chg = "+" if a1.formalcharge > 0 else "-"
            a1_lbl = a1_lbl + chg * a1.formalcharge
        if a2.formalcharge != 0:
            chg = "+" if a2.formalcharge > 0 else "-"
            a2_lbl = a2_lbl + chg * a2.formalcharge

        if a1_lbl != a2_lbl:
            left_verts.append(f'node [ id {i} label "{a1_lbl}" ]')
            right_verts.append(f'node [ id {i} label "{a2_lbl}" ]')
        else:
            ctx_verts.append(f'node [ id {i} label "{a1_lbl}" ]')

    left_edges, ctx_edges, right_edges = [], [], []

    bmap1 = build_bond_map(reactant)
    bmap2 = build_bond_map(product)
    atom_core = []

    for (src, tar), ob_bond in bmap1.items():
        if (src, tar) not in bmap2:
            left_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')
        else:
            if ob_bond != bmap2[(src, tar)]:
                if (src - 1) not in atom_core:
                    atom_core.append(src - 1)
                if (tar - 1) not in atom_core:
                    atom_core.append(tar - 1)
                left_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')
                right_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{bmap2[(src, tar)]}"]')
            else:
                ctx_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')

    for (src, tar), ob_bond in bmap2.items():
        if (src, tar) in bmap1:
            continue
        right_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')

    verts_str = "\n    ".join(left_verts + ctx_verts)
    edges_str = "\n    ".join(left_edges + ctx_edges)
    gml_str = f"""graph [
    {verts_str}
    {edges_str}
    ]
        """
    if visualize:
        f = open("gmlstring.gml", "w")
        f.write(gml_str)
        f.close()
        fig_plot('gmlstring.gml', atom_core)

    return gml_str, atom_core, energy_profiles


if __name__ == "__main__":
    gml, ac, ep = reaction_and_product_to_gml('stringfile.xyz0000', visualize=True)
    with open('stringfile.xyz0002') as fi:
        ct = fi.readlines()
    curve = read_energy_profiles(ct)
    x = root_mean_square(ep, curve)
    print(x)
