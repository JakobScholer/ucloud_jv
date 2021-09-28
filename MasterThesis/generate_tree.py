from mod import *
import openbabel.pybel as pybel
from openbabel import openbabel
from mod_to_xyz import mod_to_xyz
from igraph import *
import plotly.graph_objects as go


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
        bmap[(src,tar)] = order_map[b.GetBondOrder()]
    return bmap

'''
def product_to_gml(filename):
    with open(filename) as f:
        content = f.readlines()
    num_atoms = int(content[0])
    xyz_str: str = "".join(content[len(content) - (num_atoms + 2):])
    print(xyz_str)

    product = pybel.readstring("xyz", xyz_str)

    verts = []
    for i in range(num_atoms):
        a1 = product.atoms[i]
        symbol = openbabel.GetSymbol(a1.atomicnum)
        a1_lbl = symbol
        if a1.formalcharge != 0:
            chg = "+" if a1.formalcharge > 0 else "-"
            a1_lbl = a1_lbl + chg * a1.formalcharge

        verts.append(f'node [ id {i} label "{a1_lbl}" ]')

    edges = []

    order_map = {
        1: "-",
        2: "=",
        3: "#",
        1.5: ":"
    }
    bmap = {}
    b: openbabel.OBBond
    for b in openbabel.OBMolBondIter(product.OBMol):
        src, tar = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if src > tar:
            src, tar = tar, src
        bmap[(src, tar)] = order_map[b.GetBondOrder()]

    for (src, tar), ob_bond in bmap.items():
        edges.append(f'edge [ source {src-1} target {tar-1} label "{ob_bond}"]')

    verts_str = "\n    ".join(verts)
    edges_str = "\n    ".join(edges)
    #ruleID "{mol_reactant.formula}: {float(reaction.product.energy) - float(reaction.reactant.energy)}"
    gml_str = f"""
    graph [
    {verts_str}
    {edges_str}
    ]
    """
    return gml_str
'''


def reaction_to_gml(filename):
    with open(filename) as f:
        content = f.readlines()
    num_atoms = int(content[0])
    xyz_str_reactant: str = "".join(content[:(num_atoms + 2)])
    xyz_str_product: str = "".join(content[len(content) - (num_atoms + 2):])

    reactant = pybel.readstring("xyz", xyz_str_reactant)
    product = pybel.readstring("xyz", xyz_str_product)
    print(reactant)
    print(product)

    num_atoms: int = len(reactant.atoms)
    assert(num_atoms == len(product.atoms))

    left_verts, ctx_verts, right_verts = [], [], []
    for i in range(num_atoms):
        a1, a2 = reactant.atoms[i], product.atoms[i]
        assert(a1.idx == a2.idx and a1.atomicnum == a2.atomicnum)
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
    bond_core = []

    for (src, tar), ob_bond in bmap1.items():
        if (src, tar) not in bmap2:
            left_edges.append(f'edge [ source {src-1} target {tar-1} label "{ob_bond}"]')
        else:
            if ob_bond != bmap2[(src, tar)]:
                if (src-1) not in atom_core:
                    atom_core.append(src-1)
                if (tar-1) not in atom_core:
                    atom_core.append(tar-1)
                bond_core.append([src-1, tar-1, ob_bond])
                left_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')
                right_edges.append(f'edge [ source {src-1} target {tar-1} label "{bmap2[(src, tar)]}"]')
            else:
                ctx_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')

    for (src, tar), ob_bond in bmap2.items():
        if (src, tar) in bmap1:
            continue
        right_edges.append(f'edge [ source {src-1} target {tar-1} label "{ob_bond}"]')
    print("CORE:")
    print(atom_core)
    print(bond_core)
    left_verts_str = "\n\t    ".join(left_verts)
    ctx_verts_str = "\n\t    ".join(ctx_verts)
    right_verts_str = "\n\t    ".join(right_verts)
    left_edges_str = "\n\t    ".join(left_edges)
    ctx_edges_str = "\n\t    ".join(ctx_edges)
    right_edges_str = "\n\t    ".join(right_edges)
    #ruleID "{mol_reactant.formula}: {float(reaction.product.energy) - float(reaction.reactant.energy)}"
    gml_str = f"""
    rule [
        ruleID "{reactant.formula}: {product.energy}"

        left [
            {left_verts_str}
            {left_edges_str}
        ]
        context [
            {ctx_verts_str}
            {ctx_edges_str}
        ]
        right [
            {right_verts_str}
            {right_edges_str}
        ]
    ]
    """
    return gml_str


if __name__ == "__main__":
    #mol = smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
    #mod_to_xyz(mol)

    # Product
    p = reaction_to_gml("stringfile.xyz0000")
    print(p)

    g = Graph()

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Xe,
                             y=Ye,
                             mode='lines',
                             line=dict(color='rgb(210,210,210)', width=1),
                             hoverinfo='none'
                             ))
    fig.add_trace(go.Scatter(x=Xn,
                             y=Yn,
                             mode='markers',
                             name='bla',
                             marker=dict(symbol='circle-dot',
                                         size=18,
                                         color='#6175c1',  # '#DB4551',
                                         line=dict(color='rgb(50,50,50)', width=1)
                                         ),
                             text=labels,
                             hoverinfo='text',
                             opacity=0.8
                             ))

    #gml = graphGMLString(out)
    #mod_to_xyz(gml)
