from openbabel import pybel, openbabel
from src.stringfile_helper_functions import build_bond_map, read_stringfile_content


def stringfile_to_gml(filename: str):
    """"Takes the name of a stringfile, returns a gml rule representing it"""
    educt_str, product_str, num_atoms = read_stringfile_content(filename)
    reactant = pybel.readstring("xyz", educt_str)   # reads educt as openBabel molecule object
    product = pybel.readstring("xyz", product_str)  # reads product as openBabel molecule object

    left_verts, ctx_verts, right_verts = [], [], []
    for i in range(num_atoms):
        a1, a2 = reactant.atoms[i], product.atoms[i]
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
    for (src, tar), ob_bond in bmap1.items():
        if (src, tar) not in bmap2:
            left_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')
        else:
            if ob_bond != bmap2[(src, tar)]:
                left_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')
                right_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{bmap2[(src, tar)]}"]')
            else:
                ctx_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')

    for (src, tar), ob_bond in bmap2.items():
        if (src, tar) in bmap1:
            continue
        right_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')

    left_verts_str = "\n\t    ".join(left_verts)
    ctx_verts_str = "\n\t    ".join(ctx_verts)
    right_verts_str = "\n\t    ".join(right_verts)
    left_edges_str = "\n\t    ".join(left_edges)
    ctx_edges_str = "\n\t    ".join(ctx_edges)
    right_edges_str = "\n\t    ".join(right_edges)
    gml_str = f"""
    rule [
        ruleID "{reactant.formula}"

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
