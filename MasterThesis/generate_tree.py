from mod import *
import openbabel.pybel as pybel
from openbabel import openbabel
from mod_to_xyz import mod_to_xyz


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


if __name__ == "__main__":

    out = product_to_gml("stringfile.xyz0000")
    print(out)
    gml = graphGMLString(out)
    mod_to_xyz(gml)
