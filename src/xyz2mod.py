import openbabel.pybel as pybel
from openbabel import openbabel
from typing import List, Union, IO
import os


class XYZMolecule:
    def __init__(self,  energy: float, xyz_str: str):
        self.energy: float = energy
        self.obmol: pybel.Molecule = pybel.readstring("xyz", xyz_str)


class XYZReactionPath:
    def __init__(self, name: str, nodes: List[XYZMolecule]):
        self.name = name
        self.nodes = nodes


def read_stringfile(fn: str):
    with open(fn) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    res: List[XYZMolecule] = []
    read_head: int = 0
    while read_head < len(content):
        num_atoms: int = int(content[read_head])
        energy: float = float(content[read_head+1])
        xyz_str: str = "\n".join(content[read_head:read_head+2+num_atoms])
        res.append(XYZMolecule(energy, xyz_str))
        read_head += 2 + num_atoms
    name = fn[-4:]
    return XYZReactionPath(name, res)


def read_stringfiles(path):
    _, _, filenames = next(os.walk(path))

    res = []
    for fn in filenames:
        if not fn.startswith("stringfile"):
            continue
        res.append(read_stringfile(os.path.join(path, fn)))
    return res


def reaction_to_gml(reactant: XYZMolecule, product: XYZMolecule, rule_cleaner: []):
    num_atoms: int = len(reactant.obmol.atoms)
    assert(num_atoms == len(product.obmol.atoms))

    left_verts, ctx_verts, right_verts = [], [], []
    for i in range(num_atoms):
        a1, a2 = reactant.obmol.atoms[i], product.obmol.atoms[i]
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

    bmap1 = build_bond_map(reactant.obmol)
    bmap2 = build_bond_map(product.obmol)
    for (src, tar), ob_bond in bmap1.items():
        if (src, tar) not in bmap2:
            left_edges.append(f'edge [ source {src-1} target {tar-1} label "{ob_bond}"]')
        else:
            if ob_bond != bmap2[(src, tar)]:
                left_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')
                right_edges.append(f'edge [ source {src-1} target {tar-1} label "{bmap2[(src, tar)]}"]')
            else:
                ctx_edges.append(f'edge [ source {src - 1} target {tar - 1} label "{ob_bond}"]')

    for (src, tar), ob_bond in bmap2.items():
        if (src, tar) in bmap1:
            continue
        right_edges.append(f'edge [ source {src-1} target {tar-1} label "{ob_bond}"]')

    left_verts_str = "\n\t    ".join(left_verts)
    ctx_verts_str = "\n\t    ".join(ctx_verts)
    right_verts_str = "\n\t    ".join(right_verts)
    left_edges_str = "\n\t    ".join(left_edges)
    ctx_edges_str = "\n\t    ".join(ctx_edges)
    right_edges_str = "\n\t    ".join(right_edges)
    #ruleID "{mol_reactant.formula}: {float(reaction.product.energy) - float(reaction.reactant.energy)}"
    gml_str = f"""
    rule [
        ruleID "{reactant.obmol.formula}: {product.energy}"

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
    out = read_stringfile("stringfile.xyz0000")
    rule_str = reaction_to_gml(out.nodes[0], out.nodes[-1], [])
    print(rule_str)
