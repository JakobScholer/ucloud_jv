from openbabel import openbabel
from rdkit.Chem import BondType


def build_bond_map(mol):
    """Takes an rdkit mol object, returns a bond mapping for all atoms in mol"""
    order_map = {
        1: BondType.SINGLE,
        2: BondType.DOUBLE,
        3: BondType.TRIPLE,
        1.5: BondType.AROMATIC
    }
    bmap = {}
    b: openbabel.OBBond
    for b in openbabel.OBMolBondIter(mol.OBMol):
        src, tar = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if src > tar:
            src, tar = tar, src
        bmap[(src, tar)] = order_map[b.GetBondOrder()]
    return bmap