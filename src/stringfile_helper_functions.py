from openbabel import openbabel
from rdkit.Chem import BondType


def read_stringfile_content(filename: str):
    """"Takes the name of a stringfile, returns the educt, product and number of atoms"""
    with open(filename) as f:
        content = f.readlines()
    num_atoms = int(content[0])                             # number of atoms in molecule
    educt = "".join(content[:(num_atoms + 2)])              # reads info on all atoms in educt molecule
    product = "".join(content[len(content) - (num_atoms + 2):]) # reads info on all atoms in product molecule
    return educt, product, num_atoms


def read_energy_profiles(filename: str):
    """Takes the name of a stringfile, returns a list of all energy values within it."""
    with open(filename) as f:
        content = f.readlines()
    energy_profiles = []
    for i in range(1, len(content), (int(content[0]) + 2)):
        energy_profiles.append(float(content[i]))
    return energy_profiles


def build_bond_map(mol):
    """Takes an rdkit mol object, returns a bond mapping for all atoms in mol object"""
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