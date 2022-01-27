import openbabel.pybel as pybel
from openbabel import openbabel
from rdkit.Chem import RWMol, MolFromSmiles, Atom, Conformer
from rdkit.Geometry import Point3D
from src.stringfile_helper_functions import build_bond_map, read_stringfile_content, read_energy_profiles
from src.visualizers import visualize_rdkit_mol

def stringfile_to_rdkit(filename: str, visualize: bool = False):
    """takes a stringfile, returns a rdkit mol object, the core of the atom and the energy profile."""
    xyz_str_reactant, xyz_str_product, num_atoms = read_stringfile_content(filename)

    # find all atom coordinates and store in list
    atoms = xyz_str_reactant.split("\n")
    coordinates = []
    for atom in atoms:
        coord = atom.split()
        if len(coord) == 4:
            coord = [float(x) for x in atom.split()[1:]]
            coordinates.append(coord)

    # read all energy profiles
    energy_profiles = read_energy_profiles(filename)

    reactant = pybel.readstring("xyz", xyz_str_reactant)
    product = pybel.readstring("xyz", xyz_str_product)

    mol = RWMol(MolFromSmiles(''))

    # create rdkit atoms based on openbabel reading of stringfile
    for i in range(num_atoms):
        a1 = reactant.atoms[i]
        symbol = openbabel.GetSymbol(a1.atomicnum)
        new_atom = Atom(symbol)
        new_atom.SetFormalCharge(a1.formalcharge)
        mol.AddAtom(new_atom)

    bmap1 = build_bond_map(reactant)
    bmap2 = build_bond_map(product)
    atom_core = set()

    # find core and iterate over bonds to add them to rdkit molecule
    for (src, tar), ob_bond in bmap1.items():
        if (src, tar) in bmap2:
            if bmap1[(src, tar)] != bmap2[(src, tar)]:
                    atom_core.add(src - 1)
                    atom_core.add(tar - 1)
        else:
            atom_core.add(src - 1)
            atom_core.add(tar - 1)
        mol.AddBond((src - 1), (tar - 1), ob_bond)
    for (src, tar), ob_bond in bmap2.items():
        if (src, tar) not in bmap1:
            atom_core.add(src - 1)
            atom_core.add(tar - 1)

    # set coordinates of rdkit molecule
    conf = Conformer(mol.GetNumAtoms())
    atom_id = 0
    for c in coordinates:
        atom_position = Point3D(c[0], c[1], c[2])
        conf.SetAtomPosition(atom_id, atom_position)
        atom_id += 1
    mol.AddConformer(conf)

    if visualize:
        visualize_rdkit_mol(mol, atom_core)
    return mol, atom_core, energy_profiles
