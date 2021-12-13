import openbabel.pybel as pybel
from openbabel import openbabel
import plotly.graph_objects as go
from rdkit.Chem import RWMol, MolFromSmiles, Atom, Conformer
from rdkit.Geometry import Point3D
from src.bond_map import build_bond_map


def fig_plot(mol, core_atoms):
    """Takes an rdkit mol object and the core atoms of the molecule, displays a visual representation of the molecule"""
    num_atoms = mol.GetNumAtoms()   # number of atoms
    num_bonds = mol.GetNumBonds()   # number of bonds

    atom_labels = [mol.GetAtoms()[k].GetSymbol() for k in range(num_atoms)]     # list of all atom types
    bond_labels = [mol.GetBonds()[k].GetBondType() for k in range(num_bonds)]   # list of all bond types

    # find all core bonds
    core_bonds = []
    edge_counter = 0
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetIdx() in core_atoms and bond.GetEndAtom().GetIdx() in core_atoms:
            core_bonds.append(edge_counter)
        edge_counter += 1

    bond_x_start = [mol.GetConformer(0).GetAtomPosition(mol.GetBonds()[k].GetBeginAtom().GetIdx())[0] for k in range(num_bonds)]
    bond_y_start = [mol.GetConformer(0).GetAtomPosition(mol.GetBonds()[k].GetBeginAtom().GetIdx())[1] for k in range(num_bonds)]
    bond_z_start = [mol.GetConformer(0).GetAtomPosition(mol.GetBonds()[k].GetBeginAtom().GetIdx())[2] for k in range(num_bonds)]
    bond_x_end = [mol.GetConformer(0).GetAtomPosition(mol.GetBonds()[k].GetEndAtom().GetIdx())[0] for k in range(num_bonds)]
    bond_y_end = [mol.GetConformer(0).GetAtomPosition(mol.GetBonds()[k].GetEndAtom().GetIdx())[1] for k in range(num_bonds)]
    bond_z_end = [mol.GetConformer(0).GetAtomPosition(mol.GetBonds()[k].GetEndAtom().GetIdx())[1] for k in range(num_bonds)]
    bond_x = []
    bond_y = []
    bond_middle_x = []  # x coordinate of middle of line
    bond_middle_y = []  # y coordinate of middle of line
    for i in range(num_bonds):
        bond_x += [bond_x_start[i], bond_x_end[i], None]
        bond_y += [bond_y_start[i], bond_y_end[i], None]
        bond_middle_x.append((float(bond_x_start[i]) + float(bond_x_end[i])) / 2)
        bond_middle_y.append((float(bond_y_start[i]) + float(bond_y_end[i])) / 2)

    bond_core_x = []
    bond_core_y = []
    for bond in core_bonds:
        bond_core_x.append(bond_middle_x[bond])
        bond_core_y.append(bond_middle_y[bond])

    atom_x = [mol.GetConformer(0).GetAtomPosition(k)[0] for k in range(num_atoms)]
    atom_y = [mol.GetConformer(0).GetAtomPosition(k)[1] for k in range(num_atoms)]
    atom_z = [mol.GetConformer(0).GetAtomPosition(k)[2] for k in range(num_atoms)]

    atom_core_x = []
    atom_core_y = []
    for atom in core_atoms:
        atom_core_x.append(atom_x[atom])
        atom_core_y.append(atom_y[atom])

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
                             name='bondIDs',
                             text=list(range(1, num_bonds+1)),
                             hoverinfo='text',
                             textfont_size=1
                             ))
    fig.add_trace(go.Scatter(x=bond_middle_x,
                             y=bond_middle_y,
                             mode='text',
                             name='bondtypes',
                             text=bond_labels,
                             hoverinfo='skip',
                             textfont_size=10
                             ))
    fig.add_trace(go.Scatter(x=bond_x,
                             y=bond_y,
                             mode='lines',
                             name='bonds',
                             line=dict(color='rgb(210,210,210)', width=1),
                             text='',
                             hoverinfo='skip'
                             ))

    fig.add_trace(go.Scatter(x=atom_x,
                             y=atom_y,
                             mode='text',
                             name='atomIDs',
                             text=list(range(1, num_atoms+1)),
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
                             text=atom_labels,
                             hoverinfo='skip',
                             ))
    fig.show()


def read_energy_profiles(string_content):
    """takes the content of a stringfile, returns a list of all energy values within it."""
    energy_profiles = []
    for i in range(1, len(string_content), (int(string_content[0]) + 2)):
        energy_profiles.append(float(string_content[i]))
    return energy_profiles


def stringfile_to_rdkit(filename, visualize=False):
    """takes a stringfile, returns a rdkit mol object, the core of the atom and the energy profile."""
    with open(filename) as f:
        content = f.readlines()
    num_atoms = int(content[0])
    xyz_str_reactant: str = "".join(content[:(num_atoms + 2)])                  # string representing reactant
    xyz_str_product: str = "".join(content[len(content) - (num_atoms + 2):])    # string representing product

    # find all atom coordinates and store in list
    atoms = xyz_str_product.split("\n")
    coordinates = []
    for atom in atoms:
        coord = atom.split()
        if len(coord) == 4:
            coord = [float(x) for x in atom.split()[1:]]
            coordinates.append(coord)

    # read all energy profiles
    energy_profiles = read_energy_profiles(content)

    reactant = pybel.readstring("xyz", xyz_str_reactant)
    product = pybel.readstring("xyz", xyz_str_product)

    num_atoms: int = len(reactant.atoms)
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
        fig_plot(mol, atom_core)
    return mol, atom_core, energy_profiles


def stringfile_to_rdkit_main():
    mol, core, energy = stringfile_to_rdkit("xyz_test_files/GCD_test_files/stringfile.xyz0025", True)
    print(core)
    print(energy)
