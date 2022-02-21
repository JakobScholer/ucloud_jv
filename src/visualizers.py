from os import listdir
from os.path import isdir, isfile
import plotly.graph_objects as go
from src.stringfile_helper_functions import read_energy_profiles
from src.stringfile_tester import check_educt_to_product

def visualize_rdkit_mol(mol, core_atoms: set = None):
    """Takes an rdkit mol object and the core atoms of the molecule, displays a visual representation of the molecule"""
    if core_atoms is None:
        core_atoms = set()
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
                             hoverinfo='skip'
                             ))
    fig.add_trace(go.Scatter(x=bond_middle_x,
                             y=bond_middle_y,
                             mode='text',
                             name='bondIDs',
                             text=list(range(0, num_bonds)),
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
                             text=list(range(0, num_atoms)),
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


def visualize_cut_dag(cut_dag):
    """Takes a cut dag object, displays a visual representation of the cut dag"""
    cut_option_y_green = []
    cut_option_x_green = []
    cut_option_y_red = []
    cut_option_x_red = []
    cut_option_y_black = []
    cut_option_x_black = []
    cut_option_y = []
    cut_option_x = []
    cut_info = []
    stringfiles = []
    bond_y = []
    bond_x = []
    for layer in cut_dag.layers.keys():
        layer_length = len(cut_dag.layers.get(layer))
        for position in range(layer_length):
            if cut_dag.layers.get(layer)[position].stringfile == "NO REACTION":
                cut_option_y_black.append(layer * 10)
                cut_option_x_black.append(position * 10 - (layer_length * 10) / 2)
            elif cut_dag.layers.get(layer)[position].RMS == 1:
                cut_option_y_green.append(layer * 10)
                cut_option_x_green.append(position * 10 - (layer_length * 10) / 2)
            elif cut_dag.layers.get(layer)[position].RMS == 0:
                cut_option_y_red.append(layer * 10)
                cut_option_x_red.append(position * 10 - (layer_length * 10) / 2)
            cut_option_y.append(layer * 10)
            cut_option_x.append(position * 10 - (layer_length * 10) / 2)
            cut_info.append(cut_dag.layers.get(layer)[position].cuts)
            stringfiles.append(cut_dag.layers.get(layer)[position].stringfile)
            for child_position in cut_dag.layers.get(layer)[position].childs:
                bond_y += [layer * 10, (layer+1) * 10, None]
                bond_x += [position * 10 - (layer_length * 10)/2, child_position * 10 - (len(cut_dag.layers.get(layer+1)) * 10)/2, None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=bond_x,
                             y=bond_y,
                             mode='lines',
                             name='connections',
                             line=dict(color='rgb(210,210,210)', width=3),
                             hoverinfo='skip'
                             ))

    fig.add_trace(go.Scatter(x=cut_option_x,
                             y=cut_option_y,
                             mode='markers',
                             name='Cuts',
                             marker=dict(symbol='circle-dot',
                                         size=20,
                                         color='#800080'
                                         ),
                             text=cut_info,
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x_red,
                             y=cut_option_y_red,
                             mode='markers',
                             name='Cuts',
                             marker=dict(symbol='circle-dot',
                                         size=20,
                                         color='#FF0000'
                                         ),
                             text=cut_info,
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x_green,
                             y=cut_option_y_green,
                             mode='markers',
                             name='Cuts',
                             marker=dict(symbol='circle-dot',
                                         size=20,
                                         color='#008000'
                                         ),
                             text=cut_info,
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x_black,
                             y=cut_option_y_black,
                             mode='markers',
                             name='Cuts',
                             marker=dict(symbol='circle-dot',
                                         size=20,
                                         color='#000000'
                                         ),
                             text=cut_info,
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x,
                             y=cut_option_y,
                             mode='text',
                             name='Cut info',
                             text=cut_info,
                             hoverinfo='skip',
                             textfont_size=10
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x,
                             y=[k - 1 for k in cut_option_y],
                             mode='text',
                             name='stringfile',
                             text=stringfiles,
                             hoverinfo='skip',
                             textfont_size=10
                             ))

    fig.show()


def visualize_energy_curves(folder: str):
    """Takes a reaction folder as input, displays a visual representation of the energy curves for the stringfiles with different cuts"""
    energy_profiles = []
    energy_profiles_cuts = []
    original_energy_profile = []
    for file in listdir(folder):
        if isdir(f"{folder}/{file}"):
            if isfile(f"{folder}/{file}/stringfile.xyz0000"):
                energy_profile = read_energy_profiles(f"{folder}/{file}/stringfile.xyz0000")
                energy_profiles.append(energy_profile)
                energy_profiles_cuts.append(file)
        elif file.startswith("stringfile"):
            original_energy_profile = read_energy_profiles(f"{folder}/{file}")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(range(0, len(original_energy_profile))), y=original_energy_profile,
                             mode='lines',
                             name='original'))
    for ep, epc in zip(energy_profiles, energy_profiles_cuts):
        fig.add_trace(go.Scatter(x=list(range(0, len(ep))), y=ep,
                                 mode='lines',
                                 name=epc))
    fig.show()


def energy_curve_all_reactions(folder: str, max_energy: int = 100):
    energy_profiles = []
    energy_profiles_names = []

    reaction_folders = listdir(folder)
    reaction_folders.sort()
    for react_folder in reaction_folders:
        if isdir(f"{folder}/{react_folder}"): # take only folders
            from glob import glob
            stringfiles = glob(f"{folder}/{react_folder}/stringfile*")
            energy_profile = read_energy_profiles(stringfiles[0])
            if check_educt_to_product(stringfiles[0]) and max(energy_profile) < max_energy: # test the stringfile if a reaction happens
                energy_profiles.append(energy_profile)
                energy_profiles_names.append(react_folder)
    fig = go.Figure()
    for ep, epc in zip(energy_profiles, energy_profiles_names):
        fig.add_trace(go.Scatter(x=list(range(0, len(ep))), y=ep,
                                 mode='lines',
                                 name=epc))
    fig.show()
