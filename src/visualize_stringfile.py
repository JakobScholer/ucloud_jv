import openbabel.pybel as pybel
from openbabel import openbabel
from src.stringfile_helper_functions import build_bond_map
from rdkit.Chem import RWMol, MolFromSmiles, Atom
from rdkit.Chem.AllChem import Compute2DCoords
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolops import RemoveHs

import numpy as np
import PIL
from PIL import Image
from os import remove

def read_stringfile(strfile): # read a stringfile and return a list with energy and openbabel mol, for each step in the reaction
    # read xyz data as string
    with open(strfile) as f:
        content = f.readlines()
    num_atoms = int(content[0])

    reaction_list = []
    pointer = 0
    while pointer < len(content): # seperate all reaction steps
        reaction_step: str = "".join(content[pointer:(pointer + num_atoms + 2)]) # get data of a single step
        molecule_ob = pybel.readstring("xyz", reaction_step)
        reaction_list.append((content[pointer+1].strip("\n"),molecule_ob))
        pointer += (num_atoms + 2)  # increase pointer to next reaction step
    return reaction_list

def openbabel_to_rdkit(ob_mol): # take a openbabel mol and transform it to a RDkit mol
    num_atoms: int = len(ob_mol.atoms)
    mol = RWMol(MolFromSmiles(''))
    # create rdkit atoms based on openbabel reading of stringfile
    for i in range(num_atoms):
        a1 = ob_mol.atoms[i]
        symbol = openbabel.GetSymbol(a1.atomicnum)
        new_atom = Atom(symbol)
        new_atom.SetFormalCharge(a1.formalcharge)
        mol.AddAtom(new_atom)

    bmap = build_bond_map(ob_mol)

    # find core and iterate over bonds to add them to rdkit molecule
    for (src, tar), ob_bond in bmap.items():
        mol.AddBond((src - 1), (tar - 1), ob_bond)

    return mol # return rdkit mol


def find_core(ob_educt, ob_product): # from to openbabel mols, get the core atoms of the whole reaction
    bmap1 = build_bond_map(ob_educt)
    bmap2 = build_bond_map(ob_product)
    atom_core = set()
    bond_core = set()

    # find core and iterate over bonds to add them to rdkit molecule
    for (src, tar), ob_bond in bmap1.items():
        if (src, tar) in bmap2:
            if bmap1[(src, tar)] != bmap2[(src, tar)]:
                    atom_core.add(src - 1)
                    atom_core.add(tar - 1)
                    bond_core.add((src - 1, tar - 1))
        else:
            atom_core.add(src - 1)
            atom_core.add(tar - 1)
            bond_core.add((src - 1, tar - 1))
    for (src, tar), ob_bond in bmap2.items():
        if (src, tar) not in bmap1:
            atom_core.add(src - 1)
            atom_core.add(tar - 1)
            bond_core.add((src - 1, tar - 1))

    return atom_core, bond_core

def make_mol_png(mol, core_atoms, core_bond_pairs, hydrogens: bool, png_name, titel, size: int=500):

    # prep mol with hydrofens and coords acoordingly
    if hydrogens:
        mol = RemoveHs(mol)
    else:
        mol.UpdatePropertyCache()
        mol = RemoveHs(mol, True)
        # find all bonds in core that is stil relevant for the molecule
        core_bonds = [mol.GetBondBetweenAtoms(pair[0],pair[1]).GetIdx() for pair in core_bond_pairs if not mol.GetBondBetweenAtoms(pair[0],pair[1]) == None]
    Compute2DCoords(mol)

    # set size and id for atoms
    d = rdMolDraw2D.MolDraw2DCairo(size, size) # or MolDraw2DCairo to get PNGs
    d.drawOptions().addAtomIndices = True

    # draw and safe image of molecule
    if hydrogens:
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol, legend=titel)
    else:
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=core_atoms, highlightBonds=core_bonds, legend=titel)
    d.WriteDrawingText(png_name)

    return png_name

def combine_images(right_images, left_images, image_name):
    # get image data
    right_pil_imgs    = [ PIL.Image.open(i) for i in right_images]
    left_pil_imgs    = [ PIL.Image.open(i) for i in left_images]

    # prepare for horisont image
    right_imgs_comb = np.hstack(right_pil_imgs)
    left_imgs_comb = np.hstack(left_pil_imgs)

    right_imgs_comb = PIL.Image.fromarray( right_imgs_comb)
    right_imgs_comb.save( 'top.png' )
    left_imgs_comb = PIL.Image.fromarray( left_imgs_comb)
    left_imgs_comb.save( 'bot.png' )

    # for a vertical stacking it is simple: use vstack
    final_imgs_comb = np.vstack( [ PIL.Image.open(i) for i in ['top.png','bot.png']] )
    final_imgs_comb = PIL.Image.fromarray( final_imgs_comb)
    final_imgs_comb.save(image_name)

    # clean up after wards
    remove('bot.png')
    remove('top.png')


def visualize_2D(stringfile_path: str, image_path: str, image_name: str = "Reaction_scheme.jpg"):
    stringfile_data = read_stringfile(stringfile_path) # get list of energi levels and openbabel mols   ####### GetBondBetweenAtoms(0,1) ######

    core_atoms, core_bond_pairs = find_core(stringfile_data[0][1], stringfile_data[len(stringfile_data)-1][1]) # find core info

    left_image = "l.png"
    right_image = "r.png"
    reaction_step = 1
    left_image_list = []
    right_image_list = []
    # make left and right side image for each mol in the reaction
    for pair in stringfile_data:
        mol = openbabel_to_rdkit(pair[1]) # transform openbabel mol to rdkit mol
        left_image_list.append(make_mol_png(mol, core_atoms, core_bond_pairs, False, image_path + "/" + left_image, "Energy level: " + pair[0])) # without hydrogens
        right_image_list.append(make_mol_png(mol, core_atoms, core_bond_pairs, True, image_path + "/" + right_image, "Reaction step: " + str(reaction_step))) # with hydrogens
        left_image = "l" + left_image
        right_image = "r" + right_image
        reaction_step += 1

    # combine all images
    combine_images(right_image_list, left_image_list, image_path + "/" + image_name)

    # clean all images
    for file in (left_image_list + right_image_list):
        remove(file)
