from src.stringfile_tester import check_educt_to_product
from rdkit.Chem import RWMol, AddHs, MolFromSmiles, MolToXYZBlock, rdDepictor
from rdkit.Chem.AllChem import EmbedMolecule
from os import walk, path, listdir
from src.zstruct_and_gsm import run_zstruct_and_gsm
from src.generate_cut_dag import generate_cut_dag_main

# C=C(C)C(C(CC)CN(C(=O)OC(C)=O)C([O-])=NC(C)C(C=CC)C1CCCCC1)C2CCCCC2

def make_reactions(smiles):# tag en smiles som input
# GØR HELE BLACK BOX DELEN!
    xyz_list = []
    for string in smiles:
        mol = RWMol(MolFromSmiles(string)) # lav rdkitmol fra smiles string
        mol = AddHs(mol) # add hydrogen for good measure.
        rdDepictor.Compute2DCoords(mol) # add coordinates with a comformer
        EmbedMolecule(mol, randomSeed=0xf00d)
        xyz_list.append(MolToXYZBlock(mol)) # convert til xyz fil

    reaction_name = smiles[0]
    for i in range(1,len(smiles)):
        reaction_name = reaction_name + "_+_" + string
    # kør blackbox
    smiles_path = run_zstruct_and_gsm(xyz_list, reaction_name)
    #smiles_path = "test_folder/" # black box wannabe tester
    stringfile_path = listdir(smiles_path)
    reaction_folders = [smiles_path + s for s in stringfile_path]

    for folder in reaction_folders: # gå over hver eneste stringfile+isomer og lav en cut dag
        print(folder)
        containment = listdir(folder)
        if len(containment) > 2: # must be 2 files. if no stringfile was generated
            for file in containment:
                if "ISOMER" in file:
                    isomer_file = folder + "/" + str(file)
                elif "stringfile" in file:
                    stringfile = folder + "/" + str(file)
            print("    stringfile: " + str(stringfile))
            print("    ISOMER: " + str(isomer_file))
            if check_educt_to_product(stringfile): # if there is a reaction in the stringfile. make a cut dag!
                print("        Generate cut dag")
                generate_cut_dag_main(stringfile, isomer_file)
