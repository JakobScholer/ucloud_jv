from src.stringfile_tester import check_educt_to_product
from rdkit.Chem import RWMol, AddHs, MolFromSmiles, MolToXYZBlock
from os import walk, path, listdir

def make_reactions(smiles):# tag en smiles som input
# GØR HELE BLACK BOX DELEN!
    xyz_list = []
    for string in smiles:
        mol = RWMol(MolFromSmiles(string)) # lav rdkitmol fra smiles string
        mol = AddHs(mol) # add hydrogen for good measure.
        xyz_list.append(MolToXYZBlock(mol)) # convert til xyz fil

    reaction_name = smiles[0]
    for i in range(1,len(smiles)):
        reaction_name = reaction_name + "_+_" + string
    # kør blackbox
    #smiles_path = run_zstruct_and_gsm(xyz_list, smiles_string: reaction_name)
    smiles_path = "test_folder/" # black box wannabe tester
    stringfile_path = listdir(smiles_path)
    print(stringfile_path)
    reaction_folders = [smiles_path + s for s in stringfile_path]
    print(reaction_folders)
    for folder in reaction_folders: # gå over hver eneste stringfile+isomer og lav en cut dag
        containment = listdir(folder)
        if len(containment) > 1: # must be 2 files. if no stringfile was generated
            stringfile = folder + "/" + containment[0]
            isomer_file = folder + "/" + containment[1]
            print("succes!")
            if check_educt_to_product(stringfile): # if there is a reaction in the stringfile. make a cut dag!
                print("with a top!")
                #generate_cut_dag_main(stringfile, isomer_file, folder)
