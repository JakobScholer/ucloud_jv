from src.stringfile_tester import check_educt_to_product
from rdkit.Chem import RWMol, AddHs, MolFromSmiles, MolToXYZBlock, rdDepictor
from rdkit.Chem.AllChem import EmbedMolecule
from os import listdir
from src.blackbox import run_zstruct_and_gsm
from src.generate_cut_dag import make_cut_dag_2

# C=C(C)C(C(CC)CN(C(=O)OC(C)=O)C([O-])=NC(C)C(C=CC)C1CCCCC1)C2CCCCC2



# takes a mode and a string or list of strings as input
    # mode 0 runs blackbox and the list of strings must be the smiles for reactions
    # mode 1 reads data from a folder. string_data must be the path to the already compiled cut dag data
def make_reactions(mode: int, string_data, visual_cut_dag: bool=False, visual_stringfiles: bool=False):
# GØR HELE BLACK BOX DELEN!
    if mode == 0:
        xyz_list = []
        for string in string_data:
            mol = RWMol(MolFromSmiles(string)) # lav rdkitmol fra smiles string
            mol = AddHs(mol) # add hydrogen for good measure.
            rdDepictor.Compute2DCoords(mol) # add coordinates with a comformer
            EmbedMolecule(mol, randomSeed=0xf00d)
            xyz_list.append(MolToXYZBlock(mol)) # convert til xyz fil
        #print(xyz_list[0])
        reaction_name = string_data[0]
        for i in range(1,len(string_data)):
            reaction_name = reaction_name + "_+_" + string
        # kør blackbox
        smiles_path = run_zstruct_and_gsm(xyz_list, reaction_name)
    elif mode == 1:
        smiles_path = string_data

    stringfile_path = listdir(smiles_path)#
    reaction_folders = [smiles_path + "/" + s for s in stringfile_path]

    for folder in reaction_folders: # gå over hver eneste stringfile+isomer og lav en cut dag
        containment = listdir(folder)
        if len(containment) > 2: # must be 2 files. if no stringfile was generated
            for file in containment:
                if "ISOMER" in file:
                    isomer_file = folder + "/" + str(file)
                elif "stringfile" in file:
                    stringfile = folder + "/" + str(file)
            print("    stringfile: " + str(stringfile))
            #print("    ISOMER: " + str(isomer_file))
            if check_educt_to_product(stringfile): # if there is a reaction in the stringfile. make a cut dag!
                print("        Generate cut dag")
                make_cut_dag_2(mode, stringfile, visual_cut_dag, visual_stringfiles)
            else:
                print("        Not generating cut dag")
