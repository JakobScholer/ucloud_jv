from src.blackbox2 import run_zstruct, run_gsm_initial
from src.stringfile_tester import check_educt_to_product
from rdkit.Chem import RWMol, AddHs, MolFromSmiles, MolToXYZBlock, rdDepictor
from rdkit.Chem.AllChem import EmbedMolecule
from os import listdir
from glob import glob
from src.generate_cut_dag import make_cut_dag
from src.stringfile_helper_functions import max_energy_curve


# takes a mode and a string or list of strings as input
    # blackbox True runs blackbox and the list of strings must be the smiles for reactions
    # blackbox False reads data from a folder. string_data must be the path to the already compiled cut dag data
def make_reactions(blackbox: bool, string_data, max_energy: int=50, generate_initial_stringfiles: bool=True, visual_cut_dag: bool=False, visual_stringfiles: bool=False, debug: bool=False):
    if blackbox and generate_initial_stringfiles:
        if isinstance(string_data, list):
            xyz_list = []

            for string in string_data:
                mol = RWMol(MolFromSmiles(string)) # lav rdkitmol fra smiles string
                mol = AddHs(mol, explicitOnly=False) # add hydrogen for good measure.
                rdDepictor.Compute2DCoords(mol) # add coordinates with a comformer
                EmbedMolecule(mol, randomSeed=0xf00d)
                xyz_list.append(MolToXYZBlock(mol)) # convert til xyz fil

            reaction_name = string_data[0]
            for i in range(1,len(string_data)):
                reaction_name = reaction_name + "_+_" + string_data[i]
            # kør blackbox
            #educt, _, _ = read_stringfile_content("blackbox/output/string_radicalO2_recombination.xyz")
            #xyz_list.append(educt)
            smiles_path, isomer_count = run_zstruct(reaction_name, xyz_list, debug=False)
            run_gsm_initial(smiles_path, isomer_count)
        else:
            smiles_path = f"blackbox/output/{string_data}"
            isomer_count = len(listdir(smiles_path))
            # kør blackbox
            run_gsm_initial(smiles_path, isomer_count)
    else:
        smiles_path = string_data

    stringfile_path = listdir(smiles_path)
    reaction_folders = [smiles_path + "/" + s for s in stringfile_path]
    reaction_folders.sort()
    print("-------------------------------------------------------------------")
    # make cut dag for every reaction folder containing a stringfile
    for folder in reaction_folders:
        stringfiles = glob(f"{folder}/stringfile*")
        if stringfiles: # check if stringfile exists
            print(f"working on {folder}")
            log_data = ""

            check = check_educt_to_product(stringfiles[0])
            max = max_energy_curve(stringfiles[0], max_energy)
            if check_educt_to_product(stringfiles[0]) and max_energy_curve(stringfiles[0], max_energy): # if there is a reaction in the stringfile. make a cut dag!
                log_data += "Cut dag generated\n"
                make_cut_dag(blackbox, stringfiles[0], visual_cut_dag, visual_stringfiles, debug)
            else:
                log_data += "Cut dag not generated\n"
                log_data += f"Educt to product: {check}\n"
                log_data += f"Max energy curve: {max}\n"
            with open(f"{folder}/done.txt", "w") as f:
                f.write(log_data)

