from stringfile_tester import check_educt_to_product
from rdkit.Chem import RWMol, AddHs, MolFromSmiles, MolToXYZBlock

def make_reactions(smiles):# tag en smiles som input
# GØR HELE BLACK BOX DELEN!
    xyz_list = []
    for string in smiles:
        mol = RWMol(MolFromSmiles(string)) # lav rdkitmol fra smiles string
        mol = AddHs(mol) # add hydrogen for good measure.
        xyz_list.append(MolToXYZBlock(mol)) # convert til xyz fil

    reaction_name = smiles[0]
    for i in range(1,len(smiles))
        reaction_name = reaction_name + "_+_" + string
    #stringfile_list, isomer_list = blackbox(xyz_string) MAGIC CODE GOES HERE! ### # kør blackbox
    path = run_zstruct_and_gsm(xyz_list, smiles_string: reaction_name)
# gå over hver eneste stringfile+isomer og lav en cut dag
    # hvis der ikke sker en reaktin (educt != product)
        # skip
    # else
        # lav cut dag for den data
        generate_cut_dag_main(stringfile, isomer, folder)
