from stringfile_tester import check_educt_to_product
from rdkit.Chem import RWMol, AddHs, MolFromSmiles, MolToXYZBlock

def make_reactions(smiles):# tag en smiles som input
# GØR HELE BLACK BOX DELEN!
    mol = RWMol(MolFromSmiles(smiles)) # lav rdkitmol fra smiles string
    mol = AddHs(mol) # add hydrogen for good measure.
    xyz_string = MolToXYZBlock(mol) # convert til xyz fil
    #stringfile_list, isomer_list = blackbox(xyz_string) MAGIC CODE GOES HERE! ### # kør blackbox
# gå over hver eneste stringfile+isomer og lav en cut dag
    # hvis der ikke sker en reaktin (educt != product)
        # skip
    # else
        # lav cut dag for den data
