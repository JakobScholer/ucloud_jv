import unittest
from src.cut_molecule import find_all_cuts, make_cut_molecule
from rdkit.Chem import RWMol, MolFromSmiles

from rdkit.Chem.AllChem import Compute2DCoords
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolops import AddHs

import logging
logging.basicConfig(filename='test/logs/cut_molecule_test.log', level=logging.DEBUG)

class Test(unittest.TestCase):

    def test_make_molecule_1(self):
        smiles_string = "O=C=CCN(O)N(C(=O)O)N(CC=CO)N(O)C(O)C1CC1"
        rdk_mol = RWMol(MolFromSmiles(smiles_string))
        atom_core = {5,16,25,32}

        rdk_mol = AddHs(rdk_mol)
        Compute2DCoords(rdk_mol)

        # set size and id for atoms
        d = rdMolDraw2D.MolDraw2DCairo(500, 500) # or MolDraw2DCairo to get PNGs
        d.drawOptions().addAtomIndices = True
        rdMolDraw2D.PrepareAndDrawMolecule(d, rdk_mol, legend="derp")
        d.WriteDrawingText("derp.png")

        molecule, lookup = make_cut_molecule(rdk_mol, atom_core)

        derp = 0
        for m in molecule:
            logging.info(str(derp) + " node: " + str(m.id))
            logging.info("    children: " + str(m.children))
            derp += 1

        actual= {0}
        expected = {0}
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
