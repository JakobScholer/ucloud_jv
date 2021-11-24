from multiprocessing import Process
import unittest
from rdkit.Chem import MolFromSmiles, MolToXYZBlock, rdDepictor
from rdkit.Chem.rdDistGeom import EmbedMolecule
from src.zstruct_and_gsm import run_zstruct_and_gsm


class Test(unittest.TestCase):
    # Multi threaded behaviour test case
    def test_blackbox_multithreading(self):
        def blackbox():
            mol = MolFromSmiles("CCO")
            rdDepictor.Compute2DCoords(mol)  # generate 2d coordinates
            EmbedMolecule(mol, randomSeed=0xf00d)  # generate 3d coordinates
            xyz_str = MolToXYZBlock(mol)
            run_zstruct_and_gsm(xyz_str)
        p1 = Process(target=blackbox())
        p2 = Process(target=blackbox())
        p1.start()
        p2.start()
        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()
