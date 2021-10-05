import unittest
from src.mod_to_xyz import mod_to_xyz
from mod import smiles


class Test(unittest.TestCase):
    # Normal behaviour test case
    def test_mod_to_xyz_normal(self):
        g = smiles("CCO")
        mol = mod_to_xyz(g, toFile=False)
        actual = mol.GetNumAtoms()
        expected = 9
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
