import unittest
from src.mod_to_xyz import mod_to_xyz
from mod import smiles


class Test(unittest.TestCase):
    # Normal behaviour test case
    def test_root_mean_square_normal(self):

        g = smiles("CCO")
        x = mod_to_xyz(g, toFile=False)
        actual = x
        expected = 2.9530847717119375
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
