import unittest
from mod import graphGMLString
from src.runner import find_all_cuts, make_cut_molecule
from src.generate_tree import reaction_and_product_to_gml


class Test(unittest.TestCase):
    # Normal tree
    def test_find_cuts_tree(self):
        gml, ac, ec = reaction_and_product_to_gml('test/testfiles/stringfile_tree.xyz0000', visualize=False)
        g = graphGMLString(gml)
        m, lookup = make_cut_molecule(g, ac)
        actual = set()
        find_all_cuts(m, actual, lookup, 0)
        expected = {0}
        self.assertEqual(actual, expected)

    # ring INFINITE CURRENTLY
    def test_find_cuts_ring(self):
        gml, ac, ec = reaction_and_product_to_gml('test/testfiles/stringfile_ring.xyz0000', visualize=False)
        g = graphGMLString(gml)
        m, lookup = make_cut_molecule(g, ac)
        actual = set()
        find_all_cuts(m, actual, lookup, 0)
        expected = {1, 3, 9}
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
