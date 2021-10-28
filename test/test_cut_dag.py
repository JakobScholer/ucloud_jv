import unittest
from src.cut_molecule import cut_molecule_main, make_cut_molecule, find_all_cuts, make_cut
from src.generate_tree import generate_tree_main, reaction_and_product_to_gml, read_energy_profiles
from src.mod_to_xyz import mod_to_xyz_main, mod_to_xyz
from src.cut_dag import make_childs, make_root

class Test(unittest.TestCase):
    # check root
    def test_cut_dag_root(self):
        test = 'test/testfiles/stringfile_ring.xyz0000'
        graph = False
        dag = make_root(test, graph)

        actual_stringfile = dag.layers[0][0].stringfile
        actual_energy = dag.layers[0][0].energy
        actual_cuts = dag.layers[0][0].cuts

        expected_stringfile = "test/testfiles/stringfile_ring.xyz0000"
        expected_energy = [0.0, 7.441711, 8.316389, -19.231187, -17.474338, 5.249445, -18.4045, -45.613282, -64.595534]
        expected_cuts = set()

        self.assertEqual(actual_stringfile, expected_stringfile)
        self.assertEqual(actual_energy, expected_energy)
        self.assertEqual(actual_cuts, expected_cuts)

    # check first childs
    def test_cut_dag_layer_1(self):
        test = 'test/testfiles/stringfile_ring.xyz0000'
        graph = False
        dag = make_root(test, graph)

        make_childs(dag.layers[0][0],dag)
        actual = []
        expected = [{1},{3},{9}]
        for node in dag.layers[1]:
            actual.append(node.cuts)

        for i in range(len(expected)):
            self.assertEqual(actual[i],expected[i])

    # check first second childs
    def test_cut_dag_layer_2(self):
        test = 'test/testfiles/stringfile_ring.xyz0000'
        graph = False
        dag = make_root(test, graph)

        actual = []
        expected = [{1,9},{1,3},{3,9}]
        make_childs(dag.layers[0][0],dag)
        for node in dag.layers[1]:
            make_childs(node,dag)
        for node in dag.layers[2]:
            actual.append(node.cuts)

        for node in dag.layers[2]:
            actual.append(node.cuts)

        for i in range(len(expected)):
            self.assertEqual(actual[i],expected[i])

if __name__ == '__main__':
    unittest.main()
