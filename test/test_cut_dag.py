import unittest
from src.root_mean_square import root_mean_square
from src.cut_molecule import cut_molecule_main, make_cut_molecule, find_all_cuts, make_cut
from src.generate_tree import generate_tree_main, reaction_and_product_to_gml, read_energy_profiles
from src.mod_to_xyz import mod_to_xyz_main, mod_to_xyz

class Test(unittest.TestCase):
    # check root
    def test_cut_dag_root(self):
        test = 'test/testfiles/stringfile_ring.xyz0000'
        graph = False
        dag = make_root(test, graph)

        actual_stringfile = dag.layers[0][0].stringfile
        actual_energy = dag.layers[0][0].energy
        Actual_cuts = dag.layers[0][0].cuts

    # check first childs
    def test_cut_dag_layer_1(self):
        test = 'test/testfiles/stringfile_ring.xyz0000'
        graph = False
        dag = make_root(test, graph)

    # check first second childs
    def test_cut_dag_layer_2(self):
        test = 'test/testfiles/stringfile_ring.xyz0000'
        graph = False
        dag = make_root(test, graph)

    # Lav børn på root og chek op på dem
    print("Second test")
    make_childs(dag.layers[0][0],dag)
    for node in dag.layers[1]:
        print(node.cuts)
        make_childs(node,dag)
    print("Second test")
    for node in dag.layers[2]:
        print(node.cuts)
    # lav børn på børnene :D det skal nok blive magisk


if __name__ == '__main__':
    unittest.main()
