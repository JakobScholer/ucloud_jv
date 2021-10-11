import sys
from mod import smiles, graphGMLString

from src.root_mean_square import root_mean_square
from src.runner import runner_main, make_cut_molecule, find_all_cuts, make_cut
from src.generate_tree import generate_tree_main, reaction_and_product_to_gml, read_energy_profiles
from src.mod_to_xyz import mod_to_xyz_main, mod_to_xyz

if __name__ == '__main__':
    if len(sys.argv) == 2:
        if str(sys.argv[1]) == "runner":
            runner_main()
        elif str(sys.argv[1]) == "generate_tree":
            generate_tree_main()
        elif str(sys.argv[1]) == "mod_to_xyz":
            mod_to_xyz_main()
        elif str(sys.argv[1]) == "xyz_to_mod":
            runner_main()
    else:
        g = smiles("CCO")                       # molecule to test reaction on
        mod_to_xyz(g, to_file=True)             # convert molecule for zstruct to understand it

        # run zstruct with molecule.xyz (molecule.frozen is empty for now)

        gml_string, atom_core, energy_curve = reaction_and_product_to_gml('src/stringfile.xyz0000', visualize=True)
        g = graphGMLString(gml_string)
        molecule, lookup_dict = make_cut_molecule(g, atom_core)
        cuts = set()
        find_all_cuts(molecule, cuts, lookup_dict, 0)

        # iterate over all cuts and do:
        new_gml_string = make_cut(g, cuts, molecule)
        new_g = graphGMLString(gml_string)
        mod_to_xyz(new_g, to_file=True)
        # run zstruct with new molecule.xyz (molecule.frozen now contains all atoms that are not in the core)
        with open('src/stringfile.xyz0000') as fi:
            ct = fi.readlines()
        curve = read_energy_profiles(ct)
        x = root_mean_square(energy_curve, curve)   # based on value decide if curve is the same

