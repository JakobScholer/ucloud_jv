from os import chdir, system, listdir
import sys
from mod import smiles, graphGMLString

from src.root_mean_square import root_mean_square
from src.cut_molecule import cut_molecule_main, make_cut_molecule, find_all_cuts, make_cut
from src.generate_tree import generate_tree_main, reaction_and_product_to_gml, read_energy_profiles
from src.mod_to_xyz import mod_to_xyz_main, mod_to_xyz
from src.cut_dag import cut_dag_main, make_root
from src.generate_cut_dag import generate_cut_dag_main
from src.zstruct_and_gsm import generate_isomers, run_ssm


if __name__ == '__main__':
    if len(sys.argv) == 2:
        if str(sys.argv[1]) == "cut_molecule":
            cut_molecule_main()
        elif str(sys.argv[1]) == "generate_tree":
            generate_tree_main()
        elif str(sys.argv[1]) == "mod_to_xyz":
            mod_to_xyz_main()
        elif str(sys.argv[1]) == "xyz_to_mod":
            cut_molecule_main()
        elif str(sys.argv[1]) == "cut_dag":
            cut_dag_main()
        elif str(sys.argv[1]) == "generate_cut_dag":
            generate_cut_dag_main()
    else:
        '''
        g = smiles("CCO")                       # molecule to test reaction on
        mod_to_xyz(g, to_file=True)             # convert molecule for zstruct to understand it
        system(f"mv mod_coordinates.xyz blackbox/data")

        # run zstruct with molecule.xyz (molecule.frozen is empty for now)
        chdir("blackbox")
        isomer_count = generate_isomers("data/main_example.json")
        run_ssm(isomer_count)
        chdir("..")
        '''
        for file in listdir("blackbox/scratch/stringfiles"):
            if not file.startswith("stringfile"):
               continue
            print(file)
            tree = make_root(f'blackbox/scratch/stringfiles/{file}', False)
            print(tree)
            '''
            for cut in all_cuts:
                new_gml_string, order = make_cut(g, cuts, molecule)
                new_g = graphGMLString(gml_string)
                mod_to_xyz(new_g, to_file=True)
                # new instance of zstruct xtb with new folders
            '''

        '''
        gml_string, atom_core, energy_curve = reaction_and_product_to_gml('blackbox/scratch/stringfiles/stringfile.xyz0010', visualize=True)
        g = graphGMLString(gml_string)
        molecule, lookup_dict = make_cut_molecule(g, atom_core)
        cuts = set()
        find_all_cuts(molecule, cuts, lookup_dict, 0)

        # iterate over all cuts and do:
        new_gml_string, order = make_cut(g, cuts, molecule)
        new_g = graphGMLString(gml_string)
        mod_to_xyz(new_g, to_file=True)
        # run zstruct with new molecule.xyz (molecule.frozen now contains all atoms that are not in the core)
        with open('blackbox/scratch/stringfiles/stringfile.xyz0010') as fi:
            ct = fi.readlines()
        curve = read_energy_profiles(ct)
        x = root_mean_square(energy_curve, curve)   # based on value decide if curve is the same
        '''