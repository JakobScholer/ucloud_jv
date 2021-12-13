import sys
from src.cut_molecule import cut_molecule_main
from src.stringfile_to_gml import stringfile_to_gml
from src.stringfile_to_rdkit import stringfile_to_rdkit_main, stringfile_to_rdkit
from src.cut_dag import cut_dag_main
from src.generate_cut_dag import generate_cut_dag_main
from src.zstruct_and_gsm import zstruct_gsm_main
from src.smiles_to_reactions import make_reactions
from src.stringfile_tester import stringfile_tester_main

if __name__ == '__main__':
    if len(sys.argv) == 2:
        if str(sys.argv[1]) == "cut_molecule":
            cut_molecule_main()
        elif str(sys.argv[1]) == "generate_tree":
            stringfile_to_rdkit_main()
        elif str(sys.argv[1]) == "cut_dag":
            cut_dag_main()
        elif str(sys.argv[1]) == "generate_cut_dag":
            generate_cut_dag_main()
        elif str(sys.argv[1]) == "visualize_molecule":
            stringfile_to_rdkit("xyz_test_files/reaction0001/stringfile.xyz0001", visualize=True)
        elif str(sys.argv[1]) == "zstruct_gsm":
            zstruct_gsm_main()
        elif str(sys.argv[1]) == "smiles_to_reactions":
            make_reactions(["CCCCO"])
        elif str(sys.argv[1]) == "stringfile_tester":
            stringfile_tester_main()
        elif str(sys.argv[1]) == "gml":
            gml = stringfile_to_gml("blackbox/output/CC(C)C(C)C(C)N=C([O-])OC=O_5ae5/reaction0000/stringfile.xyz0000")
            print(gml)
        else:
            print("derp")
        #g = smiles("CCO")                       # molecule to test reaction on
        #xyz_string = mod_to_xyz(g, to_file=False)             # convert molecule for zstruct to understand it
        #run_zstruct_and_gsm(xyz_string)

        #gml_string, atom_core, energy_curve = stringfile_to_rdkit('blackbox/output/5af9b18e744943acab7bffa4d3845c4d/stringfiles/stringfile.xyz0003', visualize=True)
        '''
        g = graphGMLString(gml_string)
        molecule, lookup_dict = make_cut_molecule(g, atom_core)
        cuts = set()
        find_all_cuts(molecule, cuts, lookup_dict, 0)

        # iterate over all cuts and do:
        new_gml_string = make_cut(g, cuts, molecule)
        new_g = graphGMLString(gml_string)
        mod_to_xyz(new_g, to_file=True)
        # run zstruct with new molecule.xyz (molecule.frozen now contains all atoms that are not in the core)
        with open('blackbox/scratch/stringfiles/stringfile.xyz0000') as fi:
            ct = fi.readlines()
        curve = read_energy_profiles(ct)
        x = root_mean_square(energy_curve, curve)   # based on value decide if curve is the same
        '''
