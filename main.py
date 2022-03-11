from src.stringfile_to_gml import stringfile_to_gml
from src.generate_cut_dag import show_cut_dag
from src.smiles_to_reactions import make_reactions, make_single_reaction
from src.visualize_stringfile import visualize_2D, visualise_stringfiles
from src.visualizers import visualize_energy_curves, energy_curve_all_reactions
from src.blackbox import run_gsm_round_initial

import time
from shutil import copytree
from sys import argv

if __name__ == '__main__':
    start_time = time.time()
    if len(argv) == 2:
        # ---------------visualizers---------------visualizers---------------visualizers---------------visualizers---------------
        if str(argv[1]) == "ec":                        # show energy curves for stringfiles of reaction
            visualize_energy_curves(folder="blackbox/output/C=O_+_C(O)=CO_custom/reaction0043")
        elif str(argv[1]) == "ec_all_reactions":        # show energy curves for all reactions
            max_peak = 200 # max value for the peak of the energi curveS
            remove_no_reactions = True # remove all reactions where the educt and product are the same
            peak_end_diff_max = 5 # max diff between the peak and the end point of the curve
            max_endpoint = 50 # max value for the end point on the curve
            energy_curve_all_reactions("blackbox/output/C(O)=CO_+_C(C(CO)O)=O_decc/reaction0140", max_peak, remove_no_reactions, peak_end_diff_max, max_endpoint)
        elif str(argv[1]) == "img_all_stringfiles":     # create image for all stringfiles in folder
            visualise_stringfiles(overall_folder="blackbox/output/C(O)=CO_+_C(C(CO)O)=O_decc/reaction0140")
        elif str(argv[1]) == "img_stringfile":          # create image for specific stringfile
            visualize_2D(stringfile_path="blackbox/output/CC(N)=C(N)C(O)=C(N)C(O)=C(C)O_220d/reaction0027/stringfile.xyz", image_path="CC(N)=C(N)C(O)=C(N)C(O)=C(C)O_220d/reaction0027")
        elif str(argv[1]) == "show_cut_dag":  # create cut dag for specific stringfile
            str_file = "blackbox/output/C(O)=CO_+_C(C(CO)O)=O_decc/reaction0140/stringfile.xyz"
            show_cut_dag(stringfile=str_file, visual_cut_dag=True, visual_stringfiles=False, debug=False)

        # ---------------Calc---------------Calc---------------Calc---------------Calc---------------Calc---------------
        elif str(argv[1]) == "smiles_to_reactions_bb":  # Compute all reactions as stringfiles from smiles
            make_reactions(blackbox=True, string_data=['C=O', 'C(O)=CO'], max_energy=200, frozen=[], number_of_processes=4, coordinate_generator="openbabel", debug=True)
        elif str(argv[1]) == "make_cut_dag":
            make_single_reaction(stringfile = "blackbox/output/C(O)=CO_+_C(C(CO)O)=O_decc/reaction0140/stringfile.xyz", number_of_processes = 4, debug = False)
        elif str(argv[1]) == "smiles_to_reactions_bb_multirun":  # Compute all reactions as stringfiles from smiles
            for i in range(20):
                copytree("blackbox/output/CC(N)=C(N)C(O)=C(N)C(O)=C(C)O_220d - reaction21 - 2_4_6_8_9_10_11_12_13 - second/original", f"blackbox/output/CC(N)=C(N)C(O)=C(N)C(O)=C(C)O_220d - reaction21 - 2_4_6_8_9_10_11_12_13 - second/{i}")
                #make_single_reaction(stringfile=f"blackbox/output/test/{i}/stringfile.xyz", number_of_processes = 1, debug = False)
                run_gsm_round_initial(output_folder = f"blackbox/output/CC(N)=C(N)C(O)=C(N)C(O)=C(C)O_220d - reaction21 - 2_4_6_8_9_10_11_12_13 - second/{i}", reaction_folder = "", logfile = False)
        elif str(argv[1]) == "smiles_to_reactions_bb_continue":  # Compute all reactions as stringfiles from smiles on existing folder
            make_reactions(blackbox=True, string_data="C=O_+_C(O)=CO_custom", max_energy=200, number_of_processes=4, coordinate_generator="openbabel")
        elif str(argv[1]) == "gml":                     # create gml rule
            gml = stringfile_to_gml(filename="blackbox/output/CC(C)C(C)C(C)N=C([O-])OC=O_5ae5/reaction0000/stringfile.xyz0000")
            print(gml)
        else:
            print("input argument is unknown")
    else:
        print("input argument required")
    print("Total run time took:")
    print("--- %s seconds ---" % (time.time() - start_time))

    '''
    from src.stringfile_tester import check_product
    from src.stringfile_to_rdkit import stringfile_to_rdkit
    from src.cut_molecule import make_cut_molecule, make_cut
    from rdkit.Chem import RWMol, AddHs, MolFromSmiles, MolToXYZBlock
    from src.stringfile_tester import check_initial_file

    original_strfile = "blackbox/output/CC(N)=C(N)C(O)=C(N)C(O)=C(C)O_220d/reaction0021/stringfile.xyz"
    modified_strfile = "blackbox/output/CC(N)=C(N)C(O)=C(N)C(O)=C(C)O_220d/reaction0021/10_12/stringfile.xyz"

    rdk_mol, atom_core, energy_curve = stringfile_to_rdkit(original_strfile, visualize=False)
    molecule, lookup_dict = make_cut_molecule(rdk_mol, atom_core)
    modified_mol, order = make_cut(rdk_mol, {4,8,13}, molecule, lookup_dict)
    #print(MolToXYZBlock(modified_mol))
    #print(check_product(original_strfile, modified_strfile, {10,12}, order, molecule, lookup_dict, rdk_mol))
    for i in range(595):
        number = str(i)
        while len(number) < 4:
            number = "0" + number

        if check_initial_file(f"blackbox/output/C(O)=CO_+_C(C(CO)O)=O_decc/reaction{number}/initial0000.xyz"):
            print(number)
    '''
