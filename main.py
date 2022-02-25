from src.stringfile_to_gml import stringfile_to_gml
from src.generate_cut_dag import show_cut_dag
from src.smiles_to_reactions import make_reactions, make_single_reaction
from src.visualize_stringfile import visualize_2D, visualise_stringfiles
from src.visualizers import visualize_energy_curves, energy_curve_all_reactions

import time
from shutil import copytree
from sys import argv

if __name__ == '__main__':
    start_time = time.time()
    if len(argv) == 2:
        # ---------------visualizers---------------visualizers---------------visualizers---------------visualizers---------------
        if str(argv[1]) == "ec":                        # show energy curves for stringfiles of reaction
            visualize_energy_curves(folder="blackbox/output/nikolai_stringfile4/reaction0000")
        elif str(argv[1]) == "ec_all_reactions":        # show energy curves for all reactions
            energy_curve_all_reactions(folder="blackbox/output/CC=CC=CC=CC_not_done")
        elif str(argv[1]) == "img_all_stringfiles":     # create image for all stringfiles in folder
            visualise_stringfiles(overall_folder="blackbox/output/CC(CCC)=CC=CC=CO_04a5")
        elif str(argv[1]) == "img_stringfile":          # create image for specific stringfile
            visualize_2D(stringfile_path="blackbox/output/CC(CCC)=C(CCC)C=CO_164c/reaction0001/stringfile.xyz0001", image_path="blackbox/output/CC(CCC)=C(CCC)C=CO_164c/reaction0001")
        elif str(argv[1]) == "make_cut_dag":  # create cut dag for specific stringfile
            str_file = "blackbox/output/CC=CC=CC=CO_7498/reaction0001/stringfile.xyz"
            make_single_reaction(stringfile=str_file, number_of_processes=4, debug=False)   # computes cuts for single reaction
            show_cut_dag(stringfile=str_file, visual_cut_dag=True, visual_stringfiles=True, debug=False)    # shows cut dag for single reaction

        # ---------------Calc---------------Calc---------------Calc---------------Calc---------------Calc---------------
        elif str(argv[1]) == "smiles_to_reactions_bb":  # Compute all reactions as stringfiles from smiles
            make_reactions(blackbox=True, string_data=["CC=C"], max_energy=200, frozen=[], number_of_processes=8, debug=False)
        elif str(argv[1]) == "smiles_to_reactions_bb_multirun":  # Compute all reactions as stringfiles from smiles
            for i in range(20):
                copytree("blackbox/output/CC=CCCCC_sigmatropicForced_3/original", f"blackbox/output/CC=CCCCC_sigmatropicForced_3/{i}")
                make_reactions(blackbox=True, string_data=f"CC=CCCCC_sigmatropicForced_3/{i}", max_energy=200, frozen=[])
        elif str(argv[1]) == "smiles_to_reactions_bb_continue":  # Compute all reactions as stringfiles from smiles on existing folder
            make_reactions(blackbox=True, string_data="CC=C_4270", max_energy=200, number_of_processes=4)
        elif str(argv[1]) == "smiles_to_reactions_nb":  # Read all reactions from stringfiles
            make_reactions(blackbox=False, string_data="blackbox/output/CCCO[O]_c014", max_energy=50, visual_cut_dag=True, visual_stringfiles=True, debug=False)
        elif str(argv[1]) == "gml":                     # create gml rule
            gml = stringfile_to_gml(filename="blackbox/output/CC(C)C(C)C(C)N=C([O-])OC=O_5ae5/reaction0000/stringfile.xyz0000")
            print(gml)
        else:
            print("input argument is unknown")
    else:
        print("input argument required")
    print("Total run time took:")
    print("--- %s seconds ---" % (time.time() - start_time))
