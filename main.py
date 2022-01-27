from sys import argv
from src.stringfile_to_gml import stringfile_to_gml
from src.generate_cut_dag import make_cut_dag, visualise_stringfiles
from src.smiles_to_reactions import make_reactions
from src.visualize_stringfile import visualize_2D
from src.visualizers import visualize_energy_curves

if __name__ == '__main__':
    if len(argv) == 2:
        # ---------------visualizers---------------visualizers---------------visualizers---------------visualizers---------------
        if str(argv[1]) == "ec":                        # show energy curves for stringfiles of reaction
            visualize_energy_curves(folder="blackbox/output/CCCO[O]_c014/reaction0003")
        elif str(argv[1]) == "img_all_stringfiles":     # create image for all stringfiles in folder
            visualise_stringfiles(overall_folder="blackbox/output/CCCCCO[O-]_e08a")
        elif str(argv[1]) == "img_stringfile":          # create image for specific stringfile
            visualize_2D(stringfile_path="blackbox/output/CCCCCO[O-]_e08a/reaction0001/stringfile.xyz0001", image_path="blackbox/output/CC(C)C(C)C(C)N=C([O-])OC=O_6d1e/reaction0004")
        elif str(argv[1]) == "make_cut_dag":  # create cut dag for specific stringfile
            make_cut_dag(blackbox=False, stringfile="blackbox/output/CCCCCO[O-]_e08a/reaction0001/stringfile.xyz0001",visual_cut_dag=True , visual_stringfiles=False, DEBUG_MODE=True)

        # ---------------Calc---------------Calc---------------Calc---------------Calc---------------Calc---------------
        elif str(argv[1]) == "smiles_to_reactions_bb":  # Compute all reactions as stringfiles from smiles
            make_reactions(blackbox=True, string_data=["CCCCCO[O-]"], max_energy=200)
        elif str(argv[1]) == "smiles_to_reactions_bb_continue":  # Compute all reactions as stringfiles from smiles on existing folder
            make_reactions(blackbox=True, string_data="CCCCCO[O-]_e08a", max_energy=200)
        elif str(argv[1]) == "smiles_to_reactions_nb":  # Read all reactions from stringfiles
            make_reactions(blackbox=False, string_data="blackbox/output/CCCO[O]_c014", max_energy=50, visual_cut_dag=True, visual_stringfiles=True, debug=False)
        elif str(argv[1]) == "smiles_to_reactions_bb_skip":  # Compute all reactions from pre-computed stringfiles
            make_reactions(blackbox=True, string_data="blackbox/output/CCCO[O]_c014", generate_initial_stringfiles=False, max_energy=200)
        elif str(argv[1]) == "gml":                     # create gml rule
            gml = stringfile_to_gml(filename="blackbox/output/CC(C)C(C)C(C)N=C([O-])OC=O_5ae5/reaction0000/stringfile.xyz0000")
            print(gml)
        else:
            print("input argument is unknown")
    else:
        print("input argument required")
