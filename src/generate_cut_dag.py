from src.cut_dag import make_childs, make_root, run_blackbox
from src.energy_curve_comparison import root_mean_square
from src.stringfile_helper_functions import read_energy_profiles
from src.visualizers import visualize_cut_dag
from src.visualize_stringfile import visualize_2D
from os import listdir
from os.path import isdir, isfile

def removeDuplicates(arr): # midlertidig methode. lav core_ring_check i cut molecule for at fikse det.
    temp = []
    for e in arr:
        if e not in temp:
            temp.append(e)
    return temp

# Function run by worker processes
def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = func(*args)
        output.put(result)

# generate the empty dag
def generate_empty_dag(stringfile, DEBUG_MODE: bool=False):
    cd = make_root(stringfile, DEBUG_MODE)
    if cd is None:
        return None, 0

    # loop as while there is new nodes in the cut dag
    current_layer = [cd.layers[0][0]] # insert first the root
    next_layer = []
    layer = 0
    while len(current_layer) > 0 or len(next_layer) > 0:
        next_layer = next_layer + make_childs(cd, current_layer[0], layer) # add new nodes to next layer iteration
        current_layer.pop(0) # remove already visited node
        if len(current_layer) == 0: # go to next layer if no more nodes to work with
            current_layer = next_layer.copy()
            next_layer = []
            layer += 1
    return cd

def generate_dag_data(cd, stringfile, overall_folder, reaction_folder, DEBUG_MODE: bool=False):
    # make all tasks for blackbox
    tasks_bx = []
    for k in cd.layers.keys():
        if k > 0:
            for i in range(len(cd.layers[k])):
                data = run_blackbox(stringfile, overall_folder, cd.layers[k][i].cuts, (k,i), reaction_folder) # call black box
                node = cd.layers[k][i]
                node.stringfile = data[0]

                if data[0] == "NO REACTION": # the data is not usefull insert in no reaction list file
                    cut_reaction = ""
                    for c in sorted(node.cuts):
                        cut_reaction += str(c) + "_"
                    cut_reaction = cut_reaction[:-1]
                    print(f"making folder at {overall_folder}/{reaction_folder}/no_reaction.txt")
                    with open(f"{overall_folder}/{reaction_folder}/no_reaction.txt", 'a') as f:
                        f.write(f"{cut_reaction}\n")
                else: # teh data is usefull!
                    node.energy = read_energy_profiles(data[0])
                    node.RMS = root_mean_square(cd.layers[0][0].energy, node.energy)

def read_dag_data(cut_dag, reaction_folder_path):
    # check if any data exist
    if len(listdir(reaction_folder_path)) < 4:
        return "NO DATA"

    if isfile(f"{reaction_folder_path}/no_reaction.txt"): # get no reaction cuts from the cut dag file no_reactions.txt
        with open(f"{reaction_folder_path}/no_reaction.txt", 'r') as f:
            no_reaction_list = f.readlines() # get the entire file as a list
    else:
        no_reaction_list = []

    # gennem gå hele daggen
    for k in cut_dag.layers.keys():
        if k > 0:
            for i in range(len(cut_dag.layers[k])):
                node = cut_dag.layers[k][i] # node
                # make folder Name
                cut_folder = ""
                for c in sorted(node.cuts):
                    cut_folder += str(c) + "_"
                cut_folder = cut_folder[:-1]
                # go to folder and find stringfile
                if isfile(f"{reaction_folder_path}/{cut_folder}/stringfile.xyz0000") and cut_folder not in no_reaction_list:
                    node.stringfile = f"{reaction_folder_path}/{cut_folder}/stringfile.xyz0000"
                    node.energy = read_energy_profiles(node.stringfile)
                    node.RMS = root_mean_square(cut_dag.layers[0][0].energy, node.energy)
                else:
                    node.stringfile = "NO REACTION"

def visualise_stringfiles(overall_folder, DEBUG_MODE: bool=False):
    # go over each cut folder
    for folder in listdir(overall_folder):
        folder_name = overall_folder + "/" + str(folder)
        if isdir(folder_name): # its a folder
            for file in listdir(folder_name): # find stringfile
                if "stringfile" in file:
                    if DEBUG_MODE:
                        print("    stringfile path: " + folder_name + "/" + str(file))
                        print("    image path: " + folder_name)
                    visualize_2D(folder_name + "/" + str(file), folder_name)
        elif "stringfile" in folder: # Make a visual version of the original stringfile
            if DEBUG_MODE:
                print("    stringfile path: " + folder_name + "/" + str(file))
                print("    image path: " + folder_name)
            visualize_2D(overall_folder + "/" + str(folder), overall_folder)
            visualize_2D(overall_folder + "/" + str(folder), overall_folder)


def make_cut_dag(blackbox: bool, stringfile, visual_cut_dag: bool=False, visual_stringfiles: bool=False, DEBUG_MODE: bool = False):
    # from stringfile path, get overall path and reaction_folder
    split_path = stringfile.rsplit("/")
    overall_path = ""
    for i in range(len(split_path)-2): # make the complete overall path down to the folder before the reaction folder
        overall_path += split_path[i] + "/"
    overall_path = overall_path[:-1]
    split_path.reverse() # reverse the list to get the folders easy
    reaction_folder = split_path[1]

    # generate empty dag and check if its done
    if DEBUG_MODE:
        print("generate empty dag: Start")
    cd = generate_empty_dag(stringfile, DEBUG_MODE)
    if cd is None:
        print("ERROR no cut dag for " + str(stringfile))
        return None
    if DEBUG_MODE:
        print("generate empty dag: done")

    if blackbox: # Run black box
        if DEBUG_MODE:
            print("Blackbox run: start")
        generate_dag_data(cd, stringfile, overall_path, reaction_folder, DEBUG_MODE)
        if DEBUG_MODE:
            print("Blackbox run: done")
    else: # read data drom folder
        if DEBUG_MODE:
            print("read dag data: start")
        read_dag_data(cd, f"{overall_path}/{reaction_folder}")
        if DEBUG_MODE:
            print("read dag data: done")
    # visualise every stringfile in cut dag data
    if visual_stringfiles:
        if DEBUG_MODE:
            print("visualise stringfiles: start")
        visualise_stringfiles(f"{overall_path}/{reaction_folder}")
        if DEBUG_MODE:
            print("visualise stringfiles: done")
    # visualize the cut dag
    if visual_cut_dag:
        if DEBUG_MODE:
            print("visualise cut dag: start")
        visualize_cut_dag(cd)
        if DEBUG_MODE:
            print("visualise cut dag: Done")
    print("Generate cut dag complete.")
