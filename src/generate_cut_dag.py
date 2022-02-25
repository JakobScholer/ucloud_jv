from src.cut_dag import make_childs, make_root, run_blackbox, CutDag, CutDagNode
from src.energy_curve_comparison import root_mean_square
from src.stringfile_helper_functions import read_energy_profiles
from src.visualizers import visualize_cut_dag

from os import listdir
from os.path import isfile, isdir
from multiprocessing import Queue
from portalocker import lock, unlock, LOCK_EX


# generate the empty dag
def generate_empty_dag(stringfile: str, debug: bool=False):
    cut_dag = make_root(stringfile, debug)
    if cut_dag is None:
        return None

    # loop as while there is new nodes in the cut dag
    current_layer = [cut_dag.layers[0][0]] # insert first the root
    next_layer = []
    layer = 0
    while len(current_layer) > 0 or len(next_layer) > 0:
        next_layer = next_layer + make_childs(cut_dag, current_layer[0], layer) # add new nodes to next layer iteration
        current_layer.pop(0) # remove already visited node
        if len(current_layer) == 0: # go to next layer if no more nodes to work with
            current_layer = next_layer.copy()
            next_layer = []
            layer += 1
    return cut_dag

def generate_dag_data(task_queue: Queue, cut_dag: CutDag, stringfile: str, overall_folder: str, reaction_folder: str):
    assigned_tasks = 0
    # check if all cuts have been performed for reaction
    if isfile(f"{overall_folder}/{reaction_folder}/done.txt"):
        with open(f"{overall_folder}/{reaction_folder}/done.txt") as f:
            correct_cut_amount = int(f.readline().rstrip())
        real_cut_amount = 0
        for element in listdir(f"{overall_folder}/{reaction_folder}"):
            if isdir(f"{overall_folder}/{reaction_folder}/{element}"):
                if isfile(f"{overall_folder}/{reaction_folder}/{element}/stringfile.xyz"):
                    real_cut_amount += 1
        if correct_cut_amount == real_cut_amount:
            print(f"skipping {overall_folder}/{reaction_folder}")
            return 0
    # perform cuts
    for k in cut_dag.layers.keys():
        if k > 0:
            for i in range(len(cut_dag.layers[k])):
                assigned_tasks += 1
                node = cut_dag.layers[k][i]
                task_queue.put((dag_point_task, (stringfile, overall_folder, reaction_folder, node)))  # insert new tasks
    if not isfile(f"{overall_folder}/{reaction_folder}/done.txt"):
        with open(f"{overall_folder}/{reaction_folder}/done.txt", "w") as f:
            f.write(f"{assigned_tasks}\n")
    return assigned_tasks


def dag_point_task(stringfile: str, overall_folder: str, reaction_folder: str, node: CutDagNode):
    stringfile_path, error_code, error_message = run_blackbox(stringfile, overall_folder, node.cuts, reaction_folder) # call black box
    node.stringfile = stringfile_path

    if error_code: # the data is not usefull insert in no reaction list file
        cut_reaction = ""
        for c in sorted(node.cuts):
            cut_reaction += str(c) + "_"
        cut_reaction = cut_reaction[:-1]
        with open(f"{overall_folder}/{reaction_folder}/no_reaction.txt", 'a') as f:
            lock(f, LOCK_EX)
            f.write(f"{cut_reaction},{error_message}\n")
            unlock(f)

def read_dag_data(cut_dag, reaction_folder_path):
    # check if any data exist
    if len(listdir(reaction_folder_path)) < 4:
        return

    if isfile(f"{reaction_folder_path}/no_reaction.txt"): # get no reaction cuts from the cut dag file no_reactions.txt
        with open(f"{reaction_folder_path}/no_reaction.txt", 'r') as f:
            no_reaction_list = f.read().splitlines() # get the entire file as a list
        # split up the cut folder names and the error error_message
        no_reaction_folders = []
        no_reaction_errors = {}
        for element in no_reaction_list:
            data = element.split(",")
            no_reaction_folders.append(data[0])
            no_reaction_errors[data[0]] = data[1]
    else:
        no_reaction_folders = []
        no_reaction_errors = {}

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
                if isfile(f"{reaction_folder_path}/{cut_folder}/stringfile.xyz") and cut_folder not in no_reaction_folders:
                    node.stringfile = f"{reaction_folder_path}/{cut_folder}/stringfile.xyz"
                    node.energy = read_energy_profiles(node.stringfile)
                    node.RMS = root_mean_square(cut_dag.layers[0][0].energy, node.energy)
                else:
                    node.stringfile = no_reaction_errors.get(cut_folder)


def cut_dag_setup(stringfile: str, debug: bool=False):
    # from stringfile path, get overall path and reaction_folder
    split_path = stringfile.rsplit("/")
    overall_path = ""
    for i in range(len(split_path) - 2):  # make the complete overall path down to the folder before the reaction folder
        overall_path += split_path[i] + "/"
    overall_path = overall_path[:-1]
    split_path.reverse()  # reverse the list to get the folders easy
    reaction_folder = split_path[1]
    # generate empty dag and check if its done
    if debug:
        print("generate empty dag: Start")
    cut_dag = generate_empty_dag(stringfile, debug)
    if cut_dag is None:
        print("ERROR no cut dag for " + str(stringfile))
        return None
    if debug:
        print("generate empty dag: done")
    return cut_dag, overall_path, reaction_folder


def make_cut_dag(task_queue: Queue, stringfile: str, debug: bool = False):
    cut_dag, overall_path, reaction_folder = cut_dag_setup(stringfile, debug)
    assigned_tasks = generate_dag_data(task_queue, cut_dag, stringfile, overall_path, reaction_folder)
    return assigned_tasks


def show_cut_dag(stringfile: str, visual_cut_dag: bool=False, visual_stringfiles: bool=False, debug: bool = False):
    cut_dag, overall_path, reaction_folder = cut_dag_setup(stringfile, debug)
    if debug:
        print("read dag data: start")
    read_dag_data(cut_dag, f"{overall_path}/{reaction_folder}")
    if debug:
        print("read dag data: done")
    # visualise every stringfile in cut dag data
    if visual_stringfiles:
        if debug:
            print("visualise stringfiles: start")
        from src.visualize_stringfile import visualise_stringfiles
        visualise_stringfiles(f"{overall_path}/{reaction_folder}")
        if debug:
            print("visualise stringfiles: done")
    # visualize the cut dag
    if visual_cut_dag:
        if debug:
            print("visualise cut dag: start")
        visualize_cut_dag(cut_dag)
        if debug:
            print("visualise cut dag: Done")
