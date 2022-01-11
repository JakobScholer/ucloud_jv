from src.cut_dag import make_childs_mp, insert_childs_mp, make_root, run_blackbox
from multiprocessing import Process, Queue, freeze_support
from src.energy_curve_comparison import root_mean_square
from src.stringfile_helper_functions import read_energy_profiles
from src.stringfile_to_rdkit import stringfile_to_rdkit
from src.cut_molecule import make_cut_molecule, find_all_cuts, make_cut
from src.visualizers import visualize_cut_dag
from src.visualize_stringfile import visualize_2D
from src.stringfile_tester import check_product
import time
from os import listdir
from os.path import isdir

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

# generate the empty dag using multiprocessing.
def generate_empty_dag_mp(stringfile, DEBUG_MODE: bool=False):
    cd = make_root(stringfile, DEBUG_MODE)
    if cd is None:
        return None, 0

    freeze_support()
    NUMBER_OF_PROCESSES = 1

    # Create queues
    task_queue = Queue()
    done_queue = Queue()
    if DEBUG_MODE:
        print("    MAKE ROOT! start")
    # add first task
    root = cd.layers[0][0]
    root_task = (make_childs_mp, (stringfile, set(), (0,0)))
    task_queue.put(root_task)
    if DEBUG_MODE:
        print("    MAKE ROOT! stop")

    # Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue, done_queue)).start()

    # wait for porcesses to end
    wait_for_end = True
    tasks_sent = 1
    tasks_completed = 0
    if DEBUG_MODE:
        print("    GENERATE CUT DAG! start")
    while wait_for_end:
        if DEBUG_MODE:
            print("        TASK STATE INFO!")
            print("            tasks_sent: " + str(tasks_sent))
            print("            tasks_completed: " + str(tasks_completed))
        if tasks_sent == tasks_completed:
            wait_for_end = False
        elif done_queue.empty() == True:
            if DEBUG_MODE:
                print("        No new data recived")
            time.sleep(0.2)
        else:
            if DEBUG_MODE:
                print("        DATA REACIVED! start")
            while done_queue.empty() == False:
                child_infos = done_queue.get() # get info
                tasks_completed += 1
                if DEBUG_MODE:
                    print("            got info: " + str(removeDuplicates(child_infos[0])))
                tasks = insert_childs_mp(stringfile, cd, removeDuplicates(child_infos[0]), child_infos[1]) # insert child
                if len(tasks) > 0:
                    for t in tasks: # add new tasks to the queue
                        if DEBUG_MODE:
                            print("            sending info: " + str(t[1][1]))
                        task_queue.put(t)
                    tasks_sent += len(tasks)
                if DEBUG_MODE:
                    print("        DATA REACIVED! stop")
    if DEBUG_MODE:
        print("    GENERATE CUT DAG! stop")

    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')

    return cd, tasks_sent

def generate_dag_data_mp(cd, tasks_counter, stringfile, overall_folder, reaction_folder, DEBUG_MODE: bool=False):

    freeze_support()
    NUMBER_OF_PROCESSES = 1
    # Create queues
    task_queue = Queue()
    done_queue = Queue()
    # Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue, done_queue)).start()

    # make all tasks for blackbox
    tasks_bx = []
    for k in cd.layers.keys():
        if k > 0:
            for i in range(len(cd.layers[k])):
                task_queue.put((run_blackbox, (stringfile, overall_folder, cd.layers[k][i].cuts, (k,i), reaction_folder))) # insert new tasks

    if DEBUG_MODE:
        print("    que size: " + str(task_queue.qsize()) + " and needed tasks: " + str(tasks_counter))

    tasks_completed = 1 # root is already done as a task
    # make all tasks for the blackbox
    while tasks_completed != tasks_counter: #there is still tasks to perform
        if DEBUG_MODE:
            print("    tasks completed: " + str(tasks_completed))
            print("    tasks sent: " + str(tasks_counter))
            print("    tasks in queue: " + str(task_queue.qsize()))
        if not done_queue.empty(): # insert return data in format (stringfile, Energy, placement)
            while not done_queue.empty(): # empty the gueue
                if DEBUG_MODE:
                    print("    whuue got some BX data")
                data = done_queue.get()
                node = cd.layers[data[1][0]][data[1][1]]
                node.stringfile = data[0]
                if not data[0] == "NO REACTION":
                    node.energy = read_energy_profiles(data[0])
                    node.RMS = root_mean_square(cd.layers[0][0].energy, node.energy)
                tasks_completed += 1 # increment the number of tasks needed to be done
        else: # else wait a litle and check again
            if DEBUG_MODE:
                print("    sleep sleep")
            time.sleep(2)
            if DEBUG_MODE:
                print("    waky waky")
    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')
    if DEBUG_MODE:
        print("    done!")


def read_dag_data(cut_dag, reaction_folder):
    # check if any data exist
    if len(listdir(reaction_folder)) < 4:
        return "NO DATA"
    # gennem gå hele daggen
    for k in cut_dag.layers.keys():
        if k > 0:
            for i in range(len(cut_dag.layers[k])):
                node = cut_dag.layers[k][i] # node
                # make folder Name
                folder = ""
                for c in sorted(node.cuts):
                    folder += str(c) + "_"
                folder = folder[:-1]
                # go to folder and find stringfile
                no_stringfile = True
                for file in listdir(reaction_folder + "/" + folder):
                    if "stringfile" in file: # if stringfile insert alle data
                        # make cut molecules
                        stringfile_path = reaction_folder + "/" + folder + "/" + str(file)
                        rdk_mol, atom_core, energy_curve = stringfile_to_rdkit(stringfile_path, False)
                        molecule, lookup_dict = make_cut_molecule(rdk_mol, atom_core)
                        modified_mol, order = make_cut(rdk_mol, node.cuts, molecule, lookup_dict)

                        if check_product(cut_dag.layers[0][0].stringfile, stringfile_path, node.cuts, order, molecule, lookup_dict):
                            node.stringfile = stringfile_path
                            node.energy = read_energy_profiles(node.stringfile)
                            node.RMS = root_mean_square(cut_dag.layers[0][0].energy, node.energy)
                            no_stringfile = False
                if no_stringfile:
                    node.stringfile = "NO REACTION"
    return "done"

def visualise_stringfiles(folder, DEBUG_MODE: bool=False):
    # go over each cut folder
    for folder in listdir(folder):
        folder_name = folder + "/" + str(folder)
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
            visualize_2D(folder + "/" + str(folder), folder)


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
    cd, tasks_counter = generate_empty_dag_mp(stringfile, DEBUG_MODE)
    if cd is None:
        print("ERROR no cut dag for " + str(stringfile))
        return None
    if DEBUG_MODE:
        print("generate empty dag: done")

    if blackbox: # Run black box
        if DEBUG_MODE:
            print("Blackbox run: start")
        generate_dag_data_mp(cd, tasks_counter, stringfile, overall_path, reaction_folder, DEBUG_MODE)
        if DEBUG_MODE:
            print("Blackbox run: done")
    else: # read data drom folder
        if DEBUG_MODE:
            print("read dag data: start")
        status = read_dag_data(cd, overall_path + "/" + reaction_folder)
        if DEBUG_MODE:
            print("read dag data: " + status)
    # visualise every stringfile in cut dag data
    if visual_stringfiles:
        if DEBUG_MODE:
            print("visualise stringfiles: start")
        visualise_stringfiles(overall_path + "/" + reaction_folder)
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
