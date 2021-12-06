from src.cut_dag import CutDagNode, CutDag, make_childs_mp, insert_childs_mp, make_root, run_blackbox
from multiprocessing import Process, Queue, current_process, freeze_support
from src.root_mean_square import root_mean_square
from src.cut_molecule import make_cut, make_cut_molecule
from src.stringfile_to_rdkit import stringfile_to_rdkit, read_energy_profiles
from src.zstruct_and_gsm import run_zstruct_and_gsm
from src.stringfile_tester import check_product
#from mod import *
from igraph import *
import plotly.graph_objects as go
import time

def removeDuplicates(arr):
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

# the function that generate the full cut dag
def make_cut_dag(stringfile, overall_folder, reaction_folder, DEBUG_MODE: bool = False):
    NUMBER_OF_PROCESSES = 4

    #stringfile = "xyz_test_files/GCD_test_files/stringfile.xyz0009"
    graph = False
    cd = make_root(stringfile, graph)
    if cd is None:
        return None

    # Create queues
    task_queue = Queue()
    done_queue = Queue()
    if DEBUG_MODE:
        print("MAKE ROOT! start")
    # add first task
    root = cd.layers[0][0]
    root_task = (make_childs_mp, (stringfile, set(), (0,0)))
    task_queue.put(root_task)
    if DEBUG_MODE:
        print("MAKE ROOT! stop")

    # Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue, done_queue)).start()

    # wait for porcesses to end
    wait_for_end = True
    tasks_sent = 1
    tasks_completed = 0
    if DEBUG_MODE:
        print("GENERATE CUT DAG! start")
    while wait_for_end:
        if DEBUG_MODE:
            print("    TASK STATE INFO!")
            print("        tasks_sent: " + str(tasks_sent))
            print("        tasks_completed: " + str(tasks_completed))
        if tasks_sent == tasks_completed:
            wait_for_end = False
        elif done_queue.empty() == True:
            if DEBUG_MODE:
                print("    No new data recived")
            time.sleep(0.2)
        else:
            if DEBUG_MODE:
                print("    DATA REACIVED! start")
            while done_queue.empty() == False:
                child_infos = done_queue.get() # get info
                tasks_completed += 1
                if DEBUG_MODE:
                    print("        got info: " + str(removeDuplicates(child_infos[0])))
                tasks = insert_childs_mp(stringfile, cd, removeDuplicates(child_infos[0]), child_infos[1]) # insert child
                if len(tasks) > 0:
                    for t in tasks: # add new tasks to the queue
                        if DEBUG_MODE:
                            print("        sending info: " + str(t[1][1]))
                        task_queue.put(t)
                    tasks_sent += len(tasks)
                if DEBUG_MODE:
                    print("    DATA REACIVED! stop")
    if DEBUG_MODE:
        print("GENERATE CUT DAG! stop")

    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')

    if True:
        return cd

    # make all tasks for blackbox
    tasks_bx = []
    for k in cd.layers.keys():
        if k > 0:
            for i in range(len(cd.layers[k])):
                task_queue.put((run_blackbox, (stringfile, overall_folder, cd.layers[k][i].cuts, (k,i), reaction_folder))) # insert new tasks
    #task_counter = len(tasks_bx)

    if DEBUG_MODE:
        print("que size: " + str(task_queue.qsize()) + " and needed tasks: " + str(tasks_sent))

    tasks_completed = 1 # root is already done as a task
    # make all tasks for the blackbox
    while tasks_completed != tasks_sent: #there is still tasks to perform
        if DEBUG_MODE:
            print("tasks completed: " + str(tasks_completed))
            print("tasks sent: " + str(tasks_sent))
            print("tasks in queue: " + str(task_queue.qsize()))
        if done_queue.empty() == False: # insert return data in format (stringfile, Energy, placement)
            while done_queue.empty() == False: # empty the gueue
                if DEBUG_MODE:
                    print("whuue got some BX data")
                data = done_queue.get()
                #print("------------------------")
                #print(data)
                node = cd.layers[data[1][0]][data[1][1]]
                node.stringfile = data[0]
                if not data[0] == "NO REACTION":
                    with open(data[0]) as f:
                        data_stringfile = f.readlines()
                    node.energy = read_energy_profiles(data_stringfile)
                    node.RMS = root_mean_square(cd.layers[0][0].energy, node.energy)
                tasks_completed += 1 # increment the number of tasks needed to be done
        else: # else wait a litle and check again
            if DEBUG_MODE:
                print("sleep sleep")
            time.sleep(2)
            if DEBUG_MODE:
                print("waky waky")
    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')
    if DEBUG_MODE:
        print("done!")

    if DEBUG_MODE:
        for k in cd.layers.keys():
            for node in cd.layers[k]:
                print("layer " + str(k))
                print("node with cuts " + str(node.cuts))

    return cd


def visualizer(cut_dag, borderline_value):
    cut_option_y_green = []
    cut_option_x_green = []
    cut_option_y_red = []
    cut_option_x_red = []
    cut_option_y_black = []
    cut_option_x_black = []
    cut_option_y = []
    cut_option_x = []
    cut_info = []
    stringfiles = []
    bond_y = []
    bond_x = []
    for layer in cut_dag.layers.keys():
        layer_length = len(cut_dag.layers.get(layer))
        for position in range(layer_length):
            print(cut_dag.layers.get(layer)[position].RMS)
            if cut_dag.layers.get(layer)[position].stringfile == "NO REACTION":
                cut_option_y_black.append(layer * 10)
                cut_option_x_black.append(position * 10 - (layer_length * 10) / 2)
            elif cut_dag.layers.get(layer)[position].RMS <= borderline_value:
                cut_option_y_green.append(layer * 10)
                cut_option_x_green.append(position * 10 - (layer_length * 10) / 2)
            elif cut_dag.layers.get(layer)[position].RMS > borderline_value:
                cut_option_y_red.append(layer * 10)
                cut_option_x_red.append(position * 10 - (layer_length * 10) / 2)
            cut_option_y.append(layer * 10)
            cut_option_x.append(position * 10 - (layer_length * 10) / 2)
            cut_info.append(cut_dag.layers.get(layer)[position].cuts)
            stringfiles.append(cut_dag.layers.get(layer)[position].stringfile)
            for child_position in cut_dag.layers.get(layer)[position].childs:
                bond_y += [layer * 10, (layer+1) * 10, None]
                bond_x += [position * 10 - (layer_length * 10)/2, child_position * 10 - (len(cut_dag.layers.get(layer+1)) * 10)/2, None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=bond_x,
                             y=bond_y,
                             mode='lines',
                             name='connections',
                             line=dict(color='rgb(210,210,210)', width=3),
                             hoverinfo='skip'
                             ))

    fig.add_trace(go.Scatter(x=cut_option_x,
                             y=cut_option_y,
                             mode='markers',
                             name='Cuts',
                             marker=dict(symbol='circle-dot',
                                         size=20,
                                         color='#800080'
                                         ),
                             text=cut_info,
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x_red,
                             y=cut_option_y_red,
                             mode='markers',
                             name='Cuts',
                             marker=dict(symbol='circle-dot',
                                         size=20,
                                         color='#FF0000'
                                         ),
                             text=cut_info,
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x_green,
                             y=cut_option_y_green,
                             mode='markers',
                             name='Cuts',
                             marker=dict(symbol='circle-dot',
                                         size=20,
                                         color='#008000'
                                         ),
                             text=cut_info,
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x_black,
                             y=cut_option_y_black,
                             mode='markers',
                             name='Cuts',
                             marker=dict(symbol='circle-dot',
                                         size=20,
                                         color='#000000'
                                         ),
                             text=cut_info,
                             hoverinfo='skip',
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x,
                             y=cut_option_y,
                             mode='text',
                             name='Cut info',
                             text=cut_info,
                             hoverinfo='skip',
                             textfont_size=10
                             ))
    fig.add_trace(go.Scatter(x=cut_option_x,
                             y=[k - 1 for k in cut_option_y],
                             mode='text',
                             name='stringfile',
                             text=stringfiles,
                             hoverinfo='skip',
                             textfont_size=10
                             ))

    fig.show()

def generate_cut_dag_main(stringfile, overall_path, reaction_folder, debug_mode: bool = False):
    freeze_support()
    cut_dag = make_cut_dag(stringfile, overall_path, reaction_folder, debug_mode)

    if cut_dag is not None:
        print("SUCCES!")
        #visualizer(cut_dag, 100)
    else:
        print("ERROR not cut dag for " + str(stringfile))
