from src.cut_dag import CutDagNode, CutDag, make_childs_mp, insert_childs_mp, make_root, run_blackbox
from multiprocessing import Process, Queue, current_process, freeze_support
from src.root_mean_square import root_mean_square
from src.cut_molecule import make_cut, make_cut_molecule
from src.stringfile_to_rdkit import stringfile_to_rdkit
from src.zstruct_and_gsm import run_zstruct_and_gsm
from mod import *
from igraph import *
import plotly.graph_objects as go
import time

#
# Function run by worker processes
#

def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = func(*args)
        output.put(result)

def test_black_box(xyz, isomer, order, core):
    print(xyz)
    print(isomer)
    print(order)
    print(core)
    time.sleep(5)
    return ("random_stringfile", [1,2,3,4,5,6,7,8,9,10])

def make_cut_dag():
    NUMBER_OF_PROCESSES = 4

    # make test cut dag
    stringfile = 'xyz_test_files/GCD_test_files/stringfile.xyz0003'
    with open("xyz_test_files/GCD_test_files/ISOMERS0003", "r") as f:
        isomer = f.read()
    print(isomer)
    graph = False
    cd = make_root(stringfile, graph)

    # Create queues
    task_queue = Queue()
    done_queue = Queue()

    print("make task for root")
    # add first task
    root = cd.layers[0][0]
    root_task = (make_childs_mp, (stringfile, set(), (0,0)))
    task_queue.put(root_task)
    print("done")

    # Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue, done_queue)).start()

    # wait for porcesses to end
    wait_for_end = True
    print("staring wait for it")
    while wait_for_end:
        # make sure the quee is empty before stopping porcesses
        if done_queue.empty() and task_queue.empty():
            print("sleep time")
            time.sleep(0.2)
            print("waky waky!")
            if done_queue.empty():
                wait_for_end = False
        else: # there is childs to add the the cut dag
            while done_queue.empty() == False:
                print("whuuee got one bunch of childs")
                child_infos = done_queue.get() # get info
                tasks = insert_childs_mp(stringfile, cd, child_infos[0], child_infos[1]) # insert child
                if len(tasks) > 0:
                    print("added " + str(len(tasks)) + " to the tasks")
                    for t in tasks: # add new tasks to the queue
                        task_queue.put(t)
                else:
                    print("no new tasks added")
    print("Cut dag is generated")

    # make all tasks for blackbox
    tasks_bx = []
    task_counter = len(tasks_bx)
    for k in cd.layers.keys():
        if k > 0:
            for i in range(len(cd.layers[k])):
                #tasks_bx.append((test_blackbox, (stringfile, isomer, cd.layers[k][i].cuts, (k,i)))) # test call
                tasks_bx.append((run_blackbox, (stringfile, isomer, cd.layers[k][i].cuts, (k,i)))) # correct version

    # make all tasks for the blackbox
    while len(tasks_bx) > 0: #there is still tasks to perform
        # insert more tasks if size space
        if task_queue.qsize() < (NUMBER_OF_PROCESSES * 2):
            # insert maximum amount of tasks
            if (NUMBER_OF_PROCESSES * 2) - task_queue.qsize() < len(tasks_bx):
                insert_amount = (NUMBER_OF_PROCESSES * 2) - task_queue.qsize()
            else:
                insert_amount = len(tasks_bx)
            print("inserting " + str(insert_amount) + " new tasks")
            for i in range(insert_amount):
                task_queue.put(tasks_bx[i]) # insert new tasks
            for i in range(insert_amount):
                tasks_bx.pop(i) # remove already inserted tasks
        if done_queue.empty() == False: # insert return data in format (stringfile, Energy, placement)
            while done_queue.empty() == False: # empty the gueue
                print("whuue got some BX data")
                data = done_queue.get()
                node = cd.layers[data[1][0]][data[1][1]]
                node.stringfile = data[0]
                node.energy = read_energy_profiles(data[0])
                node.RMS = root_mean_square(cd.layers[0][0].energy, node.energy)
                task_counter -= task_counter # increment the number of tasks needed to be done
        else: # else wait a litle and check again
            print("sleep sleep")
            time.sleep(2)
            print("waky waky")
    # all tasks has been inserted

    print("no more tasks only wait for data now!")

    # wait for porcesses to end
    while task_counter > 0:
        if done_queue.empty() == False: # insert return data in format (stringfile, Energy, placement)
            while done_queue.empty() == False: # empty the gueue
                print("whuue got some BX data")
                data = done_queue.get()
                node = cd.layers[data[1][0]][data[1][1]]
                node.stringfile = data[0]
                node.energy = read_energy_profiles(data[0])
                node.RMS = root_mean_square(cd.layers[0][0].energy, node.energy)
                task_counter -= task_counter # increment the number of tasks needed to be done
        else: # else wait a litle and check again
            print("sleep sleep")
            time.sleep(2)
            print("waky waky")
    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')
    print("done!")

    for k in cd.layers.keys():
        for node in cd.layers[k]:
            print("layer " + str(k))
            print("node with cuts " + str(node.cuts))

    return cd


def visualizer(cut_dag):
    cut_option_y = []
    cut_option_x = []
    cut_info = []
    bond_y = []
    bond_x = []
    for layer in cut_dag.layers.keys():
        layer_length = len(cut_dag.layers.get(layer))
        for position in range(layer_length):
            print(cut_dag.layers.get(layer)[position].childs)
            cut_option_y.append(layer * 10)
            cut_option_x.append(position * 10 - (layer_length * 10)/2)
            cut_info.append(cut_dag.layers.get(layer)[position].cuts)
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
                                         color='#d3d3d3'
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
                             textfont_size=30
                             ))

    fig.show()

def generate_cut_dag_main():
    freeze_support()
    cut_dag = make_cut_dag()
    visualizer(cut_dag)
