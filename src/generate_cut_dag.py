import time
import random
from src.cut_dag import CutDagNode, CutDag, make_childs_mp, insert_childs_mp, make_root
from multiprocessing import Process, Queue, current_process, freeze_support
from igraph import *
import plotly.graph_objects as go

#
# Function run by worker processes
#

def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = func(*args)
        output.put(result)

def make_cut_dag():
    NUMBER_OF_PROCESSES = 4

    # make test cut dag
    stringfile = 'test/testfiles/stringfile_ring.xyz0000'
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
            print("whuuee got one bunch of childs")
            child_infos = done_queue.get() # get info
            tasks = insert_childs_mp(stringfile, cd, child_infos[0], child_infos[1]) # insert child
            if len(tasks) > 0:
                print("added " + str(len(tasks)) + " to the tasks")
                for t in tasks: # add new tasks to the queue
                    task_queue.put(t)
            else:
                print("no new tasks added")
    print("done waiting for it")

    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')
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
