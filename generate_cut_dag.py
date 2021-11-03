import time
import random
from src.cut_dag import CutDagNode, CutDag, make_childs_mp, insert_childs_mp, make_root
from multiprocessing import Process, Queue, current_process, freeze_support

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
            while done_queue.empty() == False
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

    # find root og check op på den
    print("first test")
    print(cd.layers[0][0].stringfile)
    print(cd.layers[0][0].energy)
    print(cd.layers[0][0].cuts)
    # Lav børn på root og chek op på dem
    print("Second test")
    for node in cd.layers[1]:
        print(node.cuts)
    print("Second test")
    for node in cd.layers[2]:
        print(node.cuts)
    # lav børn på børnene :D det skal nok blive magisk


if __name__ == '__main__':
    freeze_support()
    make_cut_dag()
