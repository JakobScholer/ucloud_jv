import time
import random
from src.cut_dag import CutDagNode, CutDag, make_root
from multiprocessing import Process, Queue, current_process, freeze_support

#
# Function run by worker processes
#

def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = func(*args)
        output.put(result)

#
# Function used to calculate result
#

def calculate(func, args):
    result = func(*args)
    return '%s says that %s%s = %s' % \
        (current_process().name, func.__name__, args, result)

#
# Functions referenced by tasks
#

def mul(a, b, d):
    time.sleep(0.5*random.random())
    return a * b

def plus(a, b):
    time.sleep(0.5*random.random())
    return a + b

def make_childs(node, cutdag):
    time.sleep(0.5*random.random())
    return 1

#
#
#

def test():
    NUMBER_OF_PROCESSES = 4

    # make test cut dag
    test = 'test/testfiles/stringfile_ring.xyz0000'
    graph = False
    cd = make_root(test, graph)

    TASKS1 = [(mul, (i, 7, 9)) for i in range(19)]
    TASKS2 = [(plus, (i, 8)) for i in range(10)]

    # Create queues
    task_queue = Queue()
    done_queue = Queue()

    # Submit tasks
    for task in TASKS1:
        task_queue.put(task)
    print("make task for cd")
    # add first task
    root = cd.layers[0][0]
    root_task = (make_childs, (2, 1))
    #task_queue.put(root_task)
    print("done")

    # Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue, done_queue)).start()

    # Get and print results
    print('Unordered results:')
    for i in range(1):
        print('\t', done_queue.get())

    # wait for porcesses to end
    wait_for_end = True
    print("staring wait for it")
    while wait_for_end:
        # make sure the quee is empty before stopping porcesses
        if done_queue.empty():
            print("sleep time")
            time.sleep(2)
            print("waky waky!")
            if done_queue.empty():
                wait_for_end = False
        else:
            print('\t', done_queue.get())
    print("done waiting for it")

    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')


if __name__ == '__main__':
    freeze_support()
    test()
