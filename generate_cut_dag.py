import time
import random
from multiprocessing import Process, Queue, current_process, freeze_support, Lock
from src.cut_dag import CutDagNode, CutDag, make_childs, make_root

#
# Function run by worker processes
#

def worker(input, lock):
    for func, args in iter(input.get, 'STOP'):
        result = do_task(func, args, input, lock)
        print(result)

#
# Function used to calculate result
#

def do_task(func, args, input, lock):
    result = func(*args)
    for r in result:
        new_task = (make_childs, (r, cd, lock))
        input.put(new_task)
    return '%s says that %s%s = %s' % \
        (current_process().name, func.__name__, args, result)

def test():
    NUMBER_OF_PROCESSES = 4
    cd_lock = Lock()

    # make test cut dag
    test = 'test/testfiles/stringfile_ring.xyz0000'
    graph = False
    cd = make_root(test, graph)

    # Create queues
    task_queue = Queue()

    # add first task
    root_task = (make_childs, (cd.layers[0][0], cd, cd_lock))
    task_queue.put(root_task)

    # Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue, cd_lock)).start()

    # wait for porcesses to end
    wait_for_end = True
    while wait_for_end:
        # make sure the quee is empty before stopping porcesses
        if task_queue.empty():
            time.sleep(5)
            if task_queue.empty():
                # Tell child processes to stop
                for i in range(NUMBER_OF_PROCESSES):
                    task_queue.put('STOP')
                wait_for_end = False


if __name__ == '__main__':
    freeze_support()
    test()
