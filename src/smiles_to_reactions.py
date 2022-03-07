from src.generate_cut_dag import make_cut_dag, show_cut_dag
from src.stringfile_helper_functions import max_energy_curve
from src.blackbox import run_zstruct, run_gsm_initial_multi_threaded
from src.stringfile_tester import check_educt_to_product

from rdkit.Chem import RWMol, AddHs, MolFromSmiles, MolToXYZBlock, rdDepictor, GetPeriodicTable
from rdkit.Chem.AllChem import EmbedMolecule
from os import listdir
from glob import glob
from multiprocessing import freeze_support, Queue, Process
from time import sleep

'''
takes a mode and a string or list of strings as input
    blackbox True  -  generates stringfiles, either from scratch using list of smiles strings as input
or from existing folder of reactions using name of folder
    blackbox False -  reads data from a folder. string_data must be the path to the already compiled cut dag data
'''
def make_reactions(blackbox: bool, string_data, max_energy: int=100, frozen=None, visual_cut_dag: bool=False, visual_stringfiles: bool=False, number_of_processes: int=1, debug: bool=False):
    if frozen is None:
        frozen = []

    # set up multi threading
    freeze_support()
    task_queue = Queue()
    done_queue = Queue()
    for i in range(number_of_processes):
        Process(target=worker, args=(task_queue, done_queue)).start()

    if blackbox:
        # input is list of smiles strings (runs zstruct to find reactions for them)
        if isinstance(string_data, list):
            xyz_list = []
            # convert smiles strings to xyz files
            for string in string_data:
                mol = RWMol(MolFromSmiles(string))      # make rdkit mol object from smiles string
                mol = AddHs(mol, explicitOnly=False)    # add hydrogen for good measure.
                rdDepictor.Compute2DCoords(mol)         # add 2D-coordinates to mol object
                EmbedMolecule(mol, randomSeed=0xf00d)   # convert to 3D-coordinates
                xyz_list.append(MolToXYZBlock(mol))     # convert to xyz file

            # set folder name based on input smiles strings
            reaction_name = string_data[0]
            for i in range(1, len(string_data)):
                reaction_name = reaction_name + "_+_" + string_data[i]
            # run zstruct
            smiles_path, isomer_count = run_zstruct(reaction_name, xyz_list, core=frozen, debug=False)
        else:
            # input is name of folder containing reactions
            smiles_path = f"blackbox/output/{string_data}"
            isomer_count = len(listdir(smiles_path))
        # run GSM
        run_gsm_initial_multi_threaded(task_queue, smiles_path, isomer_count)
        # wait for all gsm runs to finish
        while not done_queue.qsize() == isomer_count:
            sleep(5)
        # empty done_queue
        while not done_queue.empty():
            done_queue.get()
    else:
        smiles_path = string_data

    stringfile_path = listdir(smiles_path)
    reaction_folders = [smiles_path + "/" + s for s in stringfile_path]
    reaction_folders.sort()
    print("-------------------------------------------------------------------")
    # make cut dag for every reaction folder containing a stringfile
    stringfile_list = []
    total_tasks = 0
    for folder in reaction_folders:
        stringfiles = glob(f"{folder}/stringfile*")
        if stringfiles: # check if stringfile exists
            check_ep = check_educt_to_product(stringfiles[0])
            max_ec = max_energy_curve(stringfiles[0], max_energy)
            if check_ep and max_ec: # if there is a reaction in the stringfile. make a cut dag!
                stringfile_list.append(stringfiles[0])
                assigned_tasks = make_cut_dag(task_queue, stringfiles[0], debug)
                total_tasks += assigned_tasks
            else:
                log_data = "0\n"
                log_data += "Cut dag not generated\n"
                log_data += f"Educt to product: {check_ep}\n"
                log_data += f"Max energy curve: {max_ec}\n"
                with open(f"{folder}/done.txt", "w") as f:
                    f.write(log_data)
    while not done_queue.qsize() == total_tasks:
        sleep(5)
    # Tell child processes to stop
    for i in range(number_of_processes):
        task_queue.put('STOP')
    # show cut dags
    for stringfile in stringfile_list:
        show_cut_dag(stringfile, visual_cut_dag, visual_stringfiles, debug)


def make_single_reaction(stringfile: str, number_of_processes: int=1, debug: bool=False):
    # set up multi threading
    freeze_support()
    task_queue = Queue()
    done_queue = Queue()
    for i in range(number_of_processes):
        Process(target=worker, args=(task_queue, done_queue)).start()

    total_tasks = make_cut_dag(task_queue, stringfile, debug)
    while not done_queue.qsize() == total_tasks:
        sleep(5)
    # Tell child processes to stop
    for i in range(number_of_processes):
        task_queue.put('STOP')


# Function run by worker processes
def worker(i, o):
    for func, args in iter(i.get, 'STOP'):
        func(*args)
        o.put('DONE')
