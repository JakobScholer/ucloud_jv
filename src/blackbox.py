from multiprocessing import Queue
from os import makedirs, listdir
from os.path import isdir, isfile, exists
from shutil import move, copyfile, copytree, rmtree
from subprocess import check_call, DEVNULL, STDOUT, CalledProcessError, TimeoutExpired
from uuid import uuid4
from re import sub
from src.stringfile_tester import check_initial_file

def run_zstruct(name: str, xyz_strings: list, core=None, ordering=None, debug=False):
    if core is None:
        core = []
    if ordering is None:
        ordering = {}
    isomer_count = 0

    output_folder = f"blackbox/output/{name}_{str(uuid4().hex[:4])}"  # unique identifier for output of this smiles string
    makedirs(output_folder)

    clone_name = str(uuid4().hex)
    prepare_zstruct(clone_name, xyz_strings, ordering, core)  # make zstruct clone
    isomer_count += run_zstruct_computation(clone_name, output_folder, multiple_molecules=len(xyz_strings)-1, logfile=debug)  # run zstruct clone
    if isdir(f"blackbox/zstruct_clones/{clone_name}") and not debug:
        rmtree(f"blackbox/zstruct_clones/{clone_name}", ignore_errors=True)  # remove zstruct clone
    return output_folder, isomer_count


def prepare_zstruct(clone_name: str, xyz_strs: list, ordering: dict, core: list):
    copytree("blackbox/zstruct_clones/original", f"blackbox/zstruct_clones/{clone_name}")  # create clone of zstruct
    for i, xyz_str in enumerate(xyz_strs):
        with open(f"blackbox/zstruct_clones/{clone_name}/react{1+i}.xyz", "w") as f:              # create react file
            f.write(xyz_str)
        with open(f"blackbox/zstruct_clones/{clone_name}/frozen{1+i}.xyz", "a") as f:             # create frozen file
            for element in core:
                f.write(str(ordering.get(element, element)) + "\n")


def run_zstruct_computation(clone_name: str, output_folder: str, multiple_molecules: int, logfile: bool):
    try:
        # write to logfile or discard output
        if logfile:
            f_out = open(f"blackbox/zstruct_clones/{clone_name}/log_out", "w+")
            f_err = open(f"blackbox/zstruct_clones/{clone_name}/log_err", "w+")
        else:
            f_out = DEVNULL
            f_err = STDOUT
        check_call([f"./zstruct.exe -250 1000 1 {multiple_molecules}"], stdout=f_out, stderr=f_err, shell=True, cwd=f"blackbox/zstruct_clones/{clone_name}")  # run zstruct.exe in silent mode
    except CalledProcessError as e:
        # print(e.output)
        pass
    isomer_count = sum(filename.startswith("ISOMER") for filename in listdir(f"blackbox/zstruct_clones/{clone_name}/scratch"))  # find number of ISOMER files created
    for i in range(isomer_count): # move all ISOMER and initial files to output folder
        str_id = str(i).zfill(4)
        if multiple_molecules == 1: # if more than one molecule, check if reactions arre the same
            initial_file = f"blackbox/zstruct_clones/{clone_name}/scratch/initial{str_id}.xyz"
            if check_initial_file(initial_file): # check if the molecules are the same
                isomer_count -= 1
                continue # skip to next initial and isomer file
        makedirs(f"{output_folder}/reaction{str_id}")
        move(f"blackbox/zstruct_clones/{clone_name}/scratch/ISOMERS{str_id}", f"{output_folder}/reaction{str_id}/ISOMERS0000")
        move(f"blackbox/zstruct_clones/{clone_name}/scratch/initial{str_id}.xyz", f"{output_folder}/reaction{str_id}/initial0000.xyz")
    return isomer_count


def run_gsm_initial_multi_threaded(task_queue: Queue, output_folder: str, isomer_count: int, logfile: bool = False):
    for reaction_folder in listdir(output_folder):   # iterate over all isomers/initial pairs
        task_queue.put((run_gsm_round_initial, (output_folder, reaction_folder, logfile)))  ## compute stringfile for pair


def run_gsm_round_initial(output_folder: str, reaction_folder: str, logfile: bool = False):
    if not isdir(output_folder):
        print(f"{output_folder} not made yet")
        return "ERROR"
    if isfile(f"{output_folder}/{reaction_folder}/stringfile.xyz"):
        print(f"skipping {output_folder}/{reaction_folder}/stringfile.xyz")
        return
    print(f"working on {output_folder}/{reaction_folder}/stringfile.xyz")
    clone_name = str(uuid4().hex)
    copytree("blackbox/gsm_clones/original", f"blackbox/gsm_clones/{clone_name}")  # create clone of gsm
    copyfile(f"{output_folder}/{reaction_folder}/initial0000.xyz", f"blackbox/gsm_clones/{clone_name}/scratch/initial0000.xyz")
    copyfile(f"{output_folder}/{reaction_folder}/ISOMERS0000", f"blackbox/gsm_clones/{clone_name}/scratch/ISOMERS0000")
    try:
        # write to logfile or discard output
        if logfile:
            f_out = open(f"blackbox/gsm_clones/{clone_name}/log_out", "w+")
            f_err = open(f"blackbox/gsm_clones/{clone_name}/log_err", "w+")
        else:
            f_out = DEVNULL
            f_err = STDOUT
        check_call(["./gsm.orca"], stdout=f_out, stderr=f_err, timeout=600, cwd=f"blackbox/gsm_clones/{clone_name}")   # run gsm.orca in silent mode
    except CalledProcessError as e:
        #print(e.output)
        pass
    except TimeoutExpired as e:
        print(f"    Timeout for {output_folder}/{reaction_folder}/stringfile.xyz")
        print(e)
    # find stringfile if one was made and move to output
    if exists(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000"):
        move(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000",
             f"{output_folder}/{reaction_folder}/stringfile.xyz")
    elif exists(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g"):
        move(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g",
             f"{output_folder}/{reaction_folder}/stringfile.xyz")
    elif exists(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g1"):
        move(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g1",
             f"{output_folder}/{reaction_folder}/stringfile.xyz")
    else:
        print(f"    no stringfile generated for {output_folder}/{reaction_folder}/stringfile.xyz")
    if isdir(f"blackbox/gsm_clones/{clone_name}") and not logfile:
        rmtree(f"blackbox/gsm_clones/{clone_name}", ignore_errors=True)  # remove gsm clone
    print(f"    Completed work on {output_folder}/{reaction_folder}/stringfile.xyz")


def run_gsm_round_cuts(clone_name: str, output_folder: str, reaction_folder: str, cuts_folder: str, isomers_str: str, logfile: bool = False):
    reaction_and_cut = reaction_folder + cuts_folder[0:-1]
    with open(f"{output_folder}/{reaction_folder}/ISOMERS0000", "r") as f:
        new_isomers_str = f.read()
    if isomers_str is None or new_isomers_str == isomers_str:               # only run gsm on reaction matching pattern
        copyfile(f"{output_folder}/{reaction_and_cut}/initial0000.xyz", f"blackbox/gsm_clones/{clone_name}/scratch/initial0000.xyz")
        copyfile(f"{output_folder}/{reaction_and_cut}/ISOMERS0000", f"blackbox/gsm_clones/{clone_name}/scratch/ISOMERS0000")
        try:
            # write to logfile or discard output
            if logfile:
                f_out = open(f"blackbox/gsm_clones/{clone_name}/log_out", "w+")
                f_err = open(f"blackbox/gsm_clones/{clone_name}/log_err", "w+")
            else:
                f_out = DEVNULL
                f_err = STDOUT
            check_call(["./gsm.orca"], stdout=f_out, stderr=f_err, timeout=600, cwd=f"blackbox/gsm_clones/{clone_name}")   # run gsm.orca in silent mode
        except CalledProcessError as e:
            #print(e.output)
            pass
        except TimeoutExpired as e:
            print(f"    Timeout for {output_folder}/{reaction_folder}{cuts_folder}/stringfile.xyz")
            print(e)
        # find stringfile if one was made and move to output
        if exists(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000"):
            move(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000",
                 f"{output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz")
        elif exists(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g"):
            move(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g",
                 f"{output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz")
        elif exists(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g1"):
            move(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g1",
                 f"{output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz")
        print(f"    Completed work on {output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz")


def run_gsm_cuts(xyz_strings: list, output_folder: str, reaction_folder: str, cuts_folder: str = "", ordering=None, logfile: bool = False):
    print(f"working on {output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz")
    if not isdir(output_folder):
        print(f"{output_folder} not made yet")
        return "ERROR"
    isomer_count = 1

    makedirs(f"{output_folder}/{reaction_folder}{cuts_folder}", exist_ok=True)
    for file in listdir(f"{output_folder}/{reaction_folder}"):
        if file.startswith("ISOMERS"):
            with open(f"{output_folder}/{reaction_folder}/{file}", "r") as f:
                isomers_str = f.read()
            break
    def replacement(m):
        return str(ordering.get(int(m.group(0)), int(m.group(0))))
    modified_isomers_str = sub(r'\d+', replacement, isomers_str)  # replace numbers with dict mapping
    with open(f"{output_folder}/{reaction_folder}{cuts_folder}/ISOMERS0000", "w") as f:
        f.write(modified_isomers_str)
    with open(f"{output_folder}/{reaction_folder}{cuts_folder}/initial0000.xyz", "w") as f:
        f.write(xyz_strings[0])

    for isomer_id in range(isomer_count):        # iterate over all isomers/initial pairs
        clone_name = str(uuid4().hex)
        copytree("blackbox/gsm_clones/original", f"blackbox/gsm_clones/{clone_name}")  # create clone of gsm
        run_gsm_round_cuts(clone_name, output_folder, reaction_folder, cuts_folder, isomers_str)   # compute stringfile for pair
        if isdir(f"blackbox/gsm_clones/{clone_name}") and not logfile:
            rmtree(f"blackbox/gsm_clones/{clone_name}", ignore_errors=True)  # remove gsm clone
    if isfile(f"{output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz"):
        return f"{output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz"
    else:
        return "EMPTY"
