from os import system, chdir, makedirs, walk, path, listdir
from shutil import move, copyfile, copytree, rmtree
from uuid import uuid4
from re import sub


def run_zstruct_and_gsm(xyz_string: str, ordering=None, core=None, isomers_str=None):
    """takes an xyz string, a dictionary mapping the order of atoms, a list of core atoms and the isomer string, in return it creates output in blackbox/output"""
    if ordering is None:    # ordering is None for the first run of a compound, otherwise a dict
        ordering = {}
    if core is None:        # core is None for the first run of a compound, otherwise a list
        core = []
    if isomers_str is not None:
        isomers_str = sub(r'\d+', lambda m: ordering.get(m.group(), m.group()), isomers_str)    # replace numbers with dict mapping
    clone_name = str(uuid4().hex)  # unique identifier for folders of this process
    offset = 0
    makedirs(f"blackbox/output/{clone_name}/isomers")
    makedirs(f"blackbox/output/{clone_name}/initial")
    makedirs(f"blackbox/output/{clone_name}/stringfiles")

    prepare_zstruct(clone_name, xyz_string, ordering, core) # make zstruct clone
    offset += run_zstruct(clone_name, offset)               # run zstruct clone
    rmtree(f"blackbox/zstruct_clones/{clone_name}")         # remove zstruct clone
    run_gsm(clone_name, offset, isomers_str)                # make gsm clone and run gsm clone
    rmtree(f"blackbox/gsm_clones/{clone_name}")             # remove gsm clone

    stringfile_path = listdir(f"blackbox/output/{clone_name}/stringfiles")  # returns list of all generated stringfiles
    return [f"blackbox/output/{clone_name}/stringfiles/" + s for s in stringfile_path]


def run_zstruct(clone_name: str, offset: int):
    chdir(f"blackbox/zstruct_clones/{clone_name}")  # set current folder to zstruct.exe folder
    system("./zstruct.exe")                         # run zstruct.exe
    chdir("../../..")                               # reset current folder to ucloud_jv
    _, _, filenames = next(walk(f"blackbox/zstruct_clones/{clone_name}/scratch"))
    isomer_count = sum(fn.startswith("ISO") for fn in filenames)

    for i in range(isomer_count):
        strid = str(i).zfill(4)
        move(f"blackbox/zstruct_clones/{clone_name}/scratch/ISOMERS{strid}", f"blackbox/output/{clone_name}/isomers/ISOMERS{str(i+offset).zfill(4)}")
        move(f"blackbox/zstruct_clones/{clone_name}/scratch/initial{strid}.xyz", f"blackbox/output/{clone_name}/initial/initial{str(i+offset).zfill(4)}.xyz")
    return isomer_count


def prepare_zstruct(clone_name: str, xyz_str: str, ordering: dict, core: list):
    copytree("blackbox/zstruct_clones/original", f"blackbox/zstruct_clones/{clone_name}")  # create clone of zstruct
    with open(f"blackbox/zstruct_clones/{clone_name}/react{1}.xyz", "w") as f:              # create react file
        f.write(xyz_str)
    with open(f"blackbox/zstruct_clones/{clone_name}/frozen{1}.xyz", "a") as f:             # create frozen file
        for element in core:
            f.write(str(ordering.get(element)) + "\n")


def run_gsm_round(clone_name, i: int, isomers_str: str):
    ID = str(i).zfill(4)
    init_fn = f"initial{ID}.xyz"
    iso_fn = f"ISOMERS{ID}"
    with open(f"blackbox/output/{clone_name}/initial/{init_fn}", "r") as f:
        new_isomers_str = f.read()
    if isomers_str is None or new_isomers_str == isomers_str:               # only run gsm on reaction matching pattern
        copyfile(f"blackbox/output/{clone_name}/initial/{init_fn}", f"blackbox/gsm_clones/{clone_name}/scratch/initial0000.xyz")
        copyfile(f"blackbox/output/{clone_name}/isomers/{iso_fn}", f"blackbox/gsm_clones/{clone_name}/scratch/ISOMERS0000")

        chdir(f"blackbox/gsm_clones/{clone_name}")  # set current folder to gsm.orca folder
        system("./gsm.orca")                        # run gsm.orca
        chdir("../../..")                           # reset current folder to ucloud_jv


def run_gsm(clone_name: str, isomer_count: int, isomers_str: str):
    copytree("blackbox/gsm_clones/original", f"blackbox/gsm_clones/{clone_name}")   # create clone of gsm
    for isomer_id in range(isomer_count):                                           # iterate over all isomers/initial pairs
        run_gsm_round(clone_name, isomer_id, isomers_str)                           # compute stringfile for pair
        if path.exists(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000"):     # move stringfile to output if it was made
            move(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000",
                 f"blackbox/output/{clone_name}/stringfiles/stringfile.xyz{str(isomer_id).zfill(4)}")
