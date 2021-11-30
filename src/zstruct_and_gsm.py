import subprocess
from os import makedirs, walk, path, listdir
from shutil import move, copyfile, copytree, rmtree
from subprocess import check_call, DEVNULL, STDOUT
from uuid import uuid4
from re import sub

from rdkit.Chem import MolFromSmiles, MolToXYZBlock, AddHs
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.rdDistGeom import EmbedMolecule

from src.stringfile_to_rdkit import stringfile_to_rdkit


def run_zstruct_and_gsm(xyz_strings: list, ordering=None, core=None, isomers_str=None):
    """takes an xyz string, a dictionary mapping the order of atoms, a list of core atoms and the isomer string, in return it creates output in blackbox/output"""
    clone_name = str(uuid4().hex)  # unique identifier for folders of this process
    makedirs(f"blackbox/output/{clone_name}/isomers")
    makedirs(f"blackbox/output/{clone_name}/initial")
    makedirs(f"blackbox/output/{clone_name}/stringfiles")
    offset = 0

    if ordering is None:    # ordering is None for the first run of a compound, otherwise a dict
        ordering = {}
    if core is None:        # core is None for the first run of a compound, otherwise a list
        core = []
    if isomers_str is not None:
        isomers_str = sub(r'\d+', lambda m: ordering.get(m.group(), m.group()), isomers_str)    # replace numbers with dict mapping
        with open(f"blackbox/output/{clone_name}/isomers/ISOMERS{str(offset).zfill(4)}", "w") as f:
            f.write(isomers_str)
        with open(f"blackbox/output/{clone_name}/initial/initial{str(offset).zfill(4)}.xyz", "w") as f:
            f.write(xyz_strings[0])
        offset += 1
    else:
        prepare_zstruct(clone_name, xyz_strings, ordering, core)  # make zstruct clone
        offset += run_zstruct(clone_name, offset)  # run zstruct clone

    if path.isdir(f"blackbox/zstruct_clones/{clone_name}"):
        rmtree(f"blackbox/zstruct_clones/{clone_name}")         # remove zstruct clone
    run_gsm(clone_name, offset, isomers_str)                # make gsm clone and run gsm clone
    if path.isdir(f"blackbox/gsm_clones/{clone_name}"):
        rmtree(f"blackbox/gsm_clones/{clone_name}")             # remove gsm clone

    stringfile_path = listdir(f"blackbox/output/{clone_name}/stringfiles")  # returns list of all generated stringfiles
    return [f"blackbox/output/{clone_name}/stringfiles/" + s for s in stringfile_path]


def run_zstruct(clone_name: str, offset: int):
    try:
        check_call(["./zstruct.exe"], stdout=DEVNULL, stderr=STDOUT, cwd=f"blackbox/zstruct_clones/{clone_name}")  # run zstruct.exe in silent mode
    except subprocess.CalledProcessError as e:
        # print(e.output)
        pass
    _, _, filenames = next(walk(f"blackbox/zstruct_clones/{clone_name}/scratch"))
    isomer_count = sum(fn.startswith("ISO") for fn in filenames)

    for i in range(isomer_count):
        strid = str(i).zfill(4)
        move(f"blackbox/zstruct_clones/{clone_name}/scratch/ISOMERS{strid}", f"blackbox/output/{clone_name}/isomers/ISOMERS{str(i+offset).zfill(4)}")
        move(f"blackbox/zstruct_clones/{clone_name}/scratch/initial{strid}.xyz", f"blackbox/output/{clone_name}/initial/initial{str(i+offset).zfill(4)}.xyz")
    return isomer_count


def prepare_zstruct(clone_name: str, xyz_strs: list, ordering: dict, core: list):
    copytree("blackbox/zstruct_clones/original", f"blackbox/zstruct_clones/{clone_name}")  # create clone of zstruct
    for i, str in enumerate(xyz_strs):
        with open(f"blackbox/zstruct_clones/{clone_name}/react{1+i}.xyz", "w") as f:              # create react file
            f.write(str)
        #with open(f"blackbox/zstruct_clones/{clone_name}/frozen{1}.xyz", "a") as f:             # create frozen file
            #for element in core:
                #f.write(str(ordering.get(element)) + "\n")


def run_gsm_round(clone_name, i: int, isomers_str: str):
    ID = str(i).zfill(4)
    init_fn = f"initial{ID}.xyz"
    iso_fn = f"ISOMERS{ID}"
    with open(f"blackbox/output/{clone_name}/isomers/{iso_fn}", "r") as f:
        new_isomers_str = f.read()
    if isomers_str is None or new_isomers_str == isomers_str:               # only run gsm on reaction matching pattern
        copyfile(f"blackbox/output/{clone_name}/initial/{init_fn}", f"blackbox/gsm_clones/{clone_name}/scratch/initial0000.xyz")
        copyfile(f"blackbox/output/{clone_name}/isomers/{iso_fn}", f"blackbox/gsm_clones/{clone_name}/scratch/ISOMERS0000")
        try:
            check_call(["./gsm.orca"], stdout=DEVNULL, stderr=STDOUT, cwd=f"blackbox/gsm_clones/{clone_name}")   # run gsm.orca in silent mode
        except subprocess.CalledProcessError as e:
            #print(e.output)
            pass


def run_gsm(clone_name: str, isomer_count: int, isomers_str: str):
    copytree("blackbox/gsm_clones/original", f"blackbox/gsm_clones/{clone_name}")   # create clone of gsm
    for isomer_id in range(isomer_count):                                           # iterate over all isomers/initial pairs
        run_gsm_round(clone_name, isomer_id, isomers_str)                           # compute stringfile for pair
        # find stringfile if one was made and move to output
        if path.exists(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000"):
            move(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000",
                 f"blackbox/output/{clone_name}/stringfiles/stringfile.xyz{str(isomer_id).zfill(4)}")
        elif path.exists(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g"):
            move(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g",
                 f"blackbox/output/{clone_name}/stringfiles/stringfile.xyz{str(isomer_id).zfill(4)}")
        elif path.exists(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g1"):
            move(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g1",
                 f"blackbox/output/{clone_name}/stringfiles/stringfile.xyz{str(isomer_id).zfill(4)}")


def zstruct_gsm_main():
    mol = MolFromSmiles('CN=C([O-])N(C)C(=O)OC(C)=O')
    mol = AddHs(mol)
    Compute2DCoords(mol)  # generate 2d coordinates
    EmbedMolecule(mol, randomSeed=0xf00d)  # generate 3d coordinates
    xyz_str_1 = MolToXYZBlock(mol)
    #folders = run_zstruct_and_gsm([xyz_str_1])
    folders = ["blackbox/output/ade0008ff58c47a59cc34cc464041810\stringfiles/stringfile.xyz0009"]
    for file in folders:
        mol, atom_core, energy_profiles = stringfile_to_rdkit(file)
        if atom_core != set():
            print(file)
            print(atom_core)
        else:
            print(file)
            print("NO CORE")
