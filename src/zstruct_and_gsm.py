from os import makedirs, path, listdir
from shutil import move, copyfile, copytree, rmtree
from subprocess import check_call, DEVNULL, STDOUT, CalledProcessError
from uuid import uuid4
from re import sub
from rdkit.Chem import MolFromSmiles, MolToXYZBlock, AddHs
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.rdDistGeom import EmbedMolecule
from src.stringfile_to_rdkit import stringfile_to_rdkit


def run_zstruct_and_gsm(xyz_strings: list, smiles_string: str, ordering=None, core=None, reaction_folder: str = None, cuts_folder: str = ""):
    """takes an xyz string, a dictionary mapping the order of atoms, a list of core atoms and the isomer string, in return it creates output in blackbox/output"""
    if core is None:
        core = []
    if ordering is None:
        ordering = {}
    # if smiles_string is an existing folder make it the output folder otherwise make a folder and make it the output folder
    if not path.isdir(smiles_string):
        output_folder = "blackbox/output/" + smiles_string + "_" + str(uuid4().hex[:4])  # unique identifier for output of this smiles string
        makedirs(output_folder)
    else:
        output_folder = smiles_string
    clone_name = str(uuid4().hex)  # unique identifier for zstruct and gsm folders of this process
    offset = 0
    if reaction_folder is None:
        print("run zstruct")
        prepare_zstruct(clone_name, xyz_strings, ordering, core)  # make zstruct clone
        offset += run_zstruct(clone_name, output_folder, offset)  # run zstruct clone
        isomers_str = None
    else:
        makedirs(f"{output_folder}/{reaction_folder}{cuts_folder}", exist_ok=True)
        for file in listdir(f"{output_folder}/{reaction_folder}"):
            if file.startswith("ISOMERS"):
                with open(f"{output_folder}/{reaction_folder}/{file}", "r") as f:
                    isomers_str = f.read()
                break
        isomers_str = sub(r'\d+', lambda m: ordering.get(m.group(), m.group()), isomers_str)  # replace numbers with dict mapping
        with open(f"{output_folder}/{reaction_folder}{cuts_folder}/ISOMERS0000", "w") as f:
            f.write(isomers_str)
        with open(f"{output_folder}/{reaction_folder}{cuts_folder}/initial0000.xyz", "w") as f:
            f.write(xyz_strings[0])
        offset += 1
    if path.isdir(f"blackbox/zstruct_clones/{clone_name}"):
        rmtree(f"blackbox/zstruct_clones/{clone_name}", ignore_errors=True)         # remove zstruct clone
    run_gsm(clone_name, output_folder, offset, isomers_str, reaction_folder, cuts_folder)                         # make gsm clone and run gsm clone
    if path.isdir(f"blackbox/gsm_clones/{clone_name}"):
        rmtree(f"blackbox/gsm_clones/{clone_name}", ignore_errors=True)             # remove gsm clone
    if reaction_folder is None:
        return output_folder
    elif path.isfile(f"{output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz0000"):
        return f"{output_folder}/{reaction_folder}{cuts_folder}stringfile.xyz0000"
    else:
        return "NO REACTION"


def run_zstruct(clone_name: str, output_folder: str, offset: int):
    try:
        check_call(["./zstruct.exe"], stdout=DEVNULL, stderr=STDOUT, cwd=f"blackbox/zstruct_clones/{clone_name}")  # run zstruct.exe in silent mode
    except CalledProcessError as e:
        # print(e.output)
        pass
    isomer_count = sum(filename.startswith("ISOMER") for filename in listdir(f"blackbox/zstruct_clones/{clone_name}/scratch"))  # find number of ISOMER files created
    for i in range(isomer_count):                                                                                            # move all ISOMER and initial files to output folder
        strid = str(i).zfill(4)
        makedirs(f"{output_folder}/reaction{strid}")
        move(f"blackbox/zstruct_clones/{clone_name}/scratch/ISOMERS{strid}", f"{output_folder}/reaction{strid}/ISOMERS{str(i+offset).zfill(4)}")
        move(f"blackbox/zstruct_clones/{clone_name}/scratch/initial{strid}.xyz", f"{output_folder}/reaction{strid}/initial{str(i+offset).zfill(4)}.xyz")
    return isomer_count


def prepare_zstruct(clone_name: str, xyz_strs: list, ordering: dict, core: list):
    copytree("blackbox/zstruct_clones/original", f"blackbox/zstruct_clones/{clone_name}")  # create clone of zstruct
    for i, str in enumerate(xyz_strs):
        with open(f"blackbox/zstruct_clones/{clone_name}/react{1+i}.xyz", "w") as f:              # create react file
            f.write(str)
        with open(f"blackbox/zstruct_clones/{clone_name}/frozen{1}.xyz", "a") as f:             # create frozen file
            for element in core:
                f.write(str(ordering.get(element)) + "\n")


def run_gsm_round(clone_name: str, output_folder: str, i: int, isomers_str: str, reaction_folder: str = None, cuts_folder: str = ""):
    if reaction_folder is None:
        ID = str(i).zfill(4)
        reaction_folder = f"reaction{ID}"
        reaction_and_cut = reaction_folder
    else:
        ID = reaction_folder[-4:]
        reaction_and_cut = reaction_folder + cuts_folder[0:-1]
    init_fn = f"initial{ID}.xyz"
    iso_fn = f"ISOMERS{ID}"
    with open(f"{output_folder}/{reaction_folder}/{iso_fn}", "r") as f:
        new_isomers_str = f.read()
    if isomers_str is None or new_isomers_str == isomers_str:               # only run gsm on reaction matching pattern
        copyfile(f"{output_folder}/{reaction_and_cut}/initial0000.xyz", f"blackbox/gsm_clones/{clone_name}/scratch/initial0000.xyz")
        copyfile(f"{output_folder}/{reaction_and_cut}/ISOMERS0000", f"blackbox/gsm_clones/{clone_name}/scratch/ISOMERS0000")
        try:
            check_call(["./gsm.orca"], cwd=f"blackbox/gsm_clones/{clone_name}")   # run gsm.orca in silent mode
        except CalledProcessError as e:
            #print(e.output)
            pass
        # find stringfile if one was made and move to output
        if path.exists(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000"):
            move(f"blackbox/gsm_clones/{clone_name}/stringfile.xyz0000",
                 f"{output_folder}/{reaction_folder}{cuts_folder}/stringfile.xyz{str(i).zfill(4)}")
        elif path.exists(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g"):
            move(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g",
                 f"{output_folder}/{reaction_folder}{cuts_folder}/stringfile.xyz{str(i).zfill(4)}")
        elif path.exists(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g1"):
            move(f"blackbox/gsm_clones/{clone_name}/scratch/stringfile.xyz0000g1",
                 f"{output_folder}/{reaction_folder}{cuts_folder}/stringfile.xyz{str(i).zfill(4)}")


def run_gsm(clone_name: str, output_folder: str, isomer_count: int, isomers_str: str, reaction_folder: str = None, cuts_folder: str = ""):
    copytree("blackbox/gsm_clones/original", f"blackbox/gsm_clones/{clone_name}")   # create clone of gsm
    for isomer_id in range(isomer_count):                                           # iterate over all isomers/initial pairs
        run_gsm_round(clone_name, output_folder, isomer_id, isomers_str, reaction_folder, cuts_folder)   # compute stringfile for pair


def zstruct_gsm_main():
    #smiles_string = 'CN=C([O-])N(C)C(=O)OC(C)=O'
    smiles_string = "CCO"
    mol = MolFromSmiles(smiles_string)
    mol = AddHs(mol)
    Compute2DCoords(mol)  # generate 2d coordinates
    EmbedMolecule(mol, randomSeed=0xf00d)  # generate 3d coordinates
    xyz_str_1 = MolToXYZBlock(mol)
    folders = run_zstruct_and_gsm([xyz_str_1], smiles_string)
    #run_zstruct_and_gsm(xyz_strings=[xyz_str_1], smiles_string="CCO_518e", ordering={}, core=[], reaction_folder="reaction0001", cuts_folder="/1_2_3/")
    #print(folders)
    #folders = ["blackbox/output/ade0008ff58c47a59cc34cc464041810\stringfiles/stringfile.xyz0009"]
