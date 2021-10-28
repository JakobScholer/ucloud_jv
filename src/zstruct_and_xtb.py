import json
import os
import shutil


def run_zstruct(outdir, offset):
    os.chdir("zstruct")
    os.system("./zstruct.exe")
    _, _, filenames = next(os.walk("scratch"))
    isomer_count = sum(fn.startswith("ISO") for fn in filenames)
    os.chdir("..")

    for i in range(isomer_count):
        strid = str(i).zfill(4)
        shutil.move(f"zstruct/scratch/ISOMERS{strid}", f"{outdir}/ISOMERS{str(i+offset).zfill(4)}")
        shutil.move(f"zstruct/scratch/initial{strid}.xyz", f"{outdir}/initial{str(i+offset).zfill(4)}.xyz")

    return isomer_count


def prepare_zstruct(combination, molecules, dir_path):
    assert(len(combination) <= 2)
    os.system("rm -r zstruct/scratch")
    os.system("rm zstruct/*.xyz")
    os.mkdir('zstruct/scratch')

    for i, name in enumerate(combination):
        m = molecules[name]
        shutil.copy(f"{dir_path}/{m['xyz']}", f"zstruct/react{i+1}.xyz")
        shutil.copy(f"{dir_path}/{m['frozen']}", f"zstruct/frozen{i+1}.xyz")


def generate_isomers(path: str):
    out_dir = "scratch/isomers"
    if os.path.exists(out_dir):
        os.system(f"rm -r {out_dir}")
    os.mkdir(out_dir)

    with open(path) as f:
        data = json.load(f)

    dir_path = os.path.dirname(path)
    offset = 0
    for c in data["combinations"]:
        prepare_zstruct(c, data["molecules"], dir_path)
        offset += run_zstruct(out_dir, offset)
    return offset


def prepare_ssm(isomer_path, i: int):
    ssm_path = "ssm/scratch"
    if os.path.exists(ssm_path):
        os.system(f"rm -r {ssm_path}")
    os.makedirs(ssm_path, exist_ok=True)

    ID = str(i).zfill(4)
    init_fn = f"initial{ID}.xyz"
    iso_fn = f"ISOMERS{ID}"
    os.system(f"cp {isomer_path}/{init_fn} {ssm_path}/initial0000.xyz")
    os.system(f"cp {isomer_path}/{iso_fn} {ssm_path}/ISOMERS0000")


def execute_ssm():
    os.chdir("ssm")
    os.system("./gsm.orca")
    os.chdir("..")


def run_ssm(isomer_count: int, offset: int = 0):
    outdir = "scratch/stringfiles"
    if os.path.exists(outdir):
        os.system(f"rm -r {outdir}")
    os.makedirs(outdir, exist_ok=True)
    for i in range(isomer_count):
        isomer_id = i + offset
        prepare_ssm("scratch/isomers", isomer_id)
        execute_ssm()
        if os.path.exists("ssm/stringfile.xyz0000"):
            os.system(f"mv ssm/stringfile.xyz0000 {outdir}/stringfile.xyz{str(isomer_id).zfill(4)}")
        if os.path.exists("ssm/scratch/tsq0000.xyz"):
            os.system(f"mv ssm/scratch/tsq0000.xyz {outdir}/tsq.xyz{str(isomer_id).zfill(4)}")