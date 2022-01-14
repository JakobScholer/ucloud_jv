import unittest
from src.blackbox import prepare_zstruct, run_zstruct
from os import makedirs, listdir
from os.path import isdir, isfile
from uuid import uuid4
from shutil import rmtree


class Test(unittest.TestCase):
    # default behaviour
    def test_prepare_zstruct(self):
        # preparing zstruct
        clone_name = f"test_{str(uuid4().hex)}"
        with open("test/testfiles/stringfile.xyz0000") as f:
            xyz_string = f.read()
        prepare_zstruct(clone_name=clone_name, xyz_strs=[xyz_string], ordering={}, core=[])
        # checking result
        result = False
        msg = "zstruct not prepared correctly"
        if isdir(f"blackbox/zstruct_clones/{clone_name}"):
            msg = "zstruct executable missing"
            if isfile(f"blackbox/zstruct_clones/{clone_name}/zstruct.exe"):
                result = True
        self.assertTrue(result, msg)
        # cleanup
        rmtree(f"blackbox/zstruct_clones/{clone_name}", ignore_errors=True)  # remove zstruct clone

    # default behaviour
    def test_run_zstruct(self):
        # preparing zstruct
        clone_name = f"test_{str(uuid4().hex)}"
        with open("test/testfiles/stringfile.xyz0000") as f:
            xyz_string = f.read()
        prepare_zstruct(clone_name=clone_name, xyz_strs=[xyz_string], ordering={}, core=[])
        # running zstruct
        output_folder = f"test/blackbox/{clone_name}"
        makedirs(output_folder)
        run_zstruct(clone_name=clone_name, output_folder=output_folder, logfile=False)
        # checking number of reactions
        number_of_reactions = len(listdir(f"test/blackbox/{clone_name}"))
        self.assertEqual(155, number_of_reactions)
        # checking content of initial file
        with open("test/testfiles/stringfile.xyz0000") as f:
            expected_initial_file = [next(f) for x in range(22)]
        with open(f"test/blackbox/{clone_name}/reaction0000/initial0000.xyz") as f:
            actual_initial_file = [next(f) for x in range(22)]
        self.assertAlmostEqual(expected_initial_file, actual_initial_file, 5, "hi")
        # cleanup
        rmtree(f"blackbox/zstruct_clones/{clone_name}", ignore_errors=True)  # remove zstruct clone
        #rmtree(f"test/blackbox/{clone_name}", ignore_errors=True)  # remove output folder

if __name__ == '__main__':
    unittest.main()
