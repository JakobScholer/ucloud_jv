import unittest
from unittest.mock import patch
from src.stringfile_to_rdkit import stringfile_to_rdkit


class Test(unittest.TestCase):
    # Normal behaviour test case
    def test_root_mean_square_normal(self):

        stringfile_to_rdkit("test/testfiles/stringfile.xyz0000")
        import sys
        sys.stdout.write(str(mock_print.call_args) + '\n')
        sys.stdout.write(str(mock_print.call_args_list) + '\n')


if __name__ == '__main__':
    unittest.main()
