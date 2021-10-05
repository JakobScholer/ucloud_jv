import sys
from src.runner import runner_main
from src.generate_tree import generate_tree_main

if __name__ == '__main__':
    print(sys.argv)
    if len(sys.argv) == 2:
        print(sys.argv)
        if str(sys.argv[1]) == "runner":
            runner_main()
        elif str(sys.argv[1]) == "generate_tree":
            generate_tree_main()
        elif str(sys.argv[1]) == "mod_to_xyz":
            runner_main()
        elif str(sys.argv[1]) == "xyz_to_mod":
            runner_main()
    else:
        print("possible arguments are: 'runner', 'generate_tree', 'mod_to_xyz', 'xyz_to_mod'")
