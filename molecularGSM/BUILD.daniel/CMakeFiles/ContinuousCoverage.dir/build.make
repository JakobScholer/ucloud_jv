# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /work/DanielStuff/intel/oneapi/intelpython/python3.7/envs/nikolai/bin/cmake

# The command to remove a file.
RM = /work/DanielStuff/intel/oneapi/intelpython/python3.7/envs/nikolai/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work/DanielStuff/nikolai/molecularGSM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/DanielStuff/nikolai/molecularGSM/BUILD

# Utility rule file for ContinuousCoverage.

# Include any custom commands dependencies for this target.
include CMakeFiles/ContinuousCoverage.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ContinuousCoverage.dir/progress.make

CMakeFiles/ContinuousCoverage:
	/work/DanielStuff/intel/oneapi/intelpython/python3.7/envs/nikolai/bin/ctest -D ContinuousCoverage

ContinuousCoverage: CMakeFiles/ContinuousCoverage
ContinuousCoverage: CMakeFiles/ContinuousCoverage.dir/build.make
.PHONY : ContinuousCoverage

# Rule to build all files generated by this target.
CMakeFiles/ContinuousCoverage.dir/build: ContinuousCoverage
.PHONY : CMakeFiles/ContinuousCoverage.dir/build

CMakeFiles/ContinuousCoverage.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ContinuousCoverage.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ContinuousCoverage.dir/clean

CMakeFiles/ContinuousCoverage.dir/depend:
	cd /work/DanielStuff/nikolai/molecularGSM/BUILD && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/DanielStuff/nikolai/molecularGSM /work/DanielStuff/nikolai/molecularGSM /work/DanielStuff/nikolai/molecularGSM/BUILD /work/DanielStuff/nikolai/molecularGSM/BUILD /work/DanielStuff/nikolai/molecularGSM/BUILD/CMakeFiles/ContinuousCoverage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ContinuousCoverage.dir/depend

