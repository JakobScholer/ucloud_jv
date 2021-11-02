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

# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_SOURCE_DIR = /work/DanielStuff/nikolai/molecularGSM/cmake/tribits/core/config_tests/fmangle

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/DanielStuff/nikolai/molecularGSM/build-gcc/CMakeFiles/CMakeTmp/fmangle/LOWER/UNDER

# Include any dependencies generated for this target.
include CMakeFiles/cmangle.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cmangle.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cmangle.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cmangle.dir/flags.make

CMakeFiles/cmangle.dir/cmangle.c.o: CMakeFiles/cmangle.dir/flags.make
CMakeFiles/cmangle.dir/cmangle.c.o: /work/DanielStuff/nikolai/molecularGSM/cmake/tribits/core/config_tests/fmangle/cmangle.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/work/DanielStuff/nikolai/molecularGSM/build-gcc/CMakeFiles/CMakeTmp/fmangle/LOWER/UNDER/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/cmangle.dir/cmangle.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cmangle.dir/cmangle.c.o -c /work/DanielStuff/nikolai/molecularGSM/cmake/tribits/core/config_tests/fmangle/cmangle.c

CMakeFiles/cmangle.dir/cmangle.c.i: cmake_force
	@echo "Preprocessing C source to CMakeFiles/cmangle.dir/cmangle.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /work/DanielStuff/nikolai/molecularGSM/cmake/tribits/core/config_tests/fmangle/cmangle.c > CMakeFiles/cmangle.dir/cmangle.c.i

CMakeFiles/cmangle.dir/cmangle.c.s: cmake_force
	@echo "Compiling C source to assembly CMakeFiles/cmangle.dir/cmangle.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /work/DanielStuff/nikolai/molecularGSM/cmake/tribits/core/config_tests/fmangle/cmangle.c -o CMakeFiles/cmangle.dir/cmangle.c.s

# Object files for target cmangle
cmangle_OBJECTS = \
"CMakeFiles/cmangle.dir/cmangle.c.o"

# External object files for target cmangle
cmangle_EXTERNAL_OBJECTS =

cmangle: CMakeFiles/cmangle.dir/cmangle.c.o
cmangle: CMakeFiles/cmangle.dir/build.make
cmangle: libfmangle.a
cmangle: CMakeFiles/cmangle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/work/DanielStuff/nikolai/molecularGSM/build-gcc/CMakeFiles/CMakeTmp/fmangle/LOWER/UNDER/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable cmangle"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cmangle.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cmangle.dir/build: cmangle
.PHONY : CMakeFiles/cmangle.dir/build

CMakeFiles/cmangle.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cmangle.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cmangle.dir/clean

CMakeFiles/cmangle.dir/depend:
	cd /work/DanielStuff/nikolai/molecularGSM/build-gcc/CMakeFiles/CMakeTmp/fmangle/LOWER/UNDER && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/DanielStuff/nikolai/molecularGSM/cmake/tribits/core/config_tests/fmangle /work/DanielStuff/nikolai/molecularGSM/cmake/tribits/core/config_tests/fmangle /work/DanielStuff/nikolai/molecularGSM/build-gcc/CMakeFiles/CMakeTmp/fmangle/LOWER/UNDER /work/DanielStuff/nikolai/molecularGSM/build-gcc/CMakeFiles/CMakeTmp/fmangle/LOWER/UNDER /work/DanielStuff/nikolai/molecularGSM/build-gcc/CMakeFiles/CMakeTmp/fmangle/LOWER/UNDER/CMakeFiles/cmangle.dir/DependInfo.cmake
.PHONY : CMakeFiles/cmangle.dir/depend
