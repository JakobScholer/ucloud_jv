#
# Robust install cmake -P script that will ensure that 'cmake -P
# cmake_install.cmake' gets called in every package directory if the global
# 'cmake -P cmake_install.cmake' command fails.
#

set(PROJECT_BINARY_DIR /work/DanielStuff/nikolai/molecularGSM/BUILD)
set(TRIBITS_ENABLED_PACKAGES_BINARY_DIRS /work/DanielStuff/nikolai/molecularGSM/BUILD/GSM)

execute_process(
  COMMAND /work/DanielStuff/intel/oneapi/intelpython/python3.7/envs/nikolai/bin/cmake -P "cmake_install.cmake"
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
  RESULT_VARIABLE global_install_rtn)

if (NOT global_install_rtn EQUAL 0)

  message(SEND_ERROR "Error, full install failed!")

  message(
    "\n"
    "***\n"
    "*** The global install failed so resorting to package-by-package installs\n"
    "***\n"
    )

  foreach (pkg_binary_dir ${TRIBITS_ENABLED_PACKAGES_BINARY_DIRS})
    message("\nRunning install comamnds for '${pkg_binary_dir}/':")
    execute_process(COMMAND /work/DanielStuff/intel/oneapi/intelpython/python3.7/envs/nikolai/bin/cmake -P cmake_install.cmake
      WORKING_DIRECTORY ${pkg_binary_dir}
      RESULT_VARIABLE pkg_install_rtn)
  endforeach()

  message("\nRunning install comamnds for 'add_project_install_commands/':")
  execute_process(COMMAND /work/DanielStuff/intel/oneapi/intelpython/python3.7/envs/nikolai/bin/cmake -P cmake_install.cmake
    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/add_project_install_commands"
    RESULT_VARIABLE project_install_rtn)

  set(set_installed_group_and_permissions_file
    "${PROJECT_BINARY_DIR}/set_installed_group_and_permissions.cmake")
  if (EXISTS "${set_installed_group_and_permissions_file}")
    message("\nFixing up group and/or permissions of installed files and dirs:")
    execute_process(
      COMMAND /work/DanielStuff/intel/oneapi/intelpython/python3.7/envs/nikolai/bin/cmake -P "${set_installed_group_and_permissions_file}"
      RESULT_VARIABLE set_installed_rtn)
    if (NOT set_installed_rtn EQUAL 0)
      message(SEND_ERROR "ERROR: Fixup of group and permissions filed!")
    endif()
  endif()
  # NOTE: If the global install did not fail, then the group and permissions
  # would get fixed up since it gets run as part of the default 'install'
  # target.  Therefore, there is no need to run this if the 'install' target
  # passed.

endif()
