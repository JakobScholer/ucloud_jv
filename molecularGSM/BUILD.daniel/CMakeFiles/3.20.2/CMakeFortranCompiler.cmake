set(CMAKE_Fortran_COMPILER "/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/bin/intel64/ifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "2021.2.0.20210228")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/work/DanielStuff/intel/oneapi/vpl/2021.2.2/include;/work/DanielStuff/intel/oneapi/tbb/2021.2.0/include;/work/DanielStuff/intel/oneapi/mpi/2021.2.0/include;/work/DanielStuff/intel/oneapi/mkl/latest/include;/work/DanielStuff/intel/oneapi/ippcp/2021.2.0/include;/work/DanielStuff/intel/oneapi/ipp/2021.2.0/include;/work/DanielStuff/intel/oneapi/dpl/2021.2.0/linux/include;/work/DanielStuff/intel/oneapi/dpcpp-ct/2021.2.0/include;/work/DanielStuff/intel/oneapi/dnnl/2021.2.0/cpu_dpcpp_gpu_dpcpp/lib;/work/DanielStuff/intel/oneapi/dev-utilities/2021.2.0/include;/work/DanielStuff/intel/oneapi/dal/2021.2.0/include;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/include;/work/DanielStuff/intel/oneapi/ccl/2021.2.0/include/cpu_gpu_dpcpp;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/compiler/include/intel64;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/compiler/include/icc;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/compiler/include;/usr/local/include;/usr/lib/gcc/x86_64-linux-gnu/7/include;/usr/lib/gcc/x86_64-linux-gnu/7/include-fixed;/usr/include;/usr/include/x86_64-linux-gnu")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/work/DanielStuff/intel/oneapi/vpl/2021.2.2/lib;/work/DanielStuff/intel/oneapi/tbb/2021.2.0/lib/intel64/gcc4.8;/work/DanielStuff/intel/oneapi/mpi/2021.2.0/libfabric/lib;/work/DanielStuff/intel/oneapi/mpi/2021.2.0/lib/release;/work/DanielStuff/intel/oneapi/mpi/2021.2.0/lib;/work/DanielStuff/intel/oneapi/mkl/latest/lib/intel64;/work/DanielStuff/intel/oneapi/ippcp/2021.2.0/lib/intel64;/work/DanielStuff/intel/oneapi/ipp/2021.2.0/lib/intel64;/work/DanielStuff/intel/oneapi/dnnl/2021.2.0/cpu_dpcpp_gpu_dpcpp/lib;/work/DanielStuff/intel/oneapi/dal/2021.2.0/lib/intel64;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/lib;/work/DanielStuff/intel/oneapi/clck/2021.2.0/lib/intel64;/work/DanielStuff/intel/oneapi/ccl/2021.2.0/lib/cpu_gpu_dpcpp;/usr/lib/gcc/x86_64-linux-gnu/7;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib64;/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
