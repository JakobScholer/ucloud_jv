set(CMAKE_C_COMPILER "/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/bin/intel64/icc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "Intel")
set(CMAKE_C_COMPILER_VERSION "2021.2.0.20210228")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "GNU")
set(CMAKE_C_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_C_SIMULATE_VERSION "7.5.0")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCC )
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_C_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/work/DanielStuff/intel/oneapi/vpl/2021.2.2/include;/work/DanielStuff/intel/oneapi/tbb/2021.2.0/include;/work/DanielStuff/intel/oneapi/mpi/2021.2.0/include;/work/DanielStuff/intel/oneapi/mkl/latest/include;/work/DanielStuff/intel/oneapi/ippcp/2021.2.0/include;/work/DanielStuff/intel/oneapi/ipp/2021.2.0/include;/work/DanielStuff/intel/oneapi/dpl/2021.2.0/linux/include;/work/DanielStuff/intel/oneapi/dpcpp-ct/2021.2.0/include;/work/DanielStuff/intel/oneapi/dnnl/2021.2.0/cpu_dpcpp_gpu_dpcpp/lib;/work/DanielStuff/intel/oneapi/dev-utilities/2021.2.0/include;/work/DanielStuff/intel/oneapi/dal/2021.2.0/include;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/include;/work/DanielStuff/intel/oneapi/ccl/2021.2.0/include/cpu_gpu_dpcpp;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/compiler/include/intel64;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/compiler/include/icc;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/compiler/include;/usr/local/include;/usr/lib/gcc/x86_64-linux-gnu/7/include;/usr/lib/gcc/x86_64-linux-gnu/7/include-fixed;/usr/include;/usr/include/x86_64-linux-gnu")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "imf;svml;irng;m;ipgo;decimal;cilkrts;stdc++;gcc;gcc_s;irc;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/work/DanielStuff/intel/oneapi/vpl/2021.2.2/lib;/work/DanielStuff/intel/oneapi/tbb/2021.2.0/lib/intel64/gcc4.8;/work/DanielStuff/intel/oneapi/mpi/2021.2.0/libfabric/lib;/work/DanielStuff/intel/oneapi/mpi/2021.2.0/lib/release;/work/DanielStuff/intel/oneapi/mpi/2021.2.0/lib;/work/DanielStuff/intel/oneapi/mkl/latest/lib/intel64;/work/DanielStuff/intel/oneapi/ippcp/2021.2.0/lib/intel64;/work/DanielStuff/intel/oneapi/ipp/2021.2.0/lib/intel64;/work/DanielStuff/intel/oneapi/dnnl/2021.2.0/cpu_dpcpp_gpu_dpcpp/lib;/work/DanielStuff/intel/oneapi/dal/2021.2.0/lib/intel64;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin;/work/DanielStuff/intel/oneapi/compiler/2021.2.0/linux/lib;/work/DanielStuff/intel/oneapi/clck/2021.2.0/lib/intel64;/work/DanielStuff/intel/oneapi/ccl/2021.2.0/lib/cpu_gpu_dpcpp;/usr/lib/gcc/x86_64-linux-gnu/7;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib64;/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")