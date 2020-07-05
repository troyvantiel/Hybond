cmake_minimum_required(VERSION 2.8)

project(CHARMM C CXX Fortran)
set(CHARMM_VERSION_MAJOR 44)
set(CHARMM_VERSION_MINOR 2)

# for the intel fortran compiler
# the fortran compiler version number is necessary to
# set the correct openmp flag which dropped
# -openmp in favor of -qopenmp
# for version 18
# and the fortran compiler version number is not set
# for CMake versions below 3.3.1
if( (NOT DEFINED CMAKE_Fortran_COMPILER_VERSION) AND
    (CMAKE_Fortran_COMPILER_ID MATCHES "Intel") )
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --version
    OUTPUT_VARIABLE output
    ERROR_VARIABLE output)
  string(REGEX MATCH "ifort \\(IFORT\\) ([0-9\\.]+) [0-9]+"
    outver
    "${output}")
  set(CMAKE_Fortran_COMPILER_VERSION "${CMAKE_MATCH_1}")
  message(STATUS
    "Fortran compiler version set to ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CHARMM_SOURCE_DIR}/tool/cmake")

if(NOT DEFINED add_keywords)
  set(add_keywords "")
endif()
set(add_keywords "${add_keywords}"
  CACHE STRING "prefx keywords to add to pref.dat, comma delimited" FORCE)
string(REPLACE "," ";" add_keywords "${add_keywords}")

if(NOT DEFINED remove_keywords)
  set(remove_keywords "")
endif()
set(remove_keywords "${remove_keywords}"
  CACHE STRING "prefx keys to remove to pref.dat, comma delimited" FORCE)
string(REPLACE "," ";" remove_keywords "${remove_keywords}")

# default ON features to compile in
option(colfft "COLFFT support" ON)
option(cuda "Nvidia CUDA support" ON)
option(exafmm "ExaFMM support" ON)
option(fftw "FFTW support" ON)
option(mkl "Intel MKL support" ON)
option(mpi "MPI support" ON)
option(opencl "OpenCL support" ON)
option(openmm "OpenMM support" ON)
option(openmp "OpenMP support" ON)
option(qchem "Q-Chem support" ON)
option(quantum "QUANTUM support" ON)

# default OFF features which users must specifically request
option(abpo "ABPO support; requires MPI and activates ENSEMBLE")
option(gamus "GAMUS support")
option(ensemble "support for ENSEMBLE replicas; requires MPI")
option(in_place_install "install binaries and libraries in the source root directory")
option(lite "compile with minimal features")
option(nih "use the NIH set of extra pref.dat keywords")
option(pipf "support for Polarizable Intermolecular Potential Function")
option(repdstr "repdstr support")
option(static "link only static libs for a portable executable")
option(stringm "stringm support")
option(tsri "use the TSRI set of extra pref.dat keywords")
option(x11 "X11 graphics support")

# not available in free version
set(domdec OFF)
set(domdec_gpu OFF)

# QM/MM options begin

# mutually exclusive default OFF qm/mm options
option(squantm "replace quantum with squantum")
option(sccdftb "replace quantum with sccdftb")
option(g09 "replace quantum with Gaussian09")
option(qturbo "replace quantum with qturbo")
option(mndo97 "replace quantum with mndo97")
option(qmmmsemi "replace quantum with qmmmsemi (amber semi-empirical qmmm")

foreach(opt lite squantm sccdftb g09 qturbo mndo97 qmmmsemi gamess)
  if(quantum AND ${opt})
    message(WARNING "${opt} and quantum are incompatible; switching quantum OFF")
    set(quantum OFF)
  endif()
endforeach(opt)

foreach(opt lite squantm mndo97)
  if(qchem AND ${opt})
    message(WARNING "${opt} and qchem are incompatible; switching qchem OFF")
    set(qchem OFF)
  endif()
endforeach(opt)

function(test_option_conflict opt1 opt2)
  if(${opt1})
    message(FATAL_ERROR "${opt1} and ${opt2} are incompatible; choose one and try again")
  endif()
endfunction(test_option_conflict)

# begin sanity checks on incompatible QM builds
if(quantum)
  test_option_conflict(gamess quantum) 
  test_option_conflict(sccdftb quantum)
  test_option_conflict(mndo97 quantum)
  test_option_conflict(squantm quantum)
endif(quantum)

if(mndo97)
  test_option_conflict(gamess mndo97)
  test_option_conflict(sccdftb mndo97)
  test_option_conflict(squantm mndo97)
endif(mndo97)
# end sanity checks on incompatible QM builds
# QM/MM options end

# handle scc and colfft dependencies and conflicts
if(sccdftb AND colfft)
  message(WARNING "sccdftb and COLFFT are incompatible; switching COLFFT OFF")
  set(colfft OFF)
endif()

if(NOT colfft)
  message(WARNING "colfft required for mkl and fftw; switching mkl and fftw OFF")
  set(mkl OFF)
  set(fftw OFF)
endif()

if(NOT colfft)
  message(WARNING "colfft required for mkl and fftw; switching mkl and fftw OFF")
  set(mkl OFF)
  set(fftw OFF)
endif()

# begin set the compiler opt/arch flags
# begin optimization level
set(default_opt 2)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(default_opt 3)
endif()

# for release builds only:
# level to use with -O for most fortran files
# (0, 1, 2, 3)
if(NOT DEFINED optimization)
  set(optimization ${default_opt})
endif()
# end optimization level

# begin architecture
# begin C architecture flags
set(default_arch_flag -march=native)  # GNU & clang flag is default
if(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(default_arch_flag -xHost)
elseif(CMAKE_C_COMPILER_ID MATCHES "PGI")
  set(default_arch_flag -fastsse)
endif()

set(arch_msg "architecture flag for C compiler")
if(arch_flag_C)
  set(arch_flag_C ${arch_flag_C})
elseif(arch_flag_all)
  set(arch_flag_C ${arch_flag_all})
else()
  set(arch_flag_C ${default_arch_flag})
endif()
# end C architecture flags

# begin CXX architecture flags
set(default_arch_flag -march=native)  # GNU & clang flag is default
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(default_arch_flag -xHost)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  set(default_arch_flag -fastsse)
endif()

set(arch_msg "architecture flag for C++ compiler")
if(arch_flag_CXX)
  set(arch_flag_CXX ${arch_flag_CXX})
elseif(arch_flag_all)
  set(arch_flag_CXX ${arch_flag_all})
else()
  set(arch_flag_CXX ${default_arch_flag})
endif()
# end CXX architecture flags

# begin fortran architecture flags
set(default_arch_flag)  # GNU & clang flag is default
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(default_arch_flag -axSSE4.1)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(default_arch_flag)
endif()

set(arch_msg "architecture flag for fortran compiler")
if(arch_flag_Fortran)
  set(arch_flag_Fortran ${arch_flag_Fortran})
elseif(arch_flag_all)
  set(arch_flag_Fortran ${arch_flag_all})
else()
  set(arch_flag_Fortran ${default_arch_flag})
endif()
# end fortran architecture flags
# end architecture
# end set the compiler opt/arch flags

if(NOT static)  # Setup installation RPATH
  if(APPLE)
    set(CMAKE_MACOSX_RPATH 1)
  endif()
  
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
    "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  ENDIF("${isSystemDir}" STREQUAL "-1")
else(NOT static)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
endif(NOT static)  # END Setup installation RPATH

# ensure that lite build comes with no added features/keywords
if(NOT lite)

  if(x11)
    find_package(X11)
  endif(x11)

  if(NOT repdstr)
    if(mkl)
      find_package(MKL)
    endif(mkl)

    if(fftw AND (NOT MKL_FOUND))
      find_package(FFTW)
      find_package(FFTWF)
      if(NOT FFTWF_FOUND)
        message(WARNING "Unable to find single precision FFTW")
      endif()
    endif()
  endif(NOT repdstr)

  if(colfft AND repdstr)
    message(WARNING "repdstr is incompatible with COLFFT; switching COLFFT OFF")
    set(colfft OFF)
  endif()
  
  if(colfft AND (NOT MKL_FOUND) AND (NOT FFTW_FOUND))
    message(WARNING "FFTW or MKL required for COLFFT; switching COLFFT OFF")
    set(colfft OFF)
  endif()

  if(mpi)
    find_package(MPI)
  endif()

  if(openmp)
    find_package(OpenMP)
    find_package(OpenMP_Fortran)

    # fix openmp flags for old versions of cmake
    # and new versions of intel icc
    if(UNIX AND OPENMP_FOUND)
      if(CMAKE_C_COMPILER_ID MATCHES "Intel" AND
          (NOT (CMAKE_C_COMPILER_VERSION VERSION_LESS "15")) AND
          OpenMP_C_FLAGS STREQUAL "-openmp")
        set(OpenMP_C_FLAGS "-qopenmp")
        message(STATUS "OpenMP C flag changed to -qopenmp")
      endif()
      if(CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND
          (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15")) AND
          OpenMP_CXX_FLAGS STREQUAL "-openmp")
        set(OpenMP_CXX_FLAGS "-qopenmp")
        message(STATUS "OpenMP CXX flag changed to -qopenmp")
      endif()
    endif()

    if(APPLE AND NOT OPENMP_FOUND)
      find_package(AppleOpenMP)
      if(AppleOpenMP_FOUND)
        message(STATUS "Apple OpenMP library found at ${APPLE_OPENMP_LIBRARY}")
        message(STATUS "Apple OpenMP C and C++ flags set to ${APPLE_OPENMP_CFLAGS}")
      endif()
    endif()    
  endif(openmp)

  if(openmm)
    find_package(OpenMM COMPONENTS CPU CUDA OPENCL)
  endif()

  if(exafmm)
    find_package(ExaFMM)
  endif()
else(NOT lite)
  set(colfft OFF)
  set(cuda OFF)
  set(domdec OFF)
  set(exafmm OFF)
  set(fftw OFF)
  set(ljpme OFF)
  set(mkl OFF)
  set(mpi OFF)
  set(opencl OFF)
  set(openmm OFF)
  set(openmp OFF)
  set(qchem OFF)
  set(quantum OFF)
  message(WARNING "lite build chosen; all features are set to OFF")
endif(NOT lite)

if(EXAFMM_FOUND)
  link_directories(${ExaFMM_LIBRARIES})
endif()

if(repdstr AND domdec)
  message(WARNING "repdstr is incompatible with DOMDEC; switching DOMDEC OFF")
  set(domdec OFF)
endif()

if((lite) AND domdec)
  message(WARNING "lite is incompatible with DOMDEC; switching DOMDEC OFF")
  set(domdec OFF)
endif()

if((NOT MPI_Fortran_FOUND) AND domdec)
  message(WARNING "MPI is required for DOMDEC; switching DOMDEC OFF")
  set(domdec OFF)
endif()

if((NOT colfft) AND domdec)
  message(WARNING "COLFFT is required for DOMDEC; switching DOMDEC OFF")
  set(domdec OFF)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are Debug Release."
    FORCE)

if(cuda AND (NOT lite))
  find_package(CUDA)
endif()

set(needs_abi_zero FALSE)
if(OPENMM_FOUND)
  try_compile(can_link_omm
    "${CMAKE_BINARY_DIR}/omm_abi_test"
    "${CMAKE_SOURCE_DIR}/tool/cmake/omm_abi_test.cpp"
    LINK_LIBRARIES OpenMM
    CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES=${OPENMM_INCLUDE_DIR}"
    "-DLINK_DIRECTORIES=${OPENMM_LIBRARY_DIR}")

  if(NOT can_link_omm)
    try_compile(needs_abi_zero
      "${CMAKE_BINARY_DIR}/omm_abi_test"
      "${CMAKE_SOURCE_DIR}/tool/cmake/omm_abi_test.cpp"
      COMPILE_DEFINITIONS "-D_GLIBCXX_USE_CXX11_ABI=0"
      LINK_LIBRARIES OpenMM
      CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES=${OPENMM_INCLUDE_DIR}"
      "-DLINK_DIRECTORIES=${OPENMM_LIBRARY_DIR}")
  endif()

  if((NOT can_link_omm) AND (NOT needs_abi_zero))
    message(WARNING "test to link to OpenMM library failed"
      "\nplease check OpenMM installation and environment variables"
      "\nand try again if compilation fails")
  endif()

  set(omm_std_opt)
  foreach(opt "" "-std=c++11" "-std=c++0x")
      try_compile(can_omm
        "${CMAKE_BINARY_DIR}/omm_cxx11_test"
        "${CMAKE_SOURCE_DIR}/tool/cmake/omm_cxx11_test.cpp"
        COMPILE_DEFINITIONS "${opt}"
        CMAKE_FLAGS
        "-DINCLUDE_DIRECTORIES=${OPENMM_INCLUDE_DIR}")
      
      if(can_omm)
        set(omm_std_opt "${opt}")
        break()
      endif(can_omm)
  endforeach(opt)

  if(NOT can_omm)
    message(WARNING "test to compile to OpenMM plugins failed"
      "\nplease check OpenMM installation and environment variables"
      "\nand try again if compilation fails")
  endif()
endif(OPENMM_FOUND)

if(domdec_gpu)
  if(lite)
    message(FATAL_ERROR
      "The lite build is incompatible with DOMDEC_GPU.")
  endif()

  if(NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR
      "MPI was not found and is required for DOMDEC_GPU.")
  endif()

  if(NOT domdec)
    message(FATAL_ERROR
      "DOMDEC is OFF but is required for DOMDEC_GPU.")
  endif()

  if(NOT CUDA_FOUND)
    message(FATAL_ERROR
      "CUDA support was not found and is required for DOMDEC_GPU.")
  endif()

  set(CUDA_PROPAGATE_HOST_FLAGS OFF)

  execute_process(COMMAND ${CUDA_HOST_COMPILER} --version
    OUTPUT_VARIABLE nvcc_cxx_version)
  string(FIND "${nvcc_cxx_version}" "GCC" gnu_find)
  string(FIND "${nvcc_cxx_version}" "Intel" intel_find)
  string(FIND "${nvcc_cxx_version}" "clang" clang_find)
  string(FIND "${nvcc_cxx_version}" "PGI" pgi_find)

  set(NVCC_HOST_COMPILER_ID "UNKNOWN")
  if(${gnu_find} GREATER -1)
    set(NVCC_HOST_COMPILER_ID "GNU")
  elseif(${intel_find} GREATER -1)
    set(NVCC_HOST_COMPILER_ID "Intel")
  elseif(${clang_find} GREATER -1)
    set(NVCC_HOST_COMPILER_ID "Clang")
  elseif(${pgi_find} GREATER -1)
    set(NVCC_HOST_COMPILER_ID "PGI")
  endif()

  # begin host architecture flags
  set(default_arch_flag -march=native)  # GNU & clang flag is default
  if(NVCC_HOST_COMPILER_ID MATCHES "Intel")
    set(default_arch_flag -xHost)
  elseif(NVCC_HOST_COMPILER_ID MATCHES "PGI")
    set(default_arch_flag -fastsse)
  endif()

  set(arch_msg "architecture flag for nvcc/cuda host compiler")
  if(arch_flag_host)
    set(arch_flag_host ${arch_flag_host})
  elseif(arch_flag_all)
    set(arch_flag_CXX ${arch_flag_all})
  else()
    set(arch_flag_host ${default_arch_flag})
  endif()
  # end host architecture flags
  
  set(nvcc_host_flags)
  if(NVCC_HOST_COMPILER_ID MATCHES "GNU")
    set(nvcc_host_flags "-Xcompiler;${arch_flag_host},-fopenmp")
  elseif(NVCC_HOST_COMPILER_ID MATCHES "Clang")
    set(nvcc_host_flags "-Xcompiler;${arch_flag_host}")
    if(AppleOpenMP_FOUND)
      set(nvcc_host_flags "${nvcc_host_flags},-Xpreprocessor,-fopenmp")
    endif()
  elseif(NVCC_HOST_COMPILER_ID MATCHES "Intel")
    string(REGEX MATCH "icp?c \\(ICC\\) ([0-9\\.]+) [0-9]+"
      outver
      "${nvcc_cxx_version}")
    set(nvcc_host_ver "${CMAKE_MATCH_1}")

    if(nvcc_host_ver VERSION_LESS "15")
      set(OpenMP_NVCC_FLAGS "-openmp")
    else()
      set(OpenMP_NVCC_FLAGS "-qopenmp")
    endif()

    set(nvcc_host_flags "-Xcompiler;-mp1,${arch_flag_host}")
    if(OpenMP_NVCC_FLAGS)
      set(nvcc_host_flags "${nvcc_host_flags},${OpenMP_NVCC_FLAGS}")
    endif()
  elseif(NVCC_HOST_COMPILER_ID MATCHES "PGI")
    set(nvcc_host_flags "-Xcompiler;${arch_flag_host},-mp")
  endif()

  if(needs_abi_zero AND nvcc_host_flags)
    set(nvcc_host_flags "${nvcc_host_flags},-D_GLIBCXX_USE_CXX11_ABI=0")
  endif()

  set(local_nvcc_flags)
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    list(APPEND local_nvcc_flags -g -O0 -G)
  elseif(CMAKE_BUILD_TYPE MATCHES Release)
    list(APPEND local_nvcc_flags -O3)
  endif()

  if(NOT (nvcc_ptx_target OR nvcc_arch_target))
    list(APPEND gencode_flags
      -gencode arch=compute_30,code=compute_30)
  endif()

  if(nvcc_ptx_target)
    list(APPEND gencode_flags
      -gencode arch=compute_${nvcc_ptx_target},code=compute_${nvcc_ptx_target})
  endif(nvcc_ptx_target)

  if(nvcc_arch_target)
    list(APPEND gencode_flags
      -gencode arch=compute_${nvcc_arch_target},code=sm_${nvcc_arch_target})
  endif(nvcc_arch_target)
  
  list(APPEND local_nvcc_flags
    -fmad=true
    -use_fast_math
    ${nvcc_host_flags}
    ${gencode_flags})
  if(NOT CHARMM_NVCC_FLAGS)
    set(CHARMM_NVCC_FLAGS "${local_nvcc_flags}"
      CACHE STRING "flags for the CUDA compiler" FORCE)
  endif()
  set(CUDA_NVCC_FLAGS "${CHARMM_NVCC_FLAGS}")

  find_library(NSIGHT_LIBRARIES
    nvToolsExt
    HINTS ${CUDA_TOOLKIT_ROOT_DIR}/lib64 ${CUDA_TOOLKIT_ROOT_DIR}/lib)
  if(NOT NSIGHT_LIBRARIES)
    message(WARNING
      "Nvidia NSIGHT library could not be found.")
    message(FATAL_ERROR
      "Nvidia NSIGHT library is required for DOMDEC_GPU.")
  endif()
endif(domdec_gpu)

if(opencl AND (NOT lite))
  if(CUDA_FOUND)
    set(ENV{CUDA_PATH} ${CUDA_TOOLKIT_ROOT_DIR})
  endif()
  find_package(OpenCL QUIET)
  if(OpenCL_FOUND)
    message(STATUS "OpenCL include dir found: ${OpenCL_INCLUDE_DIR}")
    message(STATUS "OpenCL library found: ${OpenCL_LIBRARY}")
  else()
    message(STATUS "OpenCL not found")
  endif()
endif()

if(OPENMM_FOUND AND (NOT lite))
  list(APPEND CMAKE_INSTALL_RPATH ${OPENMM_LIBRARIES})

  get_filename_component(OPENMM_PLUGIN_DIR ${OPENMM_LIBRARY} PATH)
  set(OPENMM_PLUGIN_DIR "${OPENMM_PLUGIN_DIR}/plugins")

  # CHARMM_PLUGIN_INSTALL_DIR: suffix for use in sub CMakeLists.txt
  # CHARMM_PLUGIN_DIR: absolute path for use in fortran code
  set(CHARMM_PLUGIN_INSTALL_DIR "lib")
  set(CHARMM_PLUGIN_DIR
    "${CMAKE_INSTALL_PREFIX}/${CHARMM_PLUGIN_INSTALL_DIR}")
  configure_file(source/openmm/plugin_locs.f90.in plugin_locs.f90)

  # is openmm version >= 7.1? api change happend again
  find_file(OMM_BAROSTAT_HEADER MonteCarloBarostat.h
    PATHS ${OPENMM_INCLUDE_DIRS}
    PATH_SUFFIXES openmm)
  file(STRINGS ${OMM_BAROSTAT_HEADER} OPENMM_API_UPDATE
    REGEX getDefaultTemperature)
  if(OPENMM_API_UPDATE)
    set(OPENMM_API_UPDATE TRUE)
  else()
    set(OPENMM_API_UPDATE FALSE)
  endif()
endif()

if(repdstr)
  if(lite)
    message(FATAL_ERROR
      "The lite build is incompatible with REPDSTR.")
  endif()

  if(domdec OR domdec_gpu)
    message(FATAL_ERROR
      "DOMDEC is incompatible with REPDSTR.")
  endif()

  if(NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR
      "MPI was not found and is required for REPDSTR.")
  endif()
endif(repdstr)

if(stringm)
  if(lite)
    message(FATAL_ERROR
      "The lite build is incompatible with STRINGM.")
  endif()

  if(NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR
      "MPI was not found and is required for STRINGM.")
  endif()

  if(domdec)
    message(WARNING
      "STRINGM and DOMDEC can be compiled together but cannot be used simultaneously")
  endif()
endif(stringm)

find_program(GIT_COMMAND git)
find_program(SVN_COMMAND svnversion)

if(IS_DIRECTORY "${CHARMM_SOURCE_DIR}/.git" AND EXISTS ${GIT_COMMAND})
  execute_process(COMMAND ${GIT_COMMAND} rev-parse --short HEAD
    WORKING_DIRECTORY ${CHARMM_SOURCE_DIR}
    OUTPUT_VARIABLE REVISION_INFO
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(REVISION_INFO "Git commit ID ${REVISION_INFO}")
elseif(IS_DIRECTORY "${CHARMM_SOURCE_DIR}/.svn" AND EXISTS ${SVN_COMMAND})
  execute_process(COMMAND ${SVN_COMMAND}
    WORKING_DIRECTORY ${CHARMM_SOURCE_DIR}
    OUTPUT_VARIABLE REVISION_INFO
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(REVISION_INFO "SVN revision ${REVISION_INFO}")
else()
  set(REVISION_INFO "686801d1a")
endif()

configure_file(source/charmm/revid.h.in revid.h)

include(tool/cmake/prefx_keywords.cmake)

list(SORT keywords)
list(REMOVE_DUPLICATES keywords)
list(LENGTH keywords NUM_PREFX_KEYWORDS)

list(GET keywords -1 last_word)
list(REMOVE_AT keywords -1)

set(PREFX_KEYWORDS_STRINGS)
foreach(word ${keywords})
  set(PREFX_KEYWORDS_STRINGS
    "${PREFX_KEYWORDS_STRINGS}\"${word}\", &\n       ")
endforeach()
set(PREFX_KEYWORDS_STRINGS
  "${PREFX_KEYWORDS_STRINGS}\"${last_word}\"")

list(APPEND keywords ${last_word})

configure_file(source/util/keywords.inc.in keywords.inc)

if(lite)
  message(STATUS "lite build chosen")
else()
  message(STATUS "full build chosen")
endif()

configure_file(tool/cmake/pref.dat.in pref.dat)

# fortran_files.cmake holds the list charmm_src_files
# just a list of all the source/*/*.F90 files
# the configure shell script may append files to this list
include(tool/cmake/fortran_files.cmake)

if(sccdftb)
  # these files will be included later from
  # tool/cmake/sccdftb_fortran_files.cmake
  # these files do not need flags for OpenMP
  # and may not perform correctly with -O3
  list(REMOVE_ITEM charmm_src_files
    source/sccdftbint/sccdftbini.F90
    source/sccdftbint/sccpme.F90
    source/sccdftbint/sccpmeutil.F90
    source/sccdftbint/stbgho.F90
    source/sccdftbint/stbgho_ltm.F90)
endif()

if(CMAKE_GENERATOR MATCHES Ninja)
  set(charmm_fortran_files ${charmm_src_files})
else()
  set(charmm_fortran_files)

  function(process_src_file from_file to_file)
    add_custom_command(
      OUTPUT "${to_file}"
      COMMAND ${CMAKE_COMMAND} -E copy "${from_file}" "${to_file}"
      MAIN_DEPENDENCY "${from_file}")
    set_source_files_properties(${to_file} PROPERTIES GENERATED TRUE)
  endfunction(process_src_file)

  foreach(from_file ${charmm_src_files})
    get_filename_component(base ${from_file} NAME_WE)
    set(to_file "${CMAKE_CURRENT_BINARY_DIR}/${base}.F90")

    process_src_file("${CMAKE_SOURCE_DIR}/${from_file}" "${to_file}")

    list(APPEND charmm_fortran_files "${to_file}")
  endforeach()
endif()

# list of C language files named charmm_c_files
include(tool/cmake/c_files.cmake)

# list of C++ language files named charmm_cxx_files
include(tool/cmake/cxx_files.cmake)

list(FIND keywords "MTS" _index)
if(domdec AND (${_index} GREATER -1))
  message(FATAL_ERROR
    "DOMDEC and MTS are incompatible. "
    "If you are using the configure shell script "
    "try --without-domdec or remove -a MTS. "
    "If you are using cmake "
    "try -Ddomdec=OFF.")
endif()

if(domdec_gpu)
  # list of cuda files named domdec_cuda_files
  # source/domdec*/*.cu plus associated header files
  include(tool/cmake/domdec_cuda_files.cmake)
  cuda_add_library(charmm_cuda
    ${domdec_cuda_files} ${domdec_gpu_cuda_files})
endif(domdec_gpu)

if(qchem)
  add_executable(qchem tool/qchem.f90)
endif(qchem)

if(sccdftb)
  include(tool/cmake/sccdftb_fortran_files.cmake)
endif(sccdftb)

if(CMAKE_GENERATOR MATCHES Ninja)
  set(CHARMM_MAIN_FILE source/charmm/charmm_main.F90)
  list(REMOVE_ITEM charmm_fortran_files
    "source/charmm/charmm_main.F90")
else()
  set(CHARMM_MAIN_FILE charmm_main.F90)
  list(REMOVE_ITEM charmm_fortran_files
    "${CMAKE_CURRENT_BINARY_DIR}/charmm_main.F90")
endif()

if(NOT sccdftb)
  add_library(charmm_fortran ${charmm_fortran_files})
else()
  add_library(charmm_fortran ${charmm_fortran_files} ${sccdftb_fortran_files})
endif()

if(sccdftb)
  set_property(SOURCE ${sccdftb_fortran_files}
    APPEND_STRING PROPERTY COMPILE_FLAGS "-O1 ")
endif(sccdftb)

if(charmm_c_files)
  add_library(charmm_c ${charmm_c_files})
endif()

if(charmm_cxx_files)
  add_library(charmm_cxx ${charmm_cxx_files})
endif()

add_executable(charmm ${CHARMM_MAIN_FILE})

if(static)
  set_target_properties(charmm PROPERTIES LINK_SEARCH_START_STATIC 1)
  set_target_properties(charmm PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_source_files_properties(
    "source/machdep/cstuff.c"
    "source/machdep/uninf.cpp"
    PROPERTIES COMPILE_DEFINITIONS "STATIC=1")
endif(static)

target_link_libraries(charmm charmm_fortran)

if(charmm_c_files)
  target_link_libraries(charmm charmm_c)
endif()

if(domdec_gpu)
  target_link_libraries(charmm charmm_cuda)
endif()

if(domdec_gpu AND (CMAKE_Fortran_COMPILER_ID MATCHES "PGI"))
  if(NVCC_HOST_COMPILER_ID MATCHES "GNU")
    target_link_libraries(charmm gomp)
  elseif(NVCC_HOST_COMPILER_ID MATCHES "Intel")
    target_link_libraries(charmm iomp5)
  endif()
endif()

if(charmm_cxx_files)
  target_link_libraries(charmm charmm_cxx)
endif()

set_property(TARGET charmm PROPERTY LINKER_LANGUAGE Fortran)

if(gamus)
  if(lite)
    message(FATAL_ERROR
      "The lite build is incompatible with GAMUS.")
  endif()

  find_package(LAPACK)

  if(NOT LAPACK_FOUND)
    message(FATAL_ERROR
      "LAPACK support was not found and is required for GAMUS.")
  endif()

  target_link_libraries(charmm ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
endif(gamus)

if(CMAKE_GENERATOR MATCHES Ninja)
  list(APPEND LOWER_OPT_FILES
    source/dynamc/dynamc.F90
    source/energy/ecmap.F90
    source/energy/gbmvmodule.F90
    source/misc/eef1.F90)
else()
  list(APPEND LOWER_OPT_FILES
    dynamc.F90
    ecmap.F90
    eef1.F90
    gbmvmodule.F90)
endif()

# these files cause Intel 16.x to choke on -O3
if((CMAKE_BUILD_TYPE MATCHES Release)
  AND (CMAKE_Fortran_COMPILER_ID MATCHES "Intel"))
  set_source_files_properties(
    ${LOWER_OPT_FILES}
    PROPERTIES COMPILE_FLAGS -O2)
endif()

# set prefx keywords as compile definitions for fortran
# see tool/cmake/prefx_kewords.cmake for definition of 'keywords' list
foreach(keyword ${keywords})
  target_compile_definitions(charmm_fortran PRIVATE -DKEY_${keyword}=1 )
  target_compile_definitions(charmm PRIVATE -DKEY_${keyword}=1 )
endforeach()

# add GAMESS compilation for QM/MM CHARMM, if configured (configure -q)
# for parallel compile use:
# "time (make -j 8;make -j 8)", otherwise just: make
# charmm_frotran needs to be specified again after gamess
# even if it is already above
if(gamess)
  add_subdirectory( source/gamint )
  target_link_libraries(charmm gamess gmsc charmm_fortran)
endif(gamess)

# add NWChem compilation for QM/MM CHARMM, if configured (configure --nwchem)
if(nwchem)
  if (EXISTS $ENV{NWCHEM_LIB_CHARMM})
    message("-- Using NWChem libraries.")
  else (EXISTS $ENV{NWCHEM_LIB_CHARMM})
    message(FATAL_ERROR
      "NWCHEM_LIB_CHARMM not defined or nonexistent NWChem libraries.")
  endif (EXISTS $ENV{NWCHEM_LIB_CHARMM})
  target_link_libraries(charmm
    $ENV{NWCHEM_LIB_CHARMM}/../../src/nwchem.o
    $ENV{NWCHEM_LIB_CHARMM}/../../src/stubs.o
    -L$ENV{NWCHEM_LIB_CHARMM}
    -L$ENV{NWCHEM_LIB_CHARMM}/../../src/tools/install/lib
    -lnwctask -lccsd -lmcscf -lselci -lmp2 -lmoints -lstepper
    -ldriver -loptim -lnwdft -lgradients -lcphf -lesp -lddscf
    -lguess -lhessian -lvib -lnwcutil -lrimp2 -lproperty -lsolvation
    -lnwints -ldplot -lvscf -letrans -ltce -lbq -lcons -lnwcutil
    -lga -larmci -lpeigs -lperfm -lcons -lbq -lnwcutil
    -lnwclapack  -lnwcblas
    )
endif(nwchem)

######
# BEGIN set compile and link flags

# make sure that the default is a RELEASE

if(APPLE)
  SET(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols <TARGET>")
endif(APPLE)

set(link_flags)
if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(link_flags "-g")
endif()

set(fortran_compile_flags)
set(debug_fortran_compile_flags "-g -O0")
if(domdec)
  set(debug_fortran_compile_flags
    "${debug_fortran_compile_flags} ${arch_flag_Fortran}")
endif()
set(release_fortran_compile_flags "-O${optimization} ${arch_flag_Fortran}")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "4.5.0")
    set(link_flags "${link_flags} -fno-use-linker-plugin")
  endif()
  set(fortran_compile_flags
    "${fortran_compile_flags} -DGNU -ffixed-line-length-none -ffree-line-length-none -fmax-identifier-length=63")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(fortran_compile_flags
    "${fortran_compile_flags} -132 -align all -heap-arrays 256 -assume byterecl")
  set(debug_fortran_compile_flags
    "${debug_fortran_compile_flags} -u -traceback")
  set(release_fortran_compile_flags
    "${release_fortran_compile_flags} -mp1 -fp-model strict")
  if(APPLE)
    SET(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  endif(APPLE)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(link_flags "${link_flags} -lgcc_s")
  set(fortran_compile_flags
    "${fortran_compile_flags} -Mbackslash")
  set(release_fortran_compile_flags
    "${release_fortran_compile_flags} -Munroll -Mnoframe ")
endif()

if(NOT CHARMM_Fortran_FLAGS)
  set(CHARMM_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${fortran_compile_flags}"
    CACHE STRING "Fortran flags" FORCE)
endif()
set(CMAKE_Fortran_FLAGS "${CHARMM_Fortran_FLAGS}")

if(NOT CHARMM_Fortran_FLAGS_DEBUG)
  set(CHARMM_Fortran_FLAGS_DEBUG "${debug_fortran_compile_flags}"
    CACHE STRING "to use when CMAKE_BUILD_TYPE=Debug" FORCE)
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "${CHARMM_Fortran_FLAGS_DEBUG}")

if(NOT CHARMM_Fortran_FLAGS_RELEASE)
  set(CHARMM_Fortran_FLAGS_RELEASE "${release_fortran_compile_flags}"
    CACHE STRING "to use when CMAKE_BUILD_TYPE=Release" FORCE)
endif()
set(CMAKE_Fortran_FLAGS_RELEASE "${CHARMM_Fortran_FLAGS_RELEASE}")

set(c_compile_flags)
set(debug_c_compile_flags "-g -O0")
if(domdec)
  set(debug_c_compile_flags
    "${debug_c_compile_flags} ${arch_flag_C}")
endif()
set(release_c_compile_flags "-O2 ${arch_flag_C}")

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  if(APPLE)
    set(c_compile_flags "${c_compile_flags} -Wa,-q")
  endif(APPLE)
  # set(release_c_compile_flags "${release_c_compile_flags} -march=native")
elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
  # set(release_c_compile_flags "${release_c_compile_flags} -march=native")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(release_c_compile_flags "${release_c_compile_flags} -mp1")
elseif(CMAKE_C_COMPILER_ID MATCHES "PGI")
  # set(release_c_compile_flags "${release_c_compile_flags} -fastsse")
endif()

if(NOT CHARMM_C_FLAGS)
  set(CHARMM_C_FLAGS "${CMAKE_C_FLAGS} ${c_compile_flags}"
    CACHE STRING "C flags" FORCE)
endif()
set(CMAKE_C_FLAGS "${CHARMM_C_FLAGS}")

if(NOT CHARMM_C_FLAGS_DEBUG)
  set(CHARMM_C_FLAGS_DEBUG "${debug_c_compile_flags}"
    CACHE STRING "to use when CMAKE_BUILD_TYPE=Debug" FORCE)
endif()
set(CMAKE_C_FLAGS_DEBUG "${CHARMM_C_FLAGS_DEBUG}")

if(NOT CHARMM_C_FLAGS_RELEASE)
  set(CHARMM_C_FLAGS_RELEASE "${release_c_compile_flags}"
    CACHE STRING "to use when CMAKE_BUILD_TYPE=Release" FORCE)
endif()
set(CMAKE_C_FLAGS_RELEASE "${CHARMM_C_FLAGS_RELEASE}")

set(cxx_compile_flags)

if(needs_abi_zero)
  set(cxx_compile_flags "${cxx_compile_flags} -D_GLIBCXX_USE_CXX11_ABI=0")
endif()

set(debug_cxx_compile_flags "-g -O0")
if(domdec)
  set(debug_cxx_compile_flags
    "${debug_cxx_compile_flags} ${arch_flag_CXX}")
endif()
set(release_cxx_compile_flags "-O2 ${arch_flag_CXX}")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6.0"))
    set(cxx_compile_flags "${cxx_compile_flags} -fabi-version=0")
  endif()
  if(APPLE)
    set(cxx_compile_flags "${cxx_compile_flags} -Wa,-q")
  endif(APPLE)
  # set(cxx_compile_flags "${cxx_compile_flags} -march=native")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # set(cxx_compile_flags "${cxx_compile_flags} -march=native")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(cxx_compile_flags "${cxx_compile_flags} -mp1")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  # set(cxx_compile_flags "${cxx_compile_flags} -fastsse")
endif()

if(NOT CHARMM_CXX_FLAGS)
  set(CHARMM_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${cxx_compile_flags}"
    CACHE STRING "C++ flags" FORCE)
endif()
set(CMAKE_CXX_FLAGS "${CHARMM_CXX_FLAGS}")

if(NOT CHARMM_CXX_FLAGS_DEBUG)
  set(CHARMM_CXX_FLAGS_DEBUG "${debug_cxx_compile_flags}"
    CACHE STRING "to use when CMAKE_BUILD_TYPE=Debug" FORCE)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CHARMM_CXX_FLAGS_DEBUG}")

if(NOT CHARMM_CXX_FLAGS_RELEASE)
  set(CHARMM_CXX_FLAGS_RELEASE "${release_cxx_compile_flags}"
    CACHE STRING "to use when CMAKE_BUILD_TYPE=Release" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${CHARMM_CXX_FLAGS_RELEASE}")

if(NOT CMAKE_EXE_LINKER_FLAGS)
  set(CMAKE_EXE_LINKER_FLAGS "${link_flags}"
    CACHE STRING "for the linking phase" FORCE)
else()
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${link_flags}"
    CACHE STRING "for the linking phase" FORCE)
endif()

# END set compile and link flags
######

if(CMAKE_GENERATOR MATCHES Ninja)
  set_source_files_properties(
    source/util/keywords.F90
    source/charmm/iniall.F90
    PROPERTIES COMPILE_FLAGS "-I\"${CMAKE_CURRENT_BINARY_DIR}\"")
endif()

# begin ex EXPAND files
function(set_expand_inc inc_dir)
  if(CMAKE_GENERATOR MATCHES Ninja)
    set(EXPAND_PREFIX source/${inc_dir}/)
  else()
    set(EXPAND_PREFIX)
  endif()
  foreach(expand_file ${ARGN})
    set_source_files_properties("${EXPAND_PREFIX}${expand_file}"
      PROPERTIES COMPILE_FLAGS "-I\"${CMAKE_SOURCE_DIR}/source/${inc_dir}\"")
  endforeach()
endfunction(set_expand_inc)

set_expand_inc(nbonds
  colfft_util.F90
  colfft_func.F90
  colfft_kernel.F90
  enbaexp.F90
  enbonda.F90
  enbondg.F90
  nbondg.F90
  nbonda.F90)

set_expand_inc(pipf
  dpfimg.F90
  dpipf.F90
  epfinv.F90
  epipf.F90
  pfdyn.F90
  vpipf.F90)

set_expand_inc(manip fsshake_kernel.F90)
set_expand_inc(energy ediff.F90)
set_expand_inc(misc fctblock.F90 fctall.F90)
set_expand_inc(cff enbonda_cff.F90)
set_expand_inc(image nbondm.F90)
# end ex EXPAND files

if((CMAKE_C_COMPILER_ID MATCHES "GNU")
    OR (CMAKE_C_COMPILER_ID MATCHES "Intel")
    OR (CMAKE_C_COMPILER_ID MATCHES "PGI"))
  set_source_files_properties(
    "source/machdep/cstuff.c"
    "source/machdep/socket.c"
    "source/machdep/xdisp.c"
    PROPERTIES COMPILE_FLAGS "-DCHARMM_GNUC")
endif()

if(domdec)
  target_include_directories(charmm_fortran PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/source/domdec")
endif()

if(domdec_gpu AND OPENMP_FOUND
    AND (NOT (CMAKE_Fortran_COMPILER_ID MATCHES "Intel"))
    AND ((CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      OR (CMAKE_C_COMPILER_ID MATCHES "Intel")))
  target_link_libraries(charmm iomp5)
endif()

if(MPI_Fortran_FOUND)
  set_property(TARGET charmm charmm_fortran
    APPEND_STRING PROPERTY COMPILE_FLAGS
    "${MPI_Fortran_COMPILE_FLAGS} ")
  target_include_directories(charmm_fortran PRIVATE
    ${MPI_Fortran_INCLUDE_PATH})
  target_include_directories(charmm PRIVATE ${MPI_Fortran_INCLUDE_PATH})
  target_link_libraries(charmm ${MPI_Fortran_LIBRARIES})
endif(MPI_Fortran_FOUND)

if(EXAFMM_FOUND)
  target_link_libraries(charmm charmm_exafmm)
  if(TBB_FOUND)
    target_link_libraries(charmm tbb)
  endif(TBB_FOUND)
endif(EXAFMM_FOUND)

if(X11_FOUND)
  target_compile_definitions(charmm_fortran PRIVATE -Dxdisplay)
  target_include_directories(charmm_fortran PRIVATE ${X11_X11_INCLUDE_PATH})
  target_compile_definitions(charmm_c PRIVATE -Dxdisplay)
  target_include_directories(charmm_c PRIVATE ${X11_X11_INCLUDE_PATH})
  target_compile_definitions(charmm_cxx PRIVATE -Dxdisplay)
  target_include_directories(charmm_cxx PRIVATE ${X11_X11_INCLUDE_PATH})
  target_compile_definitions(charmm PRIVATE -Dxdisplay)
  target_include_directories(charmm PRIVATE ${X11_X11_INCLUDE_PATH})
  target_link_libraries(charmm ${X11_X11_LIB})
endif(X11_FOUND)

if(MKL_FOUND)
  target_include_directories(charmm_fortran PRIVATE ${MKL_INCLUDE_DIR}/fftw)
  target_include_directories(charmm PRIVATE ${MKL_INCLUDE_DIR}/fftw)

  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_target_properties(charmm
      PROPERTIES LINK_FLAGS -mkl=sequential)
  elseif(${MKL_GNU})
    target_link_libraries(charmm
      ${MKL_GNU}
      ${MKL_SEQUENTIAL}
      ${MKL_CORE})
  else()
    message(WARNING
      "MKL_GNU not found on this system."
      "Use the Intel Fortran compiler or FFTW.")
  endif()
    
  find_package(Threads REQUIRED)
  if(THREADS_HAVE_PTHREAD_ARG)
    set_property(TARGET charmm charmm_c charmm_cxx charmm_fortran
      APPEND_STRING PROPERTY COMPILE_FLAGS
      "-pthread ")
  endif()
  if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(charmm ${CMAKE_THREAD_LIBS_INIT})
  endif()
elseif(FFTW_FOUND)
  target_include_directories(charmm_fortran PRIVATE ${FFTW_INCLUDES})
  target_include_directories(charmm PRIVATE ${FFTW_INCLUDES})
  target_link_libraries(charmm ${FFTW_LIBRARIES})
  if(FFTWF_FOUND)
    target_link_libraries(charmm ${FFTWF_LIBRARIES})
  endif()
endif()

if(OPENMM_FOUND)
  add_subdirectory(source/openmm/plugins/MonteCarloBarostat2)
  add_subdirectory(source/openmm/plugins/gbsw)
  target_compile_definitions(charmm_cxx PRIVATE "-DKEY_OPENMM=1")
  target_include_directories(charmm_cxx PRIVATE ${OPENMM_INCLUDE_DIRS})
  target_include_directories(charmm_fortran PRIVATE
    ${OPENMM_INCLUDE_DIRS}
    source/openmm/plugins/MonteCarloBarostat2/wrappers
    source/openmm/plugins/gbsw/wrappers)
  target_link_libraries(charmm
    mcbarostat2
    OpenMMGBSW
    ${OPENMM_LIBRARIES})
  if(OPENMM_API_UPDATE)
    set_property(TARGET charmm_fortran
      APPEND_STRING
      PROPERTY COMPILE_FLAGS "-DOPENMM_API_UPDATE ")
  endif()
  # plugin_locs.f90 is generated in CMAKE_BINARY_DIR
  # omm_main.F90 includes plugin_locs.f90
  set_source_files_properties(source/openmm/omm_main.F90
    PROPERTIES COMPILE_FLAGS "-I\"${CMAKE_BINARY_DIR}\"")
endif()

if(OPENMP_FOUND)
  set_property(TARGET charmm_c
    APPEND_STRING PROPERTY COMPILE_FLAGS
    "${OpenMP_C_FLAGS} ")
  target_link_libraries(charmm_c ${OpenMP_C_FLAGS})

  set_property(TARGET charmm_cxx
    APPEND_STRING PROPERTY COMPILE_FLAGS
    "${OpenMP_CXX_FLAGS} ")
  target_link_libraries(charmm_cxx ${OpenMP_CXX_FLAGS})
elseif(AppleOpenMP_FOUND)
  set_property(TARGET charmm_c
    APPEND_STRING PROPERTY COMPILE_FLAGS
    "${APPLE_OPENMP_CFLAGS} ")
  target_link_libraries(charmm_c
    ${APPLE_OPENMP_CFLAGS} ${APPLE_OPENMP_LIBRARY})

  set_property(TARGET charmm_cxx
    APPEND_STRING PROPERTY COMPILE_FLAGS
    "${APPLE_OPENMP_CXXFLAGS} ")
  target_link_libraries(charmm_cxx
    ${APPLE_OPENMP_CXXFLAGS} ${APPLE_OPENMP_LIBRARY})
endif()

if(OPENMP_FORTRAN_FOUND)
  set_property(SOURCE ${charmm_fortran_files}
    APPEND_STRING PROPERTY COMPILE_FLAGS
    " ${OpenMP_Fortran_FLAGS} ")
  set_property(TARGET charmm
    APPEND_STRING PROPERTY COMPILE_FLAGS
    " ${OpenMP_Fortran_FLAGS} ")
  target_link_libraries(charmm_fortran ${OpenMP_Fortran_FLAGS})
  target_link_libraries(charmm ${OpenMP_Fortran_FLAGS})
endif()

if(domdec_gpu)
  target_link_libraries(charmm
    ${CUDA_LIBRARIES}
    ${CUDA_CUFFT_LIBRARIES})
  target_link_libraries(charmm
    ${NSIGHT_LIBRARIES})
endif(domdec_gpu)

install(CODE "
  message(\"     Installing into ${CMAKE_INSTALL_PREFIX}\") ")
install(TARGETS charmm DESTINATION bin)

if(qchem)
  install(TARGETS qchem DESTINATION bin)
endif(qchem)

if(NOT in_place_install)
  include(tool/cmake/InstallSymlink.cmake)
  InstallSymlink("${PROJECT_SOURCE_DIR}/doc"
    "${CMAKE_INSTALL_PREFIX}/doc")
  InstallSymlink("${PROJECT_SOURCE_DIR}/source"
    "${CMAKE_INSTALL_PREFIX}/source")
  InstallSymlink("${PROJECT_SOURCE_DIR}/tool"
    "${CMAKE_INSTALL_PREFIX}/tool")

  install(CODE "
  execute_process(COMMAND \"${CMAKE_COMMAND}\" -E make_directory
                                \"${CMAKE_INSTALL_PREFIX}/test\") ")

  file(GLOB test_subdirs "${PROJECT_SOURCE_DIR}/test/*")
  foreach(subdir ${test_subdirs})
    get_filename_component (new_link_name ${subdir} NAME)
    InstallSymlink("${subdir}"
      "${CMAKE_INSTALL_PREFIX}/test/${new_link_name}")
  endforeach()
endif(NOT in_place_install)
