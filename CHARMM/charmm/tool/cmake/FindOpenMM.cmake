# required returns
set(OPENMM_FOUND)
set(OPENMM_INCLUDE_DIRS)
set(OPENMM_LIBRARIES)
set(OPENMM_PLUGIN_DIR)

# components: CPU, CUDA, OPENCL
set(OPENMM_CPU_PLUGIN)
set(OPENMM_CUDA_PLUGIN)
set(OPENMM_OPENCL_PLUGIN)

# used internally
set(OPENMM_INCLUDE_DIR)
set(OPENMM_LIBRARY)
set(OPENMM_LIBRARY_DIR)

if(NOT OpenMM_FIND_COMPONENTS)
    # Assume they don't want cuda or opencl
    # the rest should be more portable
    set(OpenMM_FIND_COMPONENTS CPU)
endif()

find_path(
  OPENMM_INCLUDE_DIR
  NAMES OpenMMFortranModule.f90
  PATHS
    ENV OPENMM_INCLUDE_PATH
    ENV OPENMM_HOME
    ENV OPENMM_DIR
    $ENV{OPENMM_PLUGIN_DIR}/../..
    /usr/local/openmm
  PATH_SUFFIXES "include"
  NO_DEFAULT_PATH)

find_library(
  OPENMM_LIBRARY
  NAMES OpenMM
  PATHS
    ENV OPENMM_LIBRARY_PATH
    ENV OPENMM_HOME
    ENV OPENMM_DIR
    $ENV{OPENMM_PLUGIN_DIR}/../..
    /usr/local/openmm
  PATH_SUFFIXES "lib"
  NO_DEFAULT_PATH)

set(OPENMM_INCLUDE_DIRS "${OPENMM_INCLUDE_DIR}")
set(OPENMM_LIBRARIES "${OPENMM_LIBRARY}")
get_filename_component(OPENMM_LIBRARY_DIR "${OPENMM_LIBRARY}" PATH)
set(OPENMM_PLUGIN_DIR "${OPENMM_LIBRARY_DIR}/plugins")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  OPENMM
  DEFAULT_MSG
  OPENMM_INCLUDE_DIRS
  OPENMM_LIBRARIES
  OPENMM_PLUGIN_DIR)

foreach(component ${OpenMM_FIND_COMPONENTS})
  set(lib_suffix "${component}")
  if(component STREQUAL "OPENCL")
    set(lib_suffix "OpenCL")
  endif()

  find_library(
    OPENMM_${component}_PLUGIN
    NAMES OpenMM${lib_suffix}
    PATHS
      ENV OPENMM_PLUGIN_DIR
      ${OPENMM_PLUGIN_PATH}
      /usr/local/openmm/lib/plugins
    NO_DEFAULT_PATH)

  find_package_handle_standard_args(
    OPENMM_${component}
    DEFAULT_MSG
    OPENMM_${component}_PLUGIN)
endforeach()
