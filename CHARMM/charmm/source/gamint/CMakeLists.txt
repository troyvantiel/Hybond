cmake_minimum_required(VERSION 2.8) # same as CHARMM ? MUST TEST !!!!

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLINUX64") # just for zunix.c

set(fortcomp "gfortran") # default
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(fortcomp "ifort")
endif()
if (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(fortcomp "pgfortran")
endif()

set(openmpcomp "false")
if (OpenMP_Fortran_FLAGS MATCHES "openmp")
  set(openmpcomp "false") # currently not supported due to new code in GAMESS
endif()

# this could replace gamess's config & install.info

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/charmm_install.info
     "set chmsrc=${CMAKE_CURRENT_SOURCE_DIR}\n"
     "setenv FC ${CMAKE_Fortran_COMPILER}\n"
     "set GMS_OPENMP=${openmpcomp}\n"
     "set GMS_FORTRAN=${fortcomp}\n")

function(process_gmssrc_file from_file to_file)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${to_file}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../../tool/gmscomp ${CMAKE_CURRENT_SOURCE_DIR}/${from_file}
    )
  set_source_files_properties(${to_file} PROPERTIES GENERATED TRUE)
endfunction(process_gmssrc_file)

# modules needs to go first!!!
# this is probably OK, or it could be made into one line with the list??
process_gmssrc_file(gamess/modmnfun.src modmnfun.o)
process_gmssrc_file(gamess/mpcdatpm6.src mpcdatpm6.o)
process_gmssrc_file(gamess/mpchbond.src mpchbond.o)
process_gmssrc_file(gamess/secor.src secor.o)
process_gmssrc_file(gamess/prec.src prec.o)
process_gmssrc_file(gamess/params.src params.o)
process_gmssrc_file(gamess/grd2_consts.src grd2_consts.o)
process_gmssrc_file(gamess/mx_limits.src mx_limits.o)
process_gmssrc_file(gamess/blkint.src blkint.o)
process_gmssrc_file(gamess/ompmod_tools.src ompmod_tools.o)
process_gmssrc_file(gamess/ompmod.src ompmod.o)
process_gmssrc_file(gamess/ompgrd2.src ompgrd2.o)
process_gmssrc_file(gamess/grd2a.src grd2a.o)
#process_gmssrc_file(gamess/rimp2grd.src rimp2grd.o)
process_gmssrc_file(gamess/ryspol.src ryspol.o)
process_gmssrc_file(gamess/efpmodule.src efpmodule.o)
process_gmssrc_file(gamess/rmddat.src rmddat.o)
#process_gmssrc_file(gamess/mod_gauss_hermite.src mod_gauss_hermite.o)
#process_gmssrc_file(gamess/mod_nosp_basis.src mod_nosp_basis.o)
#process_gmssrc_file(gamess/mod_shell_tools.src mod_shell_tools.o)
#process_gmssrc_file(gamess/mod_1e_primitives.src mod_1e_primitives.o)


# make a list of all src files in gamess/ folder
file(GLOB GmsSource
  RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/gamess"
  "${CMAKE_CURRENT_SOURCE_DIR}/gamess/*.src")

file(GLOB CCSD3Source
  RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/gamess"
  "${CMAKE_CURRENT_SOURCE_DIR}/gamess/ccsd3*.src")

# do we care if we double the above files???
list(REMOVE_ITEM GmsSource modmnfun.src)
list(REMOVE_ITEM GmsSource params.src)
list(REMOVE_ITEM GmsSource efpmodule.src)
list(REMOVE_ITEM GmsSource mpcdatpm6.src)
list(REMOVE_ITEM GmsSource mpchbond.src)
list(REMOVE_ITEM GmsSource rmddat.src)
list(REMOVE_ITEM GmsSource secor.src)
list(REMOVE_ITEM GmsSource prec.src)
list(REMOVE_ITEM GmsSource mx_limits.src)
list(REMOVE_ITEM GmsSource grd2_consts.src)
list(REMOVE_ITEM GmsSource blkint.src)
list(REMOVE_ITEM GmsSource ompmod_tools.src)
list(REMOVE_ITEM GmsSource ompmod.src)
list(REMOVE_ITEM GmsSource ompgrd2.src)
list(REMOVE_ITEM GmsSource grd2a.src)
#list(REMOVE_ITEM GmsSource rimp2grd.src)
list(REMOVE_ITEM GmsSource ryspol.src)
#list(REMOVE_ITEM GmsSource mod_gauss_hermite.src)
#list(REMOVE_ITEM GmsSource mod_nosp_basis.src)
#list(REMOVE_ITEM GmsSource mod_shell_tools.src)
#list(REMOVE_ITEM GmsSource mod_1e_primitives.src)
# later make this optional:
list(REMOVE_ITEM GmsSource ${CCSD3Source})
# we don't want these two yet:
list(REMOVE_ITEM GmsSource vector.src)
list(REMOVE_ITEM GmsSource ga.src)

# file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/debug_lists.gmssource  ${GmsSource})



FOREACH(GmsFile ${GmsSource})

  STRING(REGEX REPLACE ".src\$" "" GmsO "${GmsFile}")
  process_gmssrc_file(source/gamess/${GmsFile} ${GmsO}.o)
  SET(GmsObjs ${GmsObjs} "${GmsO}.o")
  
ENDFOREACH(GmsFile)

# this is just for the correct order - use depend !???
# we need zunix.c & gukini_mod dependencies...
add_library(gmsmod params.o secor.o mpchbond.o rmddat.o efpmodule.o mpcdatpm6.o
  prec.o mx_limits.o grd2_consts.o blkint.o ompmod_tools.o ompmod.o ompgrd2.o
  modmnfun.o grd2a.o rimp2omp.o ryspol.o
  )

add_library(gmsc ${CMAKE_CURRENT_SOURCE_DIR}/gamess/zunix.c)
# this is the real library:
add_library(gamess ${GmsObjs} mpchbond.o rmddat.o mpcdatpm6.o params.o secor.o efpmodule.o
  prec.o mx_limits.o grd2_consts.o blkint.o ompmod_tools.o ompmod.o ompgrd2.o modmnfun.o
  grd2a.o rimp2omp.o ryspol.o
  )

set_property(TARGET gmsmod PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET gamess PROPERTY LINKER_LANGUAGE Fortran )
