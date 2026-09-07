#Specify minimum version of Cmake
cmake_minimum_required(VERSION 2.6...3.21)
enable_language (Fortran)
#Project name
project(Dust)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#BUG in CMake 3.7.x a policy error in find OpenMP raises unwanted warings
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()

#CMake 3.15.2 raises blocking warnings unless setting:
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

#Set folders
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)


# =========== Options ========================
SET(PRECISION "DOUBLE" CACHE STRING "Computation precision")
SET_PROPERTY(CACHE PRECISION PROPERTY STRINGS "DOUBLE" "SINGLE")
IF(PRECISION STREQUAL "SINGLE")
  ADD_DEFINITIONS(-DDUST_PRECISION=1)
ELSEIF(PRECISION STREQUAL "DOUBLE")
  ADD_DEFINITIONS(-DDUST_PRECISION=2)
ELSE()
  MESSAGE(ERROR "Unknown computational precision")
ENDIF()

SET(KERNEL "ROSENHEAD" CACHE STRING "Kernel to be used for the FMM")
SET_PROPERTY(CACHE KERNEL PROPERTY STRINGS "ROSENHEAD" "HOA")
IF(KERNEL STREQUAL "ROSENHEAD")
  ADD_DEFINITIONS(-DDUST_KERNEL=1)
ELSEIF(KERNEL STREQUAL "HOA")
  ADD_DEFINITIONS(-DDUST_KERNEL=2)
ELSE()
  MESSAGE(ERROR "Unknown computational precision")
ENDIF()

SET(WITH_PRECICE "NO" CACHE STRING "Compile with PreCICE")
SET_PROPERTY(CACHE WITH_PRECICE PROPERTY STRINGS "NO" "YES")
IF(WITH_PRECICE STREQUAL "YES")
  ADD_DEFINITIONS(-DUSE_PRECICE)
ELSEIF(WITH_PRECICE STREQUAL "NO")
  # Do nothing
ELSE()
  MESSAGE(ERROR "Unknown choice for WITH_PRECICE")
ENDIF()


# =========== COMPILER ========================
OPTION (USE_OpenMP "Use OpenMP to enamble <omp.h>" ON)

# Find OpenMP
if(APPLE AND USE_OpenMP)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

if(USE_OpenMP)
  find_package(OpenMP REQUIRED)
endif(USE_OpenMP)


find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp ")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-function -Wno-character-truncation ")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -g -fopenmp")
  set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O3 -g -fopenmp")
  #set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -ggdb -g -O0 -pedantic -Wall -Wconversion -Wconversion-extra -fcheck=all -fcheck=bounds -std=f2008ts -fall-intrinsics -finit-real=snan -ffpe-trap=invalid -finit-integer=-666 -Wintrinsics-std -fbacktrace")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -ggdb -g -O0 -fopenmp -pedantic -Wall -Wconversion -Wconversion-extra -fcheck=all -fcheck=bounds  -fall-intrinsics -finit-real=snan -ffpe-trap=invalid -finit-integer=-666 -Wintrinsics-std -fbacktrace -Wno-unused-dummy-argument")
ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel") 
  #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -standard-semantics")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -heap-arrays 10")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 -g -xHost")
  set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS}  -O2 -g -xHost")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -g -O0 -stand f08 -check all -check noarg_temp_created -warn all -fpe0 -traceback -ftrapuv -fp-stack-check -assume norealloc_lhs")
ELSE()
  message(ERROR "Unknown compiler")
ENDIF()


#> Lapack -------------------------------------------------------
OPTION(DUST_MKL "Use MKL libraries instead of BLAS/LAPACK" OFF)
IF (DUST_MKL)
#  SET(BLA_VENDOR "Intel") #search only for intel libraries
  find_package(MKL REQUIRED)
  LIST(APPEND incldirs ${MKL_INCLUDE_DIRS}) 
  #LIST(APPEND linkedlibs ${MKL_LIBRARIES} pthread dl)
  LIST(APPEND linkedlibs "-Wl,--start-group" ${MKL_LIBRARIES} "-Wl,--end-group" pthread dl)
  IF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel") 
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -mkl")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} -mkl")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -mkl")
  ENDIF() 
ELSE()
#  SET(BLA_VENDOR "OpenBLAS")
  find_package(LAPACK REQUIRED)
  LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
ENDIF()

#> HDF5 ---------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C Fortran)
LIST(APPEND linkedlibs ${HDF5_LIBRARIES})

#> CGNS ---------------------------------------------------------
# Add CGNS installation (for the moment could not make it look for it autonomously)
SET(CGNS_INC "/usr/include" CACHE STRING "Include path for cgns headers") 
SET(CGNS_LIB "/usr/lib" CACHE STRING "Path of the cgns library") 
LIST(APPEND linkedlibs cgns)

#> PreCICE ------------------------------------------------------
IF(WITH_PRECICE STREQUAL "YES")
  SET(PRECICE_INC "/usr/include"   CACHE STRING "Include path for PreCICE fortran interface") 
  SET(PRECICE_LIB "/usr/local/lib" CACHE STRING "Path of the PreCICE library") 
  LIST(APPEND linkedlibs precice)
ENDIF()

#set the includes
include_directories(${incldirs} ${PRECICE_INC} ${CGNS_INC} ${HDF5_INCLUDE_DIRS})

#set the link path
link_directories(${PRECICE_LIB} ${CGNS_LIB} ${HDF5_LIBRARY_DIRS})

#Generate the version
execute_process(COMMAND
  git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND
  git describe
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_VERS
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND
  git symbolic-ref -q --short HEAD
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_BRAN
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

configure_file( "${PROJECT_SOURCE_DIR}/src/general/mod_version.f90.in" 
                "${PROJECT_SOURCE_DIR}/src/general/mod_version.f90")


#add files
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/geo/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/general/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/aeroel/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/linsys/*.f90
                          ${PROJECT_SOURCE_DIR}/src/io/*.f90
                          ${PROJECT_SOURCE_DIR}/src/post/*.f90
                          ${PROJECT_SOURCE_DIR}/src/options/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/viscosity/*.f90
                          ${PROJECT_SOURCE_DIR}/src/fmm/*.f90)
IF (WITH_PRECICE STREQUAL "YES")
  file(GLOB sources_precice ${PROJECT_SOURCE_DIR}/src/precice/*.f90)
  list(APPEND sources ${sources_precice} )
ENDIF()


#create a library
add_library(dustlib STATIC ${sources})
# Enforce where the generated Fortran .mod files must be placed during compilation
set_target_properties(dustlib PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(dustlib PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

#add executables to be compiled
add_executable(dust  ${PROJECT_SOURCE_DIR}/src/dust.f90)
target_link_libraries(dust PRIVATE dustlib ${linkedlibs})

add_executable(dust_pre  ${PROJECT_SOURCE_DIR}/src/dust_pre.f90)
target_link_libraries(dust_pre PRIVATE dustlib ${linkedlibs})

add_executable(dust_post  ${PROJECT_SOURCE_DIR}/src/dust_post.f90)
target_link_libraries(dust_post PRIVATE dustlib ${linkedlibs})

# only for testing pipeline! 
install(TARGETS dust dust_pre dust_post DESTINATION bin)