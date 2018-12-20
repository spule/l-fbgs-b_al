# cmake version, project name, language
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(l-fbgs-b_al Fortran)

# set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# if build type not specified, default to release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "release")
endif()

# handle integer size
if(INT)
  message(STATUS "Configuring build for ${INT}-bit integers")
  add_definitions(-DINT${INT})
else()
  message(STATUS "Configuring build for 32-bit integers")
  add_definitions(-DINT32)
endif()

# handle real size
if(REAL)
  message(STATUS "Configuring build for ${REAL}-bit reals")
  add_definitions(-DREAL${REAL})
else()
  message(STATUS "Configuring build for 32-bit reals")
  add_definitions(-DREAL32)
endif()

# compiler flags for gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)

  if(SERIAL)
    message(STATUS "Configuring to build with -fcoarrays-single")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcoarray=single")
  endif()
    
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -C -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math")

endif()

# compiler flags for ifort
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -coarray=shared -fpp -assume byterecl,realloc_lhs -heap-arrays")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -C -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

# library to archive (liblbfgsb_al.a)
add_library(lbfgsb_al src/lib/lbfgsb.f90 src/lib/lbfgsb_al.f90)

# tests
enable_testing()
foreach(execid driver1 driver2 driver3)
  add_executable(test_${execid} src/tests/test_${execid}.f90)
  target_link_libraries(test_${execid} lbfgsb_al)
  add_test(test_${execid} bin/test_${execid})
endforeach()
  
#foreach(execid mnist simple sine)
#  add_executable(example_${execid} src/tests/example_${execid}.f90)
#  target_link_libraries(example_${execid} lbfgsb_al)
#  add_test(example_${execid} bin/example_${execid})
#endforeach()
