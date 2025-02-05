cmake_minimum_required(VERSION 3.16)
project(source_tests)

#Set minimum C++ required and CXX compiler options
set(CMAKE_CXX_STANDARD 17)
if(CMAKE_CXX_STANDARD LESS 17)
  message(FATAL_ERROR "C++ standard must be set to at least 117")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Use compiler extensions")


#Enable debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g")

#LAMMPS configuration
#LAMMPS (SHARED) library configuration
if(NOT DEFINED LAMMPS_LIB_DIR)
    message(FATAL_ERROR "LAMMPS library (i.e., liblammps.so, typically located in the build folder) directory not provided. Please specify the path in your cmake command using -DLAMMPS_LIB_DIR=/path/to/lammps/shared/library/")
endif()

find_library(LAMMPS_LIBRARY
    NAMES lammps
    HINTS ${LAMMPS_LIB_DIR}
)

if(LAMMPS_LIBRARY)
    message(STATUS "Found LAMMPS library: ${LAMMPS_LIBRARY}")
else()
    message(FATAL_ERROR "LAMMPS library not found in directory: ${LAMMPS_LIB_DIR}")
endif()

#LAMMPS includes configuration
if(NOT DEFINED LAMMPS_SRC_DIR)
    message(FATAL_ERROR "LAMMPS source directory not provided. Please specify the path in your cmake command using -DLAMMPS_SRC_DIR=/path/to/lammps/src/")
endif()

include_directories(${LAMMPS_SRC_DIR})

#PAPRECA configuration
if(NOT DEFINED PAPRECA_SRC_DIR)
    message(FATAL_ERROR "PAPRECA source directory not provided. Please specify the path in your cmake command using -DPAPRECA_SRC_DIR=/path/to/PAPRECA/source/")
endif()

set(PAPRECA_LIB_DIR ${PAPRECA_SRC_DIR}/libraries/PAPRECA/)
include_directories( ${PAPRECA_LIB_DIR} )
file(GLOB PAPRECA_LIB_SOURCES ${PAPRECA_LIB_DIR}/*.cpp)

#Setup rules for PAPRECA static library
file(GLOB PAPRECA_LIB_SOURCES ${PAPRECA_LIB_DIR}/*.cpp)
add_library(PAPRECA_LIBRARY STATIC ${PAPRECA_LIB_SOURCES})
target_link_libraries(PAPRECA_LIBRARY PRIVATE ${LAMMPS_LIBRARY})

#Link the source_tests executable against the PAPRECA, and LAMMPS libraries. Also, include the necessary header directories
add_executable(source_tests ../source_tests.cpp)
target_link_libraries(source_tests PRIVATE ${LAMMPS_LIBRARY})
target_link_libraries(source_tests PRIVATE PAPRECA_LIBRARY)


message(STATUS "Configuring additional packages..." )

find_package(MPI REQUIRED)
if(MPI_FOUND)
	# do not include the (obsolete) MPI C++ bindings which makes
	# for leaner object files and avoids namespace conflicts
	set(MPI_CXX_SKIP_MPICXX TRUE)
	target_link_libraries(source_tests PUBLIC MPI::MPI_CXX)
	target_link_libraries(PAPRECA_LIBRARY PUBLIC MPI::MPI_CXX)
else()
	message(FATAL_ERROR "MPI not found. Please make sure MPI is correctly installed/configured.")
endif()
