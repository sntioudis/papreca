cmake_minimum_required(VERSION 3.16)
project(papreca)

#Set minimum C++ required and CXX compiler options
set(CMAKE_CXX_STANDARD 11)
if(CMAKE_CXX_STANDARD LESS 11)
  message(FATAL_ERROR "C++ standard must be set to at least 11")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Use compiler extensions")

if(BUILD_DEBUG)
	set(CMAKE_BUILD_TYPE Debug)
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
else()
	set(CMAKE_BUILD_TYPE Release)
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

#LAMMPS configuration
#LAMMPS library configuration
if(NOT DEFINED LAMMPS_LIB_DIR)
    message(FATAL_ERROR "LAMMPS library (i.e., liblammps.a, typically located in the build folder) directory not provided. Please specify the path in your cmake command using -DLAMMPS_LIB_DIR=/path/to/lammps/library/")
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

#Setup rules for PAPRECA static library
set(PAPRECA_LIB_DIR ../../source/libraries/PAPRECA/)
file(GLOB PAPRECA_LIB_SOURCES ${PAPRECA_LIB_DIR}/*.cpp)
add_library(PAPRECA_LIBRARY STATIC ${PAPRECA_LIB_SOURCES})
target_link_libraries(PAPRECA_LIBRARY PRIVATE ${LAMMPS_LIBRARY})

#Link papreca against PAPRECA and LAMMPS libraries. Also, include the necessary header directories
add_executable(papreca ../../source/papreca.cpp)
target_link_libraries(papreca PRIVATE ${LAMMPS_LIBRARY})

target_link_libraries(papreca PRIVATE PAPRECA_LIBRARY)
target_include_directories(papreca PRIVATE ${PAPRECA_LIB_DIR})


message(STATUS "Configuring additional packages..." )
message( WARNING "PAPRECA does not use OpenMP, JPEG, and PNG packages. However, the LAMMPS library uses them and therefore, we need to compile against those libraries to avoid linkage errors. Depending on how YOUR LAMMPS library was built (i.e., based on the LAMMPS packages you enabled) you might need to enable additional libraries in the PAPRECA CMakeLists.txt file.")

find_package(MPI REQUIRED)
if(MPI_FOUND)
	# do not include the (obsolete) MPI C++ bindings which makes
	# for leaner object files and avoids namespace conflicts
	set(MPI_CXX_SKIP_MPICXX TRUE)
	target_link_libraries(papreca PUBLIC MPI::MPI_CXX)
	target_link_libraries(PAPRECA_LIBRARY PUBLIC MPI::MPI_CXX)
else()
	message(FATAL_ERROR "MPI not found. Please make sure MPI is installed.")
endif()

find_package(OpenMP COMPONENTS CXX)
if(OpenMP_FOUND)
	target_link_libraries(papreca PRIVATE OpenMP::OpenMP_CXX)
	target_link_libraries(PAPRECA_LIBRARY PRIVATE OpenMP::OpenMP_CXX)
else()
	message(FATAL_ERROR "Cannot find the OpenMP package required by the LAMMPS library.")
endif()

# tweak jpeg library names to avoid linker errors with MinGW cross-compilation
set(JPEG_NAMES libjpeg libjpeg-62)
find_package(JPEG QUIET)
option(WITH_JPEG "Enable JPEG support" ${JPEG_FOUND})
if(WITH_JPEG)
	find_package(JPEG REQUIRED)
	target_link_libraries(papreca PRIVATE JPEG::JPEG)
	target_link_libraries(PAPRECA_LIBRARY PRIVATE JPEG::JPEG)
endif()

find_package(PNG QUIET)
find_package(ZLIB QUIET)
if(PNG_FOUND AND ZLIB_FOUND)
	option(WITH_PNG "Enable PNG support" ON)
else()
	option(WITH_PNG "Enable PNG support" OFF)
endif()
if(WITH_PNG)
	find_package(PNG REQUIRED)
	find_package(ZLIB REQUIRED)
	target_link_libraries(papreca PRIVATE PNG::PNG ZLIB::ZLIB)
	target_link_libraries(PAPRECA_LIBRARY PRIVATE PNG::PNG ZLIB::ZLIB)
endif()
