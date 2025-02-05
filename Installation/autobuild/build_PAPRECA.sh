#!/bin/bash

#Read number of cores
read -p "Clean previous build (yes/no)?: " clean_build
# Convert input to lowercase to handle cases like "Yes", "YES", etc.
clean_build=$(echo "$clean_build" | tr '[:upper:]' '[:lower:]')

#Clean build if user selects it
if [[ "$clean_build" == "yes" || "$clean_build" == "y" ]]; then
	rm -rf ../../build/
fi


#Install LAMMPS with required packages and fix papreca
bash build_LAMMPS.sh

#Navigate to build folder
cd ../../build

#Get library and source directories
lammps_lib_path="./LAMMPS/build/"
lammps_lib_path=$(realpath "$lammps_lib_path") #Resolve real paths to avoid directory conflicts
lammps_src_path="./LAMMPS/src/"
lammps_src_path=$(realpath "$lammps_src_path")

mkdir PAPRECA; cd PAPRECA

#Install PAPRECA on top of LAMMPS
cmake ../../Installation/CMake -DLAMMPS_SRC_DIR="${lammps_src_path}" -DLAMMPS_LIB_DIR="${lammps_lib_path}"
make -j

