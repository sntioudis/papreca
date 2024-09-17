#!/bin/bash

#User defined installation packages for LAMMPS
package_args=" -DPKG_EXTRA-DUMP=on -DPKG_MOLECULE=on -DPKG_RIGID=on -DPKG_QEQ=on -DPKG_REAXFF=on -DPKG_REPLICA=on -DPKG_KSPACE=on -DPKG_KOKKOS=on -DPKG_MEAM=on"


#Create build directories
cd ../../
mkdir build; cd build

#Pull LAMMPS from LAMMPS repository
git clone --depth 1 --branch patch_29Aug2024 https://github.com/lammps/lammps.git LAMMPS #clone LAMMPS with tag patch_17Apr2024

#Copy source files for fix PAPRECA
cd LAMMPS
cp ../../source/libraries/LAMMPS/fix_papreca.h ./src/
cp ../../source/libraries/LAMMPS/fix_papreca.cpp ./src/


#Compile LAMMPS with package options and fix PAPRECA
mkdir build; cd build
cmake -DBUILD_LIB=on -DBUILD_SHARED_LIBS=on $package_args ../cmake
make -j

#Delete unwanted installation files (keep only build and src folders required for PAPRECA installation)
cd ..
rm -rf ./bench/ ./cmake/ ./doc/ ./examples/ ./fortran/ ./.git/ ./.github/ ./lib/ ./potentials/ ./python/ ./tools ./unittest/ ./CITATION.cff ./.gitattributes ./.gitignore ./.lgtm.yml ./LICENSE ./README ./SECURITY.md
