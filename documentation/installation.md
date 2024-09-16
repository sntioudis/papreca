# Download and Installation

\anchor installation

\section download PAPRECA repository and download

The source code of %PAPRECA is stored and can be downloaded/pulled from this [GitHub repository](https://github.com/sntioudis/papreca).

To clone the latest %PAPRECA repository on your machine execute the following commands:

```bash
git clone -b release https://github.com/sntioudis/papreca.git mypapreca #Add the latest PAPRECA repository to a folder named mypapreca
```

\section beforeInstall Read before you install: Important notes and prerequisites

&bull; Copy fix_papreca.h and fix_papreca.cpp files from the %PAPRECA source files (located in ./source/libraries/LAMMPS/ of the pulled %PAPRECA directory) to the LAMMPS source files (located in ./src/ of the pulled LAMMPS directory).
This will install LAMMPS with the fix papreca command option which is necessary to run %PAPRECA (see \ref FIX_papreca for more information).

&bull; [Build LAMMPS](https://docs.lammps.org/Install.html) as a **shared** library (with any optional LAMMPS packages that you intend to use) before you build %PAPRECA. The minimum working LAMMPS version (tag) is: patch_15Sep2022. Version patch_15Sep2022 and any later release version have been tested with %PAPRECA and are expected to function properly. Nevertheless, users are advised to use the latest working LAMMPS version (patch_29Aug2024).

The following snippet demonstrates (briefly) how LAMMPS can be built as a shared library with a few optional packages:

```bash
git clone --depth 1 --branch patch_29Aug2024 https://github.com/lammps/lammps.git mylammps #clone LAMMPS with tag patch_17Apr2024 to a folder named mylammps
cd mylammps
mkdir build; cd build
cmake -DPKG_EXTRA-DUMP=on -DPKG_MOLECULE=on -DPKG_RIGID=on -DPKG_QEQ=on -DPKG_REAXFF=on -DBUILD_LIB=on -DBUILD_SHARED_LIBS=on ../cmake #Configure LAMMPS, build with some optional package, and enable lib building as well as shared libraries
cmake --build . #make -j N (where N is the number of processes) can also be used for multicore installation
```

> **Note 1:**
> To run all the examples in the ./Examples/ folder (see \ref examples) you must build your LAMMPS library with the following packages: **EXTRA-DUMP**, **MOLECULE**, **RIGID**, and **QEQ**.


&bull; The current version of %PAPRECA (1.0) runs only on LINUX-based systems. Cross-platform compatibility will be available in future versions. At the moment, non-LINUX users can run %PAPRECA on a virtual machine (e.g., [VirtualBox](https://www.virtualbox.org/) or Windows Subsystem for Linux).

&bull; An MPI/C++ compiler (e.g., mpicxx or mpiCC) that is at least compatible with the C++11 standard is required to build %PAPRECA.

&bull; After configuring the LAMMPS build, kindly ensure that the MPI package has been correctly configured. If the MPI package is not configured, cmake will build a non-MPI LAMMPS version which will lead to %PAPRECA linkage errors. If the installed MPI compiler in your system is not detected by cmake, you should add the following argument to your cmake configuration command: -DCMAKE_CXX_COMPILER=/path/to/your/cxx/compiler. Hint: the relevant cxx compiler path can be obtained by running: 'which compiler_name' (e.g., 'which mpicxx'). It is advised that the same flag (i.e., -DCMAKE_CXX_COMPILER) is also used for the configuration of LAMMPS.

&bull; Exactly the same package (e.g., MPI protocol, MPI/C++ compiler, openMP) versions have to be used when building LAMMPS and %PAPRECA. You will probably encounter linking errors if there is a package version mismatch.

\section build Building PAPRECA

After building LAMMPS as a library, %PAPRECA can be built using a Traditional Make approach or a [CMake](https://cmake.org/) approach. However, mixed Traditional Make/CMake builds should be avoided! If you installed %PAPRECA using a traditional Make approach you should remove all associated libraries, dependency files, and executables before re-building with CMake (and vice versa).

\subsection autobuild %PAPRECA autobuild

The quickest way to install %PAPRECA is by using the autobuild script located in ./Installation/autobuild. Running the autobuild script can be as simple as executing the following bash command form the ./Installation/autobuild/ directory of the %PAPRECA project:

```bash
bash build_PAPRECA.sh
```

The autobuild script pulls LAMMPS from the official Github repository, builds LAMMPS with the fix papreca command, and finally, builds %PAPRECA. 

> **Note:**
> At the moment, building with autobuild will generate a LAMMPS library that includes the following packages: **EXTRA-DUMP**, **MOLECULE**, **RIGID**, **QEQ**, **REAXFF**, **REPLICA**, **KSPACE**, **KOKKOS**, **MEAM**. Modify the **package_args** variable of the build_LAMMPS.sh script (located in the ./Installation/autobuild/ directory of the %PAPRECA project) to build with different packages, if necessary.


\subsection cmake CMake

Building both LAMMPS and %PAPRECA with CMake is generally recommended since it does not require the manual specification of installation packages (e.g., MPI/C++ compiler). However, some modifications
to the CMakelists.txt file might have to be made if your machine includes many different installation package versions.

To build %PAPRECA with CMake you should execute the following commands in the cloned %PAPRECA repository:

```bash
mkdir build; cd build
cmake ../Installation/CMake -DLAMMPS_SRC_DIR=/path/to/LAMMPS/source -DLAMMPS_LIB_DIR=/path/to/LAMMPS/library #Replace paths with YOUR source (./src LAMMPS directory) and library (typically ./build LAMMPS directory) paths.
cmake --build . #make -j N (where N is the number of processes) can also be used for multicore installation
```

Where, you must replace the relevant paths with your source (./src LAMMPS directory) and library (typically ./build LAMMPS directory) LAMMPS . 

> **Note:**
> To build %PAPRECA with debug symbols you must also include the following flag: "-DBUILD_DEBUG=ON". If such flag is ignored, CMake will build %PAPRECA without debug symbols and with the following C++ compiler optimization flag: -O3.

\subsection traditionalMake Traditional Make

An example traditional "MakeFile" together with a "MakeConfig" file are located in this project directory: Installation/Traditional Make/. Before building the project the user should modify the "MakeConfig" file.
Specifically, the following paths must be updated:

&bull; The "BUILDDIR" variable should point at the %PAPRECA build directory. Make sure the directory exists before you build.

&bull; The "LAMMPS" variable should point at the main-parent LAMMPS directory. Note that, this directory neither the ./src or ./build paths (i.e., your input directory is the folder containing both build and src folders).

After configuring the relevant paths in "MakeConfig" the "traditionalMakeCompileLIBandSRC.sh" script (included in the Installation/Traditional Make/ directory) can be executed to build the %PAPRECA libraries and the %papreca executable (in the specified build folder):

```bash
bash traditionalMakeCompileLIBandSRC.sh
```

This is what the "traditionalMakeCompileLIBandSRC.sh" script does in detail:

```bash
#Install PAPRECA library
cd ../../source/libraries/PAPRECA/
make -f MakeFile clean
make -f MakeFile

cd ../../../Installation/Traditional\ Make/
mkdir ../../TraditionalMakeBuild

make -f MakeFile clean
make -f MakeFile
```

To remove all associated libraries, dependency files, and executables of a Traditional Make build you can use the "traditionalMakeClean.sh" (included in the Installation/Traditional Make/ directory).

```bash
bash traditionalMakeClean.sh
```

This is what the script "traditionalMakeClean.sh" does in detail:

```bash
#Uninstall PAPRECA library
cd ../../source/libraries/PAPRECA/
make -f MakeFile clean
rm libKMC.a

cd ../../../Installation/Traditional\ Make/
rm -R ../../TraditionalMakeBuild
```
