Source code tests
========

This folder contains a c++ code named tests.cpp. Within the test.cpp code, tests can be defined and executed to examine if any PAPRECA function (i.e., part of the PAPRECA library) operates as expected.

Compiling source code tests
========
The code has to be compiled against the **papreca** executable as well as the relevant libraries (i.e., **LAMMPS** and **PAPRECA**).

The tests.cpp code can be compiled as follows:

Firstly, create a "build" folder within the ./source\ tests/ directory:

```bash
mkdir build
cd ./build/
```

Then, run cmake:

```bash
cmake .. -DLAMMPS_SRC_DIR=/path/to/your/LAMMPS/src/ -DLAMMPS_LIB_DIR=/path/to/your/LAMMPS/library/
```

Finally, build with cmake:

```bash
cmake --build .
```

Note that, the compilation procedure is similar to that of the papreca source code. Once again, the path of the LAMMPS source files and library (usually located in the LAMMPS src and build folders, respectively) have to be provided to cmake.
Additionally, the path to the PAPRECA source folder has to be provided (through the PAPRECA_SRC_DIR variable) to enable linking against the PAPRECA source files.

Running source code tests
========

The code can be executed directly from the ./source\ tests/ directory. Please use the relevant (PAPRECA and LAMMPS) input files provided in the ./input\ files/ directory. Such input files set up the simulation domain and initialize
the PaprecaConfig object.

```bash
mpirun ../build/source_tests -in ./input_files/in_kmc.lmp ./input_files/in_kmc.ppc
```

or

```bash
mpiexec ../build/source_tests -in ./input_files/in_kmc.lmp ./input_files/in_kmc.ppc
```
