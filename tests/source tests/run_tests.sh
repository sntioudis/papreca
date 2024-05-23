#!/bin/bash
paprecatests_dir=$1 #Retrieve path from script invocation

mpiexec "${paprecatests_dir}/source_tests" -in ./input_files/in_kmc.lmp ./input_files/in_kmc.ppc
