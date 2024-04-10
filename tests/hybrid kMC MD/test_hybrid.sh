#!/bin/bash
papreca_dir=$1 #Retrieve path from script invocation

#Run test
mpiexec "${papreca_dir}/papreca" -in in_kmc.lmp in_kmc.ppc
