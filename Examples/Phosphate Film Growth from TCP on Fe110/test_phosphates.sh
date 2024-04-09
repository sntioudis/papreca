#!/bin/bash
papreca_dir="../..//build"

#Run test
mpiexec -np 16 "${papreca_dir}/papreca" -in in_kmc.lmp in_kmc.ppc
