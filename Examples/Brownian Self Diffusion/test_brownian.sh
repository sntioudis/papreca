#!/bin/bash
papreca_dir="../../build"

#Run test
mpiexec "${papreca_dir}/papreca" -in in_kmc.lmp in_kmc.ppc

#Move results to results folder
mkdir results
mv ./*.log ./results

#Run dedicated python script to plot distributions (change the input to the python script if you want to plot other distributions from the distributions folder).
#A distributions.jpg file will be generated in the parent directory with the requested plot
python3 plot_distributions.py ./results 10 4000 10000 20000 100000

#Remember to delete the results folder if you wish to rerun
