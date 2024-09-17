#!/bin/bash
papreca_dir="../../build/PAPRECA/"

#Run test
mpiexec "${papreca_dir}/papreca" -in in_kmc.lmp in_kmc.ppc

#Move results to results folder
mkdir distributions
find . -maxdepth 1 -name "*.log" ! -name "papreca.log" ! -name "execTimes.log" -print | xargs -I {} mv {} ./distributions/ #Move all log files expect for the papreca.log to results folder

#Run dedicated python script to plot distributions (change the input to the python script if you want to plot other distributions from the distributions folder).
#A distributions.jpg file will be generated in the parent directory with the requested plot
python3 plot_distributions.py ./distributions 1000 4000 10000 30000 100000 200000

#Remember to delete the results folder if you wish to rerun
