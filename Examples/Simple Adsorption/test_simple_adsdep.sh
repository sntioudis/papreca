#!/bin/bash

papreca_dir="../../build"


#resolve the real path of papreca_dir (this can prevent conflicts if the user provides a relative path).
papreca_dir=$(realpath "$papreca_dir")

#List ra and rd values
ra_values=(1 3 5)
rd_values=(1 3 5)

for ra in "${ra_values[@]}"; do
    for rd in "${rd_values[@]}"; do
    
        # Create the folder if it doesn't exist
        folder="ra${ra}rd${rd}"
        mkdir "./$folder"
        
        # Copy template files in the created folder
        cp -R ./template/* "${folder}/"
        
        # Change directory to created folder
        cd "./${folder}"
        
        # Change ra and rd in template PAPRECA input files
        sed -i "s/ADSRATEPLACEHOLDER/${ra}/" in_kmc.ppc
        sed -i "s/DESRATEPLACEHOLDER/${rd}/" in_kmc.ppc
        
	# Run PAPRECA simulation
	mpiexec "${papreca_dir}/papreca" -in in_kmc.lmp in_kmc.ppc
	
	#Go back to the parent directory and call the Python script with the relevant arguments
	cd ..
        python3 compare_plots.py $ra $rd "${folder}/surface_coverage.log"
        
    done
done
