#!/bin/bash

#The random see and papreca paths are provided to the test_all_events.sh script

random_seed=$1 #Retrieve seed, path, and python version from script invocation
papreca_path=$2
python_version=$3


run_test( ){

	# Check if exactly three arguments are provided
	if [ "$#" -ne 3 ]; then
		echo "Error: run_test function requires exactly 3 arguments (python_version, python_script name, and papreca_path)."
		return 1
	fi
	
	#Retrieve paths and python script name
	local python_version=$1
	local python_script=$2
	local papreca_path=$3
	
	#Run test
	"$python_version" "$python_script" "$papreca_path"
	
	#Retrieve exit code and abort if necessary
	exit_code=$?
	if [ $exit_code -eq 1 ]; then
        	echo "Test failed...aborting!"
        	exit 1
    	fi
	
}
echo -e "Running ALL predefined event tests with a random seed of "$random_seed" \n \n"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running FORMATIONS test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ./formation\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "$python_version" test_formations.py "$papreca_path"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running BREAKING test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ../breaking\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "$python_version" test_breakings.py "$papreca_path"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running DEPOSITION test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ../deposition\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "$python_version" test_depositions.py "$papreca_path"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running DIFFUSION test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ../diffusion\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "$python_version" test_diffusions.py "$papreca_path"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running MONOATOMIC DESORPTION test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ../monodesorption\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "$python_version" test_monodesorptions.py "$papreca_path"

echo -e "ALL TESTS WERE SUCCESSFUL!" #If you reach this point it means that all tests exited with 100% success rate
