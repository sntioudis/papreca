#!/bin/bash

#The random see and papreca paths are provided to the test_all_events.sh script

random_seed=$1 #Retrieve seed, path, and python version from script invocation
papreca_path=$2
python_version=$3

echo -e "Running ALL predefined event tests with a random seed of "$random_seed" \n \n"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running FORMATIONS test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ./formation\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
"$python_version" test_formations.py "$papreca_path"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running BREAKING test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ../breaking\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
"$python_version" test_breakings.py "$papreca_path"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running DEPOSITION test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ../deposition\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
"$python_version" test_depositions.py "$papreca_path"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running DIFFUSION test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ../diffusion\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
"$python_version" test_diffusions.py "$papreca_path"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n"
echo -e "Running MONOATOMIC DESORPTION test! \n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x \n \n \n"
cd ../monodesorption\ events/
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
"$python_version" test_monodesorptions.py "$papreca_path"
