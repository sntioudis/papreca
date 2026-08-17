#!/bin/bash

# The random seed, PAPRECA path, and Python executable are provided by the
# workflow/test caller.
#
# Diagnostic behavior:
#   * Every predefined event test is attempted, even if an earlier test fails.
#   * Each event directory receives:
#         test_runner.log
#         test_exit_code.txt
#   * The script exits with status 1 at the end if one or more tests failed.

random_seed=$1
papreca_path=$2
python_version=$3

failed_tests=()

run_test() {
    if [ "$#" -ne 4 ]; then
        echo "Error: run_test requires exactly 4 arguments:"
        echo "       test_name, python_version, python_script, papreca_path"
        return 1
    fi

    local test_name=$1
    local python_version=$2
    local python_script=$3
    local papreca_path=$4

    echo
    echo "Running ${test_name} Python test..."
    echo "Python script: ${python_script}"
    echo "Working directory: $(pwd)"
    echo

    "$python_version" "$python_script" "$papreca_path" 2>&1 | tee test_runner.log
    local exit_code=${PIPESTATUS[0]}

    echo "$exit_code" > test_exit_code.txt

    if [ "$exit_code" -ne 0 ]; then
        echo
        echo "${test_name} test FAILED with exit code ${exit_code}."
        echo "Continuing so that the remaining event tests also produce logs."
        failed_tests+=("$test_name")
        return 1
    fi

    echo
    echo "${test_name} test PASSED."
    return 0
}

echo -e "Running ALL predefined event tests with a random seed of ${random_seed}\n\n"

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n"
echo -e "Running FORMATIONS test!\n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n\n\n"

cd ./formation\ events/ || exit 1
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "FORMATIONS" "$python_version" test_formations.py "$papreca_path" || true

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n"
echo -e "Running BREAKING test!\n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n\n\n"

cd ../breaking\ events/ || exit 1
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "BREAKING" "$python_version" test_breakings.py "$papreca_path" || true

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n"
echo -e "Running DEPOSITION test!\n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n\n\n"

cd ../deposition\ events/ || exit 1
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "DEPOSITION" "$python_version" test_depositions.py "$papreca_path" || true

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n"
echo -e "Running DIFFUSION test!\n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n\n\n"

cd ../diffusion\ events/ || exit 1
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "DIFFUSION" "$python_version" test_diffusions.py "$papreca_path" || true

echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n"
echo -e "Running MONOATOMIC DESORPTION test!\n"
echo -e "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x\n\n\n"

cd ../monodesorption\ events/ || exit 1
cp in_kmc_template.ppc in_kmc.ppc
sed -i "s/SEEDHOLDER/${random_seed}/" in_kmc.ppc
run_test "MONOATOMIC DESORPTION" "$python_version" test_monodesorptions.py "$papreca_path" || true

echo
echo "================================================================"
echo "EVENT TEST SUMMARY"
echo "================================================================"

if [ "${#failed_tests[@]}" -gt 0 ]; then
    echo "The following event tests failed:"
    for test_name in "${failed_tests[@]}"; do
        echo "  - ${test_name}"
    done
    echo
    echo "All remaining event tests were still executed for diagnostics."
    exit 1
fi

echo "ALL TESTS WERE SUCCESSFUL!"
exit 0
