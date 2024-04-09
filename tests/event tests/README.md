Predefined Event Tests
========

Each folder contains tests to ensure that the predefined events in PAPRECA work as expected.
The test logic is the same for all tests. Firstly, a PAPRECA run is initiated from a Python script.
PAPRECA runs for a few kMC steps and executes a number of predefined events of the same class (i.e., diffusion, deposition, desorption, and reaction steps)
Then, the Python script collects the event information as exported by PAPRECA and compares it with the system configuration as exported by LAMMPS.
In general, this guarantees that PAPRECA (i.e., the kMC stage) communicates properly with LAMMPS and that the system state is updated accordingly.
For example, if PAPRECA selects to execute a bond formation event between atoms with ids 1 and 5, then, such bond should appear in the output bonds file as exported by LAMMPS.

Running event tests
========

You can run an individual event test (from within the relevant test folder) by executing:

```bash
python testname.py /path/to/your/papreca/executable/
```

or

```bash
python3 testname.py /path/to/your/papreca/executable/
```

Note that, to run individual tests you would have to modify the in_kmc.ppc files so that they contain a valid random_seed.


Alternatively, you can run all tests by executing the "test_all_events.sh" bash script:

```bash
bash test_all_events.sh random_seed /path/to/your/papreca/executable/ python_version
```

Where the random_seed is an integer number between 0 and 900000000 and python_version is the full name of your python package (e.g., python3).

Example usage:

```bash
bash test_all_events.sh 49123 /home/sn120/Desktop/Codes/mypapreca/build python3
```
