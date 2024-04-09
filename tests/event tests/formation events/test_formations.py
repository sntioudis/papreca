import numpy as np
import re
import sys
import os
import subprocess


def readBetweenLines(file_path, start_marker, end_marker):
    reading = False
    lines = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() == start_marker:
                reading = True
                continue #To skip the label line
            if line.strip() == end_marker:
                reading = False
            if reading:
                lines.append(line.strip())
    return np.array( lines )


def bondExistsInBondsArrays( atom1 , atom2 , atom1_arr , atom2_arr ):
    
    for i in range(len(atom1_arr)):
        if( atom1 == atom1_arr[i] and atom2 == atom2_arr[i] ):
            return True
    return False
    
def extractDataFromLines( lines ):

    btype_arr = []
    atom1_arr = []
    atom2_arr = []
    
    for line in lines:
        data = line.split()
        btype = int(data[1])
        atom1 = int(data[2])
        atom2 = int(data[3])

        
        if( not bondExistsInBondsArrays( atom1 , atom2 , atom1_arr , atom2_arr ) ):
            btype_arr.append(btype)
            atom1_arr.append(atom1)
            atom2_arr.append(atom2)
    return np.array( btype_arr ) , np.array( atom1_arr ) , np.array( atom2_arr )
    

def compareArraysAndPrintStats( btype_PAPRECA , btype_LAMMPS , atom1_PAPRECA , atom1_LAMMPS , atom2_PAPRECA , atom2_LAMMPS ):

    btype_success = 0
    atom1_success = 0
    atom2_success = 0
    
    for i in range(len(btype_PAPRECA)):
        if( btype_PAPRECA[i] == btype_LAMMPS[i] ):
            btype_success += 1
        if( atom1_PAPRECA[i] == atom1_LAMMPS[i] ):
            atom1_success += 1
        if( atom2_PAPRECA[i] == atom2_LAMMPS[i] ):
            atom2_success += 1
            
        print( "PAPRECA formed the following bond: " , "byte= " , btype_PAPRECA[i] , " atom1_id=" , atom1_PAPRECA[i] , " atom2_id=" , atom2_PAPRECA[i] , " in step " , i )
        print( "The following bond: " , "byte=" , btype_LAMMPS[i] , " atom1_id=" , atom1_LAMMPS[i] , " atom2_id=" , atom2_LAMMPS[i] , " exists in the LAMMPS bonds lists, in step " , i )

    btype_success = 100 * float( btype_success ) / len(btype_PAPRECA)
    atom1_success = 100 * float( atom1_success ) / len(btype_PAPRECA)
    atom2_success = 100 * float( atom2_success ) / len(btype_PAPRECA)
    
    print( " " )
    print( "PRINTING TEST SUMMARY" )
    print( "----------------------------------------------------------------" )
    print( "Bond formation test results after comparing the bond types, atom1_ids, and atom2_ids of " , len(btype_PAPRECA) , " PAPRECA events..." )
    print( "Successful btypes: " , btype_success )
    print( "Successful atom1_ids: " , atom1_success )
    print( "Successful atom2_ids: " , atom2_success )
    
    if( btype_success == 100.0 and atom1_success == 100.0 and atom2_success ):
        print( "The test was 100% successful!" )
    print( "----------------------------------------------------------------  \n \n \n" )    

def main():

    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python3(or python) test_depositions.py path/to/papreca/executable")
        return
        
    # Get the path from the command-line argument
    papreca_path = sys.argv[1]
    
    #Run PAPRECA with MPI and send screen output to a file called papreca_full.log
    print( "Running PAPRECA...")
    command = "mpiexec " + papreca_path + "/papreca -in in_kmc.lmp in_kmc.ppc > papreca_full.log"
    mpi_return = subprocess.run(command, shell=True )
    
    if mpi_return.returncode != 0:
        print( "Error: PAPRECA did not finish successfully! Please check your papreca executable path" )
        sys.exit(1)
    else:
        print( "Papreca finished successfully...initiating formations test!")
        print( " " )


    #Open the papreca_full.log file and get the bond formations events in the order they have been executed
    btype_PAPRECA = []
    atom1_PAPRECA = []
    atom2_PAPRECA = []
    
    # Open the file for reading
    with open('papreca_full.log', 'r') as file:
        # Iterate through each line in the file
        for line in file:
            # Check if the line starts with "Executing bond formation event"
            if line.startswith(' Executing bond formation event'):
                # Extract the numbers using regular expressions
                bond_type = re.search(r'BOND_TYPE=(\d+)', line)
                atom1_id = re.search(r'ATOM1_ID=(\d+)', line)
                atom2_id = re.search(r'ATOM2_ID=(\d+)', line)
                
                if bond_type and atom1_id and atom2_id:
                    # Append the numbers to their respective lists
                    btype_PAPRECA.append(int(bond_type.group(1)))
                    atom1_PAPRECA.append(int(atom1_id.group(1)))
                    atom2_PAPRECA.append(int(atom2_id.group(1)))

    btype_PAPRECA = np.array( btype_PAPRECA )
    atom1_PAPRECA = np.array( atom1_PAPRECA )
    atom2_PAPRECA = np.array( atom2_PAPRECA )
    
    #Now read the bonds_full.log file and see if exactly the same bonds are present and in the correct order (i.e., bond with id 1 as executed by PAPRECA has to have the correct corresponding atom IDs)
    file_path = 'bonds_full.log'
    start_marker = 'ITEM: ENTRIES index c_bondsinfo[1] c_bondsinfo[2] c_bondsinfo[3]'
    end_marker = 'ITEM: TIMESTEP'
    
    #Collect lines and unpack data for btype, atom1, and atom2
    lines = readBetweenLines(file_path, start_marker, end_marker)
    btype_LAMMPS , atom1_LAMMPS , atom2_LAMMPS = extractDataFromLines( lines )
    
    #Now compare results and print stats
    compareArraysAndPrintStats( btype_PAPRECA , btype_LAMMPS , atom1_PAPRECA , atom1_LAMMPS , atom2_PAPRECA , atom2_LAMMPS )

if __name__ == "__main__":
    main()
