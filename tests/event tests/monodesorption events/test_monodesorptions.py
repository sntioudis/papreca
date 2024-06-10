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

def getDesorbedIdsFromLammps( lines ):
    #To get the desoprtion IDs we first add all the ids in a set. Those are the initial (original) atom IDs in the system before the execution of any monoatomic desorption events
    
    index_next = 0
    initial_ids = set()
    for i in range (len(lines)):
        data = lines[i].split()
        id_current = int(data[0])
        
        initial_ids.add(id_current)
        
        if( i < len(lines) - 1 ): #Ensure that you will not try to access an invalid memory block
            data_next = lines[i+1].split()
            id_next = int( data_next[0] )
            if( id_next < id_current ): #This means that we have reached the end of the current
                index_next = i+1
                break
                
    #Now loop through the remaining IDs, identify the deleted (desorbed) atom of each step and form the id_LAMMPS array
    #Note: A new set is created to avoid duplicate insertions to the id_LAMMPS array
    id_LAMMPS = [] #array holding the ids of the desorbed atoms (with the correct events order)
    step_ids = set() #temporary set holding the ids of the atoms on a given set
    id_LAMMPS_set = set() #auxiliary set to ensure that no duplicates are added in the id_LAMMPS array
    
    step_done = False
    for i in range (index_next,len(lines)):
        data = lines[i].split()
        id_current = int(data[0])
        
        step_ids.add( id_current )
        
        #Check if step is over
        if( i < len(lines) - 1 ):
            data_next = lines[i+1].split()
            id_next = int( data_next[0] )
            if( id_next < id_current ):
                step_done = True
        elif( i == len(lines) -1 ): #This is an edge case where you reach the end of the lines array.
            step_done = True
        
        #if the step is over compare the current step ids with the initial step ids
        if( step_done ):
            id_diff = initial_ids - step_ids
            
            #Now loop through the difference (i.e., the missing elements and determine if you have to add them in the id_LAMMPS array with the help of the id_LAMMPS_set container
            for id in id_diff:
                if id not in id_LAMMPS_set:
                    id_LAMMPS.append( id )
                    id_LAMMPS_set.add( id )
            
            #In any case clear the temporary set IDs set and reset the step_done flag
            step_ids.clear()
            step_done = False
        
        
    return np.array( id_LAMMPS )

def compareArraysAndPrintStats( id_PAPRECA , id_LAMMPS ):

    success = 0
    #Compare trimmed ids between LAMMPS and PAPRECA lists
    
    for i in range(len(id_PAPRECA)):
        if( id_PAPRECA[i] == id_LAMMPS[i] ):
            success += 1
        
        print( "PAPRECA desorbed atom ID for step " , i , ": " , id_PAPRECA[i] )
        print( "LAMMPS desorbed atom ID for step " , i , ": " , id_LAMMPS[i] )
        
    success = 100 * float( success ) / len(id_PAPRECA)
    
    print( " " )
    print( "PRINTING TEST SUMMARY" )
    print( "----------------------------------------------------------------" )
    print( "Monoatomic desorption tests after comparing the ids of " , len(id_PAPRECA) , " events..." )
    print( "Successful id comparisons: " , success)
    print( "The test success rate was: " + str( success ) )
    print( "----------------------------------------------------------------  \n \n \n" )
    
    return success

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
        print( "Papreca finished successfully...initiating monoatomic desorptions test!")
        print( " " )


    #Collect all output data from LAMMPS
    file_path = 'lammps_full.dat'
    start_marker = 'ITEM: ATOMS id'
    end_marker = 'ITEM: TIMESTEP'

    #Collect lines
    lines = readBetweenLines(file_path, start_marker, end_marker)
    id_LAMMPS = getDesorbedIdsFromLammps( lines )

    # Now get the desorbed IDs from the PAPRECA file
    id_PAPRECA = []
    with open('papreca_full.log', 'r') as file:
        # Iterate through each line in the file
        for line in file:
            # Check if the line starts with "Executing"
            if line.startswith(' Executing '):
                # Extract the numbers using regular expressions
                atom_id = re.search(r'ATOM_ID=(\d+)', line)
                
                if atom_id:
                    # Append the numbers to the list
                    id_PAPRECA.append(int(atom_id.group(1)))


    id_PAPRECA = np.array( id_PAPRECA )


    #Compare ID results
    success = compareArraysAndPrintStats( id_PAPRECA , id_LAMMPS )
    
    #Exit with the relevant code
    if( success < 100.0 ):
        sys.exit(1) #1 means failed test, while 0 means successful test. Those codes are handled by the caller bash script to abort prematurely
    else:
        sys.exit(0)
        


if __name__ == "__main__":
    main()
