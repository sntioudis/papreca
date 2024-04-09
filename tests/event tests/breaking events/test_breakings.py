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

def extractDataFromLines( lines ):


    btype_arr = [] #2D arrays holding bonds from all steps
    atom1_arr = []
    atom2_arr = []
    
    btype_step = [] #1D array holding values of specific step
    atom1_step = []
    atom2_step = []
    
    append = False
    
    for i in range(len(lines)):
    
       	data = lines[i].split()
       	
        btype_step.append( int(data[1]) )
        atom1_step.append( int(data[2]) )
        atom2_step.append( int(data[3]) )
        

        if( i < len(lines) - 1 ): #If this is not the final line, check the id of the next bond. If it is 1, then it means you have reached the end of the current step (because the next bond id is 1 and corresponds to a new step
            data_next = lines[i+1].split()
            bondid_next = int( data_next[0] )
            if( bondid_next == 1 ):
                append = True
       	elif( i == len(lines)-1 ): #If this is the last line that you read, append it
       	    append = True
       	    
       	if( append ):
       	    btype_arr.append( btype_step )
       	    atom1_arr.append( atom1_step )
       	    atom2_arr.append( atom2_step )
       	    
       	    #And (re)define the step arrays to clear the relevant data
       	    btype_step = []
       	    atom1_step = []
       	    atom2_step = []
       	    
       	    #Finally, reset the append flag
       	    append = False

    return btype_arr , atom1_arr , atom2_arr
    

def bondExistsInBondsList( btype_LAMMPS , btype_PAPRECA , atom1_LAMMPS , atom1_PAPRECA , atom2_LAMMPS , atom2_PAPRECA ):
    
    if( btype_LAMMPS == btype_PAPRECA and atom1_LAMMPS == atom1_PAPRECA and atom2_LAMMPS == atom2_PAPRECA ):
        return True
    
    return False

def compareArraysAndPrintStats( btype_PAPRECA , btype_LAMMPS , atom1_PAPRECA , atom1_LAMMPS , atom2_PAPRECA , atom2_LAMMPS):
    
    print( "STEP SUMMARY" )
    print( "----------------------------------------------------------------" )
    print( "btype  atom1_id  atom2_id \n" )
    for i in range( len(btype_LAMMPS) ):
        print( "THIS IS STEP " , i )
        for j in range( len(btype_LAMMPS[i]) ):
            print( btype_LAMMPS[i][j] , " " , atom1_LAMMPS[i][j] , " " , atom2_LAMMPS[i][j] )
            
    print( "END OF STEP SUMMARY" )
    print( "----------------------------------------------------------------" )
    print( " " )
    
    
    step_success = 0
        
    #Scan LAMMPS lists. Skip the first one as it contains the starting bonds list (with no broken bonds)
    for i in range( 1 , len(btype_LAMMPS) ):
    
        print( "In step " , i-1 , " PAPRECA broke bondtype=" , btype_PAPRECA[i-1] , " between atoms " , atom1_PAPRECA[i-1] , " and " , atom2_PAPRECA[i-1] )
        print( "Checking if this bond exists in the LAMMPS bondlist of exactly the next step..." )
        for j in range( len(btype_LAMMPS[i]) ): #Loop if the bond that was recently broken exists in the LAMMPS bondlist
            
            if( bondExistsInBondsList( btype_LAMMPS[i][j] , btype_PAPRECA[i-1] , atom1_LAMMPS[i][j] , atom1_PAPRECA[i-1] , atom2_LAMMPS[i][j] , atom1_PAPRECA[i-1] ) ):
                continue
        step_success += 1
        print( "Success! The bond is not in the LAMMPS bondlist of the relevant step!" )#This line is reached only if you do not continue
        print( " " );
    
    
    #len(btype_LAMMPS) has to be larger than the number of breaks (i.e., len(btype_PAPRECA) ). This happens because the bonds_full.log contains
    #information about timestep 0 (before any event has been executed). Hence, if len(btype_LAMMPS) == len(btype_PAPRECA), then the last timestep has NO bonds
    #and therefore the last bond breaking event was successful (all bonds were broken). See the bonds_full.log file for more information
    if( len(btype_LAMMPS) == len(btype_PAPRECA) ):
        step_success += 1
        print( "In step " , len(btype_PAPRECA) -1 , " PAPRECA broke bondtype=" , btype_PAPRECA[len(btype_PAPRECA)-1] , " between atoms " , atom1_PAPRECA[len(btype_PAPRECA)-1] , " and " , atom2_PAPRECA[len(btype_PAPRECA)-1] )
        print( "Success! The final timestep contains no bonds, so this bond: (btype,atom1,atom2)= " , btype_PAPRECA[-1] , " " , atom1_PAPRECA[-1] , " " , atom2_PAPRECA[-1] , " is not present in the final timestep! " )
        
    print( " " )
    print( "PRINTING TEST SUMMARY" )
    print( "----------------------------------------------------------------" )
    print( "Stats after testing " , len(btype_PAPRECA) , " sequential bond-breaking events..." )
    print( 100 * float( step_success ) / len( btype_PAPRECA ) , "% success rate!" )
    print( "----------------------------------------------------------------" )
    
        

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
        print( "Papreca finished successfully...initiating breakings test!")
        print( " " )


    #Open the papreca_full.log file and get the bond formations events in the order they have been executed
    btype_PAPRECA = []
    atom1_PAPRECA = []
    atom2_PAPRECA = []
    
    # Open the file for reading
    with open('papreca_full.log', 'r') as file:
        # Iterate through each line in the file
        for line in file:
            # Check if the line starts with "Executing"
            if line.startswith(' Executing '):
                # Extract the numbers using regular expressions
                bond_type = re.search(r'bond_type=(\d+)', line)
                atom1_id = re.search(r'atom1_id = (\d+)', line)
                atom2_id = re.search(r'atom2_id = (\d+)', line)
                
                if bond_type and atom1_id and atom2_id:
                    # Append the numbers to their respective lists
                    btype_PAPRECA.append(int(bond_type.group(1)))
                    atom1_PAPRECA.append(int(atom1_id.group(1)))
                    atom2_PAPRECA.append(int(atom2_id.group(1)))

    btype_PAPRECA = np.array( btype_PAPRECA )
    atom1_PAPRECA = np.array( atom1_PAPRECA )
    atom2_PAPRECA = np.array( atom2_PAPRECA )
    
    #Now read the bonds_full.log file and see if the correct bonds have been broken on each step
    file_path = 'bonds_full.log'
    start_marker = 'ITEM: ENTRIES index c_bondsinfo[1] c_bondsinfo[2] c_bondsinfo[3]'
    end_marker = 'ITEM: TIMESTEP'
    
    #Collect lines and unpack data for btype, atom1, and atom2
    lines = readBetweenLines(file_path, start_marker, end_marker)
    
    #For this test, the btype_LAMMPS, atom1_LAMMPS, and atom2_LAMMPS are 2D, the first dimension represents the step and the second the bonds list
    btype_LAMMPS , atom1_LAMMPS , atom2_LAMMPS = extractDataFromLines( lines )
    
    #Now compare results and print stats
    compareArraysAndPrintStats( btype_PAPRECA , btype_LAMMPS , atom1_PAPRECA , atom1_LAMMPS , atom2_PAPRECA , atom2_LAMMPS )

if __name__ == "__main__":
    main()
