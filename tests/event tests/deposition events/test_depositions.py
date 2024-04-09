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

    mol = []
    x = []
    y = []
    z = []
    
    for line in lines:
        data = line.split()
        
        mol.append(int(data[0]))
        x.append(float(data[1]))
        y.append(float(data[2]))
        z.append(float(data[3]))
        
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    
    return mol, x, y, z 

def getMolIndexesList( mol_array , mol_id ):
    
    return [i for i, value in enumerate(mol_array) if value == mol_id]
    
def getMolCoordsFromMolIndexesList( mol_indexes , x , y , z ):
    x_mol = [x[i] for i in mol_indexes]
    y_mol = [y[i] for i in mol_indexes]    
    z_mol = [z[i] for i in mol_indexes]
    
    return np.array( x_mol ), np.array( y_mol ), np.array( z_mol )

def getCOGCoordsOfPoints( x , y , z ):
    
    return np.mean( x ) , np.mean( y ), np.mean( z )

def compareArraysAndPrintStats( X_PAPRECA , X_LAMMPS , Y_PAPRECA , Y_LAMMPS , Z_PAPRECA , Z_LAMMPS ):
    x_success = 0
    y_success = 0
    z_success = 0
    
    for i in range(len(X_PAPRECA)):
        if( X_PAPRECA[i] == X_LAMMPS[i] ):
           x_success += 1
        if( Y_PAPRECA[i] == Y_LAMMPS[i] ):
            y_success += 1
        if( Z_PAPRECA[i] == Z_LAMMPS[i] ):
            z_success += 1
        
        print( "PAPRECA deposit COG for test " , i , ": " , X_PAPRECA[i] , " " , Y_PAPRECA[i] , Z_PAPRECA[i] )
        print( "LAMMPS deposit COG for test " , i , ": " , X_LAMMPS[i] , " " , Y_LAMMPS[i] , Z_LAMMPS[i] )

    x_success = 100 * float( x_success ) / len(X_PAPRECA)
    y_success = 100 * float( y_success ) / len(X_PAPRECA)
    z_success = 100 * float( z_success ) / len(X_PAPRECA)
    
    print( " " )
    print( "PRINTING TEST SUMMARY" )
    print( "----------------------------------------------------------------" )
    print( "Deposition test results after comparing the COGs of " , len(X_PAPRECA) , " deposition events..." )
    print( "Successful X COGs: " , x_success )
    print( "Successful Y COGs: " , y_success )
    print( "Successful Z COGs: " , z_success )
    
    if( x_success == 100.0 and y_success == 100.0 and z_success == 100.0 ):
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
        print( "Papreca finished successfully...initiating depositions test!")
        print( " " )


    #Collect output of LAMMPS. It should be two snapshots (the first one and the last one) including the unwrapped coordinates of atoms along with their mol ids
    file_path = 'lammps_full.dat'
    start_marker = 'ITEM: ATOMS mol xu yu zu'
    end_marker = 'ITEM: TIMESTEP'

    #Collect lines and unpack data for mol, and for x, y, and z coordinates
    lines = readBetweenLines(file_path, start_marker, end_marker)
    mol, x, y, z = extractDataFromLines( lines )

    #For all deposition events (i.e., for mol ids going from 1 to Ndeposition events), extract the relevant mol indexes, get the coordinates of molecules, and find their center of gravity.
    #For this test we performed 5 deposition events, so we will be checking 5 molecules

    X_LAMMPS = np.empty(5) #CAPITAL letters are the arrays storing all the coordinates of all tested molecules.
    Y_LAMMPS = np.empty(5)
    Z_LAMMPS = np.empty(5)

    for i in range(len(X_LAMMPS)):
        #Get mol coords of specific mol id
        mol_indexes = getMolIndexesList( mol , i+1 )
        x_mol , y_mol , z_mol = getMolCoordsFromMolIndexesList( mol_indexes , x , y , z ) #small letters are single values storing (temporarily) the COG coordinates of the current molecule id

        #Get COG of specific mol ID
        x_cog , y_cog , z_cog = getCOGCoordsOfPoints( x_mol , y_mol , z_mol )
        
        #Append to coordinates array
        X_LAMMPS[i] = x_cog
        Y_LAMMPS[i] = y_cog
        Z_LAMMPS[i] = z_cog


    #Now open and read the papreca output file and see if the deposition coordinates as calculated by PAPRECA match the coordinates as reported by LAMMPS.
    #This means that the inserted molecule is placed exactly where intended, since the deposition center (as calculated by PAPRECA in the kMC stage) matches the center of gravity of the deposited molecule (as in the LAMMPS file)

    # Initialize empty lists to store the final three numbers
    X_PAPRECA = []
    Y_PAPRECA = []
    Z_PAPRECA = []

    # Open the file for reading
    with open('papreca_full.log', 'r') as file:
        # Iterate through each line in the file
        for line in file:
            # Check if the line starts with "Executing"
            if line.startswith(' Executing'):
                # Extract the final three numbers from the line
                numbers = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                if len(numbers) >= 3:
                    # Append the numbers to their respective lists
                    X_PAPRECA.append(float(numbers[-3]))
                    Y_PAPRECA.append(float(numbers[-2]))
                    Z_PAPRECA.append(float(numbers[-1]))


    X_PAPRECA = np.array( X_PAPRECA ) 
    Y_PAPRECA = np.array( Y_PAPRECA )
    Z_PAPRECA = np.array( Z_PAPRECA )


    #Round those coming from the PAPRECA run because the PAPRECA and LAMMPS report with different accuracy
    X_PAPRECA = np.round( X_PAPRECA , 4 )
    X_LAMMPS = np.round( X_LAMMPS , 4 )

    Y_PAPRECA = np.round( Y_PAPRECA , 4 )
    Y_LAMMPS = np.round( Y_LAMMPS , 4 )

    Z_PAPRECA = np.round( Z_PAPRECA , 4 )
    Z_LAMMPS = np.round( Z_LAMMPS , 4 )


    #Compare COG results
    compareArraysAndPrintStats( X_PAPRECA , X_LAMMPS , Y_PAPRECA , Y_LAMMPS , Z_PAPRECA , Z_LAMMPS )


if __name__ == "__main__":
    main()
