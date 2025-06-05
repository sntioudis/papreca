/*----------------------------------------------------------------------------------------
PAPRECA hybrid off-lattice kinetic Monte Carlo/Molecular dynamics simulator.
Copyright (C) 2024 Stavros Ntioudis, James P. Ewen, Daniele Dini, and C. Heath Turner

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
----------------------------------------------------------------------------------------*/

/// \file ///
///@brief Independent c++ file testing core functionality of the PAPRECA software

//System Headers

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>
#include <cstdio>
#include <chrono>
#include <ctime>
#include <cmath>
#include <limits>
#include <mpi.h>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <algorithm>

//LAMMPS Headers
#include "lammps.h"
/// \cond
#include "input.h"
#include "atom.h"
#include "pair.h"
#include "thermo.h"
#include "output.h"
#include "library.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "compute.h"
#include "pointers.h"
#include "group.h"
#include "random_mars.h"
#include "molecule.h"
#include "math_extra.h"
/// \endcond

//PAPRECA kMC Headers
#include "papreca.h"

using namespace PAPRECA;
using namespace LAMMPS_NS;
using namespace std;

void initializeTests( int *narg , char ***arg , int *nprocs , int *proc_id , LAMMPS **lmp , PaprecaConfig &papreca_config ){
	
	/// Intializes MPI, LAMMPS, and %PAPRECA for source tests
	/// @param[in] narg number of command-line arguments passed to the main function (i.e., the papreca executable) during the program invocation from the terminal.
	/// @param[in] arg array containing the char* passed to the main function during the program invocation from the terminal.
	/// @param[in,out] nprocs number of MPI processes.
	/// @param[in,out] proc_id ID of current MPI process.
	/// @param[in,out] lmp pointer to LAMMPS instance.
	/// @param[in,out] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
	
	string lmp_input = (*arg)[ 2 ];
	const char *papreca_input = (*arg)[ 3 ];
	
	//MPI setup
	setupMPI( narg , arg , nprocs , proc_id );
	
	//LMP setup
	initializeLMP( lmp );
	readLMPinput( lmp_input , *lmp );
	
	//PAPRECA setup
	readInputAndInitPaprecaConfig( *lmp , *proc_id , papreca_input , papreca_config );
	papreca_config.setupExportFiles( *proc_id );
}

void finalizeTests( LAMMPS **lmp ){
	
	/// Finalizes MPI and lmp for source tests
	/// @param[in,out] lmp pointer of LAMMPS object
	
	delete *lmp; //Always delete this object last, otherwise you get a segmentation fault.
	MPI_Finalize();
	
}

void resetLAMMPS( LAMMPS **lmp , char ***arg , const int &proc_id ){
	
	/// Deletes the previously instantiated LAMMPS object and creates a new one to reset the system and perform another source code test. 
	/// @param[in,out] lmp pointer
	/// @param[in] arg array containing the char* passed to the main function during the program invocation from the terminal.
	/// @param[in] proc_id ID of current MPI process.
	/// @note This operation does not affect the global variables of %PAPRECA stored in the PAPRECA::PaprecaConfig object.
	
	//Delete lammps object
	delete *lmp;
	
	
	//Instantiate a new LAMMPS object and pointer lmp pointer to that instance.
	string lmp_input = (*arg)[ 2 ];
	
	if( proc_id == 0 ){
		printf( "\n \nPAPRECA MESSAGE: RESETTING LAMMPS OBJECT TO INITIAL SYSTEM STATE (i.e., THE ONE DEFINED IN THE LAMMPS INPUT FILE)... \n \n \n" );
	}
	
	initializeLMP( lmp );
	readLMPinput( lmp_input , *lmp ); //Here, we use the same LAMMPS input file
	
	
}

void testMolCoords( LAMMPS *lmp , PaprecaConfig &papreca_config , const int &proc_id ){
	
	///This check ensures that the molecule coordinates calculated by PAPRECA (and used for collision tests) are identical to the inserted molecule coordinates in the system (after executing the deposition event).
	/// @param[in,out] lmp pointer to LAMMPS instance.
	/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
	/// @param[in] proc_id ID of current MPI process.
	
	
	// First retrieve a molecule template (e.g., template set up for atoms of type 1).
	PredefinedDeposition *depo_template = papreca_config.getDepositionFromParentAtomType( 1 );

	
	// Get mol information from mol name
	double **mol_dx = depo_template->getCoords( );
	int *mol_atomtype = depo_template->getAtomTypes( );
	const int mol_natoms = depo_template->getAtomsNum( );
		
	//Initialize mol_xyz array storing molecule coordinates
	double **mol_xyz = NULL;
	initMolCoordsArr( &mol_xyz , mol_natoms );
	
	
	// PAPRECA predicts the molecule coordinates from the geometric center of the molecule (via PAPRECA::getMolCoords() functions).
	// Let's say the geometric center of the molecule is at (x,y,z) = (0.0,0.0,15.0).
	double candidate_center[3] = { lmp->domain->boxlo[0] , lmp->domain->boxlo[1] , 15.0 }; //We test an edge case here and place the atom on the corner of the periodic box to check that periodic boundary conditions are treated properly by the PAPRECA::getMolCoords() function.	
	getMolCoords( lmp , mol_xyz , mol_dx , mol_natoms , candidate_center );
	
	//Insert the molecule with the given center (and with no rotation).
	double rot_pos[3] = {1.0,0.0,0.0}; //Rotation vector does not really affect deposition (because the rotation angle is zero in PAPRECA::insertMolecule())
	insertMolecule( lmp , candidate_center , rot_pos , 0.0 , 0 , depo_template->getAdsorbateName( ).c_str( ) );
	resetMobileAtomsGroups( lmp , papreca_config );
	runLammps( lmp , 0 ); //Run 0 to update neighbor lists
	
	
	//Scan all atoms on all procs and if they are not type 1, you know they correspond to recently inserted atoms (via PAPRECA::insertMolecule()).
	//No need to gather atoms here, since we can perform scans on individual atoms owned by MPI processes and success counts later
	tagint *id = ( tagint *)lammps_extract_atom( lmp , "id" );
	int *type = (int *) lammps_extract_atom( lmp , "type" );
	int natoms = *(int *) lammps_extract_global( lmp , "nlocal" );
	double **pos = ( double **)lammps_extract_atom( lmp , "x" );


	int success_local = 0 , success_global = 0 , tests_local = 0 , tests_global = 0;
	double epsilon = 1.0e-16;
	for( int i = 0; i < natoms; ++i ){
		
		if( type[i] != 1 ){
			++tests_local;
			for( int j = 0; j < mol_natoms; ++j ){
				if( fabs( pos[i][0] - mol_xyz[j][0] ) < epsilon ){ //Compare (naively) system coords with molecular coordinates (as predicted by PAPRECA, via getMolCoords function)
					if( fabs( pos[i][1] - mol_xyz[j][1] ) < epsilon ){
						if( fabs( pos[i][2] - mol_xyz[j][2] ) < epsilon ){
							++success_local;
						}
					}
				}
			}
		}
		
			
	}
	
	MPI_Reduce( &success_local , &success_global , 1 , MPI_INT , MPI_SUM , 0 , MPI_COMM_WORLD ); //Reduce all success values on the master process (i.e., proc_id == 0 0) to determine successful comparisons.
	MPI_Reduce( &tests_local , &tests_global , 1 , MPI_INT , MPI_SUM , 0 , MPI_COMM_WORLD );

	double success_rate = 0.0;
	if( proc_id == 0 ){
		
		printf( "\n \nPRINTING MOLECULE COORDINATES TEST SUMMARY \n" );
		printf( "---------------------------------------------------------------- \n" );
		printf( "System atoms: %ld \n" , lmp->atom->natoms );
		printf( "Molecule atoms: %d \n" , mol_natoms );
		printf( "Total coordinate comparisons: %d \n" , tests_global );
		
		success_rate = 100.0 * static_cast<double>( success_global ) / tests_global;
		printf( "SUCCESS RATE: %f %% \n" , success_rate );
		printf( "----------------------------------------------------------------\n \n \n \n" );
		
	}

	MPI_Bcast( &success_rate , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); //Cast calculated success rate to all other processes
	if( success_rate < 100.0 ){ allAbortWithMessage( MPI_COMM_WORLD , "testMolCoords function in source_tests.cpp failed!" ); }

	

}


void testCollisions( LAMMPS *lmp , PaprecaConfig &papreca_config , const int &proc_id ){
	
	///Checks collisions between interfering atoms are predicted accurate by PAPRECA (during event detection).
	/// @param[in,out] lmp pointer to LAMMPS instance.
	/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
	/// @param[in] proc_id ID of current MPI process.
	
	//From papreca config retrieve the sigma values between 2 atom types (say 2 and 3 ).
	double sigma = papreca_config.getSigmaFromAtomTypes( 2 , 3 ); //Note: for types 2 and 3 in the current system, the sigma value should be 3.47299 Angstroms (as defined in the LAMMPS input file).
	if( sigma != 3.47299 ){
		if( proc_id == 0 ){
			printf( "TEST WARNING: SIGMA VALUE EXPECTED: %f BUT OBTAINED %f \n" , 3.47299 , sigma );
		}
	}
	
	//Create atoms of different atom types
	double atom1_pos[3] = { lmp->domain->boxlo[0] , lmp->domain->boxlo[1] , 20.0 }; //Create atom on a corner to check for edge case (i.e., see if periodic boundary conditions are treated properly).
	remap3DArrayInPeriodicBox( lmp , atom1_pos ); 
	createAtom( lmp , atom1_pos , 2 );
	
	double atom2_pos[3] = { lmp->domain->boxlo[0] , lmp->domain->boxlo[1] , 20.0 + 3.47 }; //Create another atom of type 3
	remap3DArrayInPeriodicBox( lmp , atom2_pos ); 
	createAtom( lmp , atom2_pos , 3 );
	
	double atom3_pos[3] = { lmp->domain->boxhi[0] -10, lmp->domain->boxhi[1] , 30.0 }; //use %PAPRECA wrapper remap3DArrayInPeriodicBox to remap the atom in the periodic box.
	remap3DArrayInPeriodicBox( lmp , atom3_pos ); 
	createAtom( lmp , atom3_pos , 2 );
	
	double atom4_pos[3] = { lmp->domain->boxhi[0] , lmp->domain->boxhi[1] -0.1 , 50.0 };
	remap3DArrayInPeriodicBox( lmp , atom4_pos ); 
	createAtom( lmp , atom4_pos , 3 );
	
	
	resetMobileAtomsGroups( lmp , papreca_config );
	runLammps( lmp , 0 ); //Run 0 to update neighbor lists
	

	//We created 2 atoms, but only 2 atoms interfere (i.e., atom1 and atom2 because their interatomic distance is smaller than sigma).
	//Now, we will search the neighbors lists of atoms in the systems and use the atomsCollide function.
	//ONLY 1 COLLISION SHOULD BE DETECTED IN A SUCCESSFUL TEST.
	
	//We will do the search through the HALF neighbors list, so retrieve list information
	int neiblist_id = lammps_find_fix_neighlist( lmp , "papreca" , 2 ); //Get neighbors list with ID 2 (half list as in the papreca fix)
	if( neiblist_id == -1 ){ allAbortWithMessage( MPI_COMM_WORLD , "Lammps could not find neib list with name " + papreca_config.getFullNeibListName( ) + ". Either the list does not exist or there is a spelling error in your PAPRECA input file." ); }
	int atoms_num = lammps_neighlist_num_elements( lmp , neiblist_id );
	
	//Retrieve local (per MPI proccess) system information
	tagint *id = ( tagint *)lammps_extract_atom( lmp , "id" );
	int *type = (int *) lammps_extract_atom( lmp , "type" );
	double **pos = ( double **)lammps_extract_atom( lmp , "x" );
		
	int iatom = -1, neighbors_num = - 2, *neighbors = NULL;
	//Loop over  full list
	int collisions_local = 0 , collisions_global = 0;
	for ( int i = 0; i < atoms_num; ++i ){
			
		lammps_neighlist_element_neighbors( lmp , neiblist_id , i , &iatom , &neighbors_num , &neighbors ); //get local atom index (iatom), number of neighbors of iatom, and indexes of iatom neighbors
		
		if( type[iatom] != 1 ){ //Only check for collisions types other than 1 (types 1 are Fe atoms on a lattice and their distance can be smaller than sigma)
			for( int j = 0; j < neighbors_num; ++j ){
				
				int jneib = getMaskedNeibIndex( neighbors , j );
				
				if( type[jneib] != 1 ){
					
					
					if( atomsCollide( lmp , papreca_config , pos[iatom] , type[iatom] , pos[jneib] , type[jneib] ) ){ //Only 2 atoms collide and we are testing on a half neighbor list (i.e., each pair of atoms is included once). Hence, ONLY ONE COLLISION SHOULD BE DETECTED FOR A SUCCESSFUL TEST.
						++collisions_local;
					}
				}
		
			}
		}
	}
	
	MPI_Allreduce( &collisions_local , &collisions_global , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD ); //Get collisions from all procs on master proc
	
	if( proc_id == 0 ){
		
		printf( "\n \nPRINTING COLLISIONS TEST SUMMARY \n" );
		printf( "---------------------------------------------------------------- \n" );
		printf( "Total atom insertions: 4 \n" );
		printf( "Total detected collisions: %d \n" , collisions_global );
		
		if( collisions_global == 1 ){
			printf( "The test was SUCCESSFUL \n" );
		}else{
			printf( "The test was UNSUCCESSFUL \n" );
		}
		printf( "----------------------------------------------------------------\n \n \n \n" );
		
	}

	if( collisions_global != 1 ){ allAbortWithMessage( MPI_COMM_WORLD , "testCollisions function in source_tests.cpp failed!" ); }
	
}

void testRandomNumberGenerator( PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs ){
	
	/// Tests if the random number generator can produce a unique sequence of numbers. The user has to select the test limit.
	/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
	/// @param[in] proc_id ID of current MPI process.
	/// @param[in] nprocs total number of MPI processes.
	/// @note Random numbers are produced on the master MPI proc. There is no need to draw numbers on all procs (since the same random number is used and the same collisions will be detected).
	/// @note The test checks if the RanMars class of LAMMPS (check RanMars.h header of LAMMPS) can produce sequences of numbers with few repetitions.
	
	
	if( proc_id != 0 ){ return; } //Only run this test for the master proc
	
	int test_limit = 1.0e6;
	std::unordered_set<double> generated;

	int repetitions = 0;
	for( long int i = 0; i < test_limit; ++i ){

		const double rnum = papreca_config.getUniformRanNum( );
			
		if( !elementIsInUnorderedSet( generated , rnum ) ){
			generated.insert( rnum );
		}else{
			++repetitions;
		}
		
	}
		
	printf( "\n \n RANDOM NUMBERS TEST SUMMARY \n" );
	printf( "---------------------------------------------------------------- \n" );
	printf( "A total of %d tests were performed...\n" , test_limit );
	printf( "A total of %d random number repetitions were detected (%f %%) \n" , repetitions , 100 * static_cast<double>(repetitions) / ( test_limit ) );
	printf( "----------------------------------------------------------------\n \n \n \n" );


}


int main( int narg , char **arg ){

	/// Driver function for source tests.
	/// @param[in] narg number of command-line arguments passed to the source_tests executable function during the program invocation from the terminal.
	/// @param[in] arg array containing the char* passed during the program invocation from the terminal.
	

	//Declare variables
	PaprecaConfig papreca_config;
	int proc_id , nprocs;
	LAMMPS *lmp = NULL;
	
	initializeTests( &narg , &arg , &nprocs , &proc_id , &lmp , papreca_config );
	

	//Mol coords test
	testMolCoords( lmp , papreca_config , proc_id );
	resetLAMMPS( &lmp , &arg , proc_id );
	
	//Collisions test
	testCollisions( lmp , papreca_config , proc_id );
	resetLAMMPS( &lmp , &arg , proc_id );
	
	//Random Numbers Test
	testRandomNumberGenerator( papreca_config , proc_id , nprocs );
	
	
	finalizeTests( &lmp );
	
	return 0;
	
}
