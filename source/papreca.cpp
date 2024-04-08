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
///@brief Driver function. Initializes MPI, LAMMPS, and PAPRECA. Then runs kMC/MD loops.

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

//Forward Declarations
using namespace LAMMPS_NS;
using namespace std;
using namespace PAPRECA;

void initialize( int *narg , char ***arg , int *nprocs , int *proc_id , LAMMPS **lmp , PaprecaConfig &papreca_config ){
	
	/// Intializes MPI, LAMMPS, and PAPRECA
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


void finalize( LAMMPS **lmp , PaprecaConfig &papreca_config , const int &proc_id ){
	
	/// Called at the end of the PAPRECA simulation. This function 1) closes all open files (if any), 2) deletes the previously instantiated LAMMPS object (in initializeLMP() function), and 3) calls MPI_Finalize() to ensure proper termination of all MPI processes.
	/// @param[in,out] lmp pointer to pointer of LAMMPS object (used to delete the instantiated LAMMPS object in initializeLMP() function).
	/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
	/// @param[in] proc_id ID of current MPI process.
	/// @see PAPRECA::initializeLMP()
	/// @note No need to delete anything related to the PAPRECA::PaprecaConfig object (i.e., the object carrying global parameters and simulation settings). At the moment, all variables and containers initialized for the PAPRECA::PaprecaConfig object are allocated on stack and they will be deleted after the object destructor invokation (i.e., at the end of the main() function).
	
	papreca_config.closeExportFiles( proc_id );
	delete *lmp; //Always delete this object last, otherwise you get a segmentation fault.
	MPI_Finalize();
	
}


int main( int narg , char **arg ){
	
	/// Driver function running the main PAPRECA simulation loop. The function sets up the MPI protocol, initializes all (LAMMPS and PAPRECA) variables, and performs the requested (by the user, in the PAPRECA input file) PAPRECA simulation KMC steps. On each PAPRECA simulation KMC step, each atom on every MPI process is scanned and PAPRECA::Events are discovered. Then, an event is executed on an MPI processes (the executed event as well as the MPI process firing the event are chosen based on the N-FOLD way).
	/// @param[in] narg number of command-line arguments passed to the main function (i.e., the papreca executable) during the program invocation from the terminal.
	/// @param[in] arg array containing the char* passed to the main function during the program invocation from the terminal.
	/// @see PAPRECA::setupMPI(), PAPRECA::initializeLMP(), PAPRECA::readLMPinput(), PAPRECA::readInputAndInitPaprecaConfig(), PAPRECA::Bond::initAtomID2BondsMap(), PAPRECA::loopAtomsAndIdentifyEvents(), PAPRECA::selectAndExecuteEvent(), PAPRECA::deleteAndClearLocalEvents(), PAPRECA::equilibrate(), PAPRECA::finalize()
	/// @note Example execution of PAPRECA from UNIX terminal: mpiexec papreca -in in_kmc.lmp in_kmc.ppc. CAUTION: Always provide the LAMMPS input file first and the PAPRECA input file second, otherwise the code will exit with an error. 
	/// @note See paper this paper for more information regarding the classic N-FOLD way and the event selection process: https://www.sciencedirect.com/science/article/pii/S0927025623004159
	if ( narg != 4 ) { allAbortWithMessage( MPI_COMM_WORLD , "Syntax Error. Input command should be in the following form: mpirun -np N main -in in.lammps in.papreca." );}
	
	//Declare variables
	PaprecaConfig papreca_config;
	int proc_id , nprocs;
	LAMMPS *lmp = NULL;
	
	initialize( &narg , &arg , &nprocs , &proc_id , &lmp , papreca_config );
	
	//kMC run preparation
	char event_type[50] = "NONE";
	double time = 0.0 , film_height = 0.0;
	int zero_rate = 0; //Tells you whether the total rate on a given step is zero (if it is you should equilibrate immediately)
	vector<Event*> events_local; //using a unique pointer to switch effortlessly between the child classes of Event
	events_local.reserve( 10 ); //Reserve 10 events per proc. Obviously, the vector will resize if necessary 
	ATOM2BONDS_MAP atomID2bonds;
	
	//Main loop
	for( int i = 1; i <= papreca_config.getKMCsteps( ); ++i ){
		
		//Initial timestamp for execution time measurement
		papreca_config.setHybridStartTimeStamp4ExecTimeFile( i );
		
		//Init atomID2bonds 
		PAPRECA::Bond::initAtomID2BondsMap( lmp , proc_id , atomID2bonds );
		
		//KMC Operations
		loopAtomsAndIdentifyEvents( lmp , proc_id , nprocs , i , papreca_config , events_local , atomID2bonds , film_height );
		zero_rate = selectAndExecuteEvent( lmp , i , time , event_type , proc_id , nprocs , papreca_config , events_local , atomID2bonds , film_height );
		Event::deleteAndClearLocalEvents( lmp , events_local );	
		
		//LAMMPS Equilibration
		equilibrate( lmp , proc_id , nprocs , time , papreca_config , film_height , zero_rate , i , atomID2bonds );
		
		//Reset atomID2bonds
		atomID2bonds.clear( );
		
		//Export Files
		papreca_config.dumpLAMMPSRestart( lmp , i );
		papreca_config.calcHybridAndKMCTimes4ExecTimeFile( nprocs , i );
		papreca_config.appendExportFiles( lmp , proc_id , time , event_type , film_height , i );
		
		//Test if target ending time is exceeded to exit prematurely. Ensure that 
		if( time >= papreca_config.getTimeEnd( ) ){ break; }
		
	}
	
	finalize( &lmp , papreca_config , proc_id );	
}
