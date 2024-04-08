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
///@brief Definitions for mpi_wrappers.h

#include "mpi_wrappers.h"
	
namespace PAPRECA{
	
		void setupMPI( int *narg , char ***arg , int *nprocs , int *proc_id ){
		
		/// Initializes MPI environment using command-line arguments passed during the program invocation from the terminal.
		/// @param[in] narg number of command-line arguments.
		/// @param[in] arg array containing command-line arguments.
		/// @param[in,out] nprocs number of MPI processes.
		/// @param[in,out] proc_id identifier of current MPI process.
		
		MPI_Init( narg , arg );
		MPI_Comm_rank( MPI_COMM_WORLD , proc_id );
		MPI_Comm_size( MPI_COMM_WORLD , nprocs );
	}

	const int getMPIRank( MPI_Comm communicator ){
			
		/// Wrapper around MPI_Comm_rank. Returns the MPI process ID of the specific process calling the function.
		/// @param[in] communicator MPI communicator.
		/// @return ID of MPI process.
		
		int proc_id;
		MPI_Comm_rank( communicator , &proc_id );
		return proc_id;
		
	}
	
}//End of PAPRECA namespace