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
///@brief Definitions for papreca_error.h

#include "papreca_error.h"

namespace PAPRECA{
	
	void warnOne( MPI_Comm communicator , const std::string &message ){
		
		/// Throws a warning on a specific MPI process.
		/// @param[in] communicator MPI communicator.
		/// @param[in] message std::string of the warning message.
		
		std::cout << "PAPRECA WARNING on proc " << getMPIRank( MPI_COMM_WORLD ) << " " << message << std::endl;
		
	}
	
	void warnAll( MPI_Comm communicator , const std::string &message ){
		
		/// Throws a warning to all MPI processes.
		/// @param[in] communicator MPI communicator.
		/// @param[in] message std::string of the warning message.
		
		int proc_id = getMPIRank( MPI_COMM_WORLD );
		if( proc_id == 0 ){ std::cout << "PAPRECA WARNING! " << message << std::endl; }
		
	}
	
	void allAbort( MPI_Comm communicator ){
		
		/// Aborts all MPI processes.
		/// @param[in] communicator communicator
		
		int proc_id = getMPIRank( MPI_COMM_WORLD );
		if( proc_id == 0 ){ std::cout << "FATAL PAPRECA ERROR! Code exited with an error. Look for warnings to understand what went wrong." << std::endl; }
		
		MPI_Finalize( );
		exit(1);
		
	}
	
	
	void allAbortWithMessage( MPI_Comm communicator , const std::string &message ){
		
		/// Aborts all MPI processes and throws a message on the screen.
		/// @param[in] communicator MPI communicator.
		/// @param[in] message std::string of the error message.
		
		int proc_id = getMPIRank( MPI_COMM_WORLD );;
		
		if( proc_id == 0 ){ std::cout << "FATAL PAPRECA ERROR! " << message << std::endl; }
		
		MPI_Finalize( );
		exit(1);
		
	}
	
}