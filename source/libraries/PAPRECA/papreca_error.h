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
///@brief Functions that enable the communication of error/warning messages to the terminal and coordinate the safe termination of MPI processes.

#ifndef PAPRECA_ERROR_H
#define PAPRECA_ERROR_H


//System Headers
#include <iostream>
#include <string>
#include <mpi.h>

//PAPRECA Headers
#include "mpi_wrappers.h"


namespace PAPRECA{
	
	void warnOne( MPI_Comm communicator , const std::string &message );
	void warnAll( MPI_Comm communicator , const std::string &message );
	void allAbort( MPI_Comm communicator );
	void allAbortWithMessage( MPI_Comm communicator , const std::string &message );
	
}//end of namespace PAPRECA


#endif