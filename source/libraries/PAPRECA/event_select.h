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
/// @brief Declarations for functions that (stochastically) select events for execution
///

#ifndef EVENT_SELECT_H
#define EVENT_SELECT_H

//System Headers
#include <vector>
#include <mpi.h>


//LAMMPS headers
#include "lammps.h"
/// \cond
#include "pointers.h"
/// \endcond

//KMC headers
#include "event.h"
#include "event_list.h"
#include "papreca_config.h"
#include "utilities.h"

namespace PAPRECA{

	double getLocalRate( std::vector< Event* > &events_local , PaprecaConfig &papreca_config );
	void fillAndSortIndexedRatesVec( double *arr , const int &arr_size , DOUBLE2INTPAIR_VEC &rates_indexed );
	int selectProcessStochastically( double *arr , const int &arr_size , double &rnum , double &rates_sum );
	
}//end of PAPRECA namespace 


#endif
