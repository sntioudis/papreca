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
/// @brief Declarations for functions that (stochastically) advance the simulation time
///

#ifndef SIM_CLOCK_H
#define SIM_CLOCK_H


//LAMMPS headers
#include "lammps.h"
/// \cond
#include "pointers.h"
/// \endcond

//KMC headers
#include "papreca_config.h"


namespace PAPRECA{
	
	void advanceSimClockFromKMC( PaprecaConfig &papreca_config , const double &proc_rates_sum , double &time );
	void advanceSimClockFromLAMMPS( PaprecaConfig &papreca_config , double &time );
	
}//end of PAPRECA namespace 


#endif
