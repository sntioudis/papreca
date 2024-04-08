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
/// @brief Definitions for sim_clock.h

#include "sim_clock.h"

namespace PAPRECA{

	void advanceSimClockFromKMC( PaprecaConfig &papreca_config , const double &proc_rates_sum , double &time ){

		/// Stochastically advances the simulation clock by a KMC step.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] proc_rates_sum total (gathered/sum from all MPI processes) event rate.
		/// @param[in,out] time time on current PAPRECA step.
		/// @note This function should not be confused with PAPRECA::advanceSimClockFromLAMMPS(). PAPRECA::advanceSimClockFromLAMMPS() advances the clock forward by timestep*dt (LAMMPS variables).
		/// @note See this paper for more information regarding the classic N-FOLD way and the advancement of the simulation clock: https://www.sciencedirect.com/science/article/pii/S0927025623004159
		
		double rnum = papreca_config.getUniformRanNum( );
		double dt = -log( rnum ) / proc_rates_sum;
		time += dt;

		
	}
	
	void advanceSimClockFromLAMMPS( PaprecaConfig &papreca_config , double &time ){ 
		
		/// Advances the simulation clock by timestep*dt (as defined in the LAMMPS and PAPRECA inputs).
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in,out] time current time.
		/// @note This function should not be confused with PAPRECA::advanceSimClockFromKMC(). PAPRECA::advanceSimClockFromKMC() advances the clock forward in the KMC/N-FOLD way.
		
		time += papreca_config.getCtimeConvert( ) * papreca_config.getTrajDuration( );
		
	}
	
} //End of PAPRECA Namespace
