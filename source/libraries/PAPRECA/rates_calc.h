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
/// @brief Includes functions used to calculate the rates of predefined events. All functions return rates in Hz (i.e., 1/s). However, not all inputs to these functions are in SI units. Please see the specific function descriptions for more information. PAPRECA reports time in SI units (i.e., seconds), however, all other units are in complete accordance with LAMMPS units (e.g., REAL, SI, etc.) as defined in the LAMMPS units file (through the units command).

#ifndef RATES_CALC_H
#define RATES_CALC_H

//system headers
#include <cmath>

namespace PAPRECA{

	double getRateFromArrhenius( const double &activation_energy , const double &attempt_freq , const double &temperature );
	double getDepoRateFromHertzKnudsen( const double &pressure , const double &ads_area , const double &ads_mass , const double &temperature );
	double getDesorptionRate( const double &activation_energy , const double &temperature );
	
}//end of PAPRECA namespace

#endif
