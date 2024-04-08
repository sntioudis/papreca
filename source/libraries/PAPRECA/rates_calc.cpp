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
///@brief Function definitions for rates_calc.h

#include "rates_calc.h"


namespace PAPRECA{

	double getRateFromArrhenius( const double &activation_energy , const double &attempt_freq , const double &temperature ){
	
		/// This function uses the classic Hertzian equation to output the rate of a rare event.
		/// @param[in] activation_energy Activation energy in [kcal/mol] units.
		/// @param[in] attempt_freq Attempt frequency in [Hz (1/s)] units.
		/// @param[in] temperature in [K] units.
		/// @return rate in [Hz] (1/s) units.
		/// @note Hardcoded parameters: Universal gas constant in [kcal/(mol*K)] units.
		
		double gas_constant = 1.98720425864083e-3; //added e-3 in the end to convert from cal/(mol*K) to kcal/(mol*K);
		return attempt_freq * std::exp( -activation_energy / ( gas_constant * temperature ) );
	}
	
	double getDepoRateFromHertzKnudsen( const double &pressure_in , const double &ads_area_in , const double &ads_mass_in , const double &temperature_in ){
		
		/// This function outputs the deposition base rate using the Hertz-Knudsen equation (kinetic theory of gases). 
		/// @param[in] pressure_in Partial Pressure in [Bar] units.
		/// @param[in] ads_area_in Adsorption Site Area in [Angstrom^2] units.
		/// @param[in] ads_mass_in Molecule mass in [g/mol] units.
		/// @param[in] temperature_in temperature in [K].
		/// @return rate in [Hz] (1/s) units
		/// @note See: https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Surface_Science_(Nix)/02%3A_Adsorption_of_Molecules_on_Surfaces/2.03%3A_Kinetics_of_Adsorption.
		/// @note Disclaimer: This function assumes that the sticking coefficients is 1. However, you can still use that function even if the sticking coefficient is non-unity or if the sticking coefficient is a function of surface coverage. Simply multiply the rate obtained from the function with your sticking coefficient of choice. This is automatically done in the PaprecaConfig class.		
		
		double c_boltzmann = 1.380649e-23; //[J/K]
		double c_avocadro = 6.02214076e23;
		double pi = 3.14159265359;
		
		const double pressure = 1.0e5 * pressure_in; //convert bar to Pa (kg/(m*s^2))
		const double ads_area = 1.0e-20 * ads_area_in; //convert angstroms^2 to m^2;
		const double ads_mass = 1.0e-3 * ads_mass_in / c_avocadro; //convert g/mol to kg
		
		double rate = ( pressure * ads_area ) / std::sqrt( 2 * pi * ads_mass * c_boltzmann * temperature_in );		
	
		return rate;
		
	}
	
	double getDesorptionRate( const double &activation_energy , const double &temperature ){
		
		/// Calculates a desorption rate based on an Arrhenius-type equation.
		/// @param[in] activation_energy Activation Energy in [kcal/mol] units.
		/// @param[in] temperature Temperature in [K] units.
		/// @return rate in [Hz] (1/s) units
		/// @note See: https://pubs.acs.org/doi/epdf/10.1021/acs.jpcc.8b06909.
		/// @note The pre-exponential is calculated from the boltzmann and plank constants, at a given temperature. The exponential term can be calculated from the universal gas constant.
		
		double c_boltzmann = 1.380649e-23; //[J/K]
		double c_plank = 6.62607015e-34; //J*s
		double gas_constant = 1.98720425864083e-3; //added e-3 in the end to convert from cal/(mol*K) to kcal/(mol*K);
		
		double rate = ( (c_boltzmann * temperature ) / c_plank ) * std::exp( -activation_energy / ( gas_constant *temperature ) );
		
		return rate;
		
	}
	
} //End of namespace PAPRECA
