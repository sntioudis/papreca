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
/// @brief Definitions for event_detect.h.

#include "event_select.h"

namespace PAPRECA{
	
	double getLocalRate( std::vector< Event* > &events_local , PaprecaConfig &papreca_config ){
	
		/// Calculates the total local (on the current MPI process) rate by summing the rates of the detected local events (contained in the PAPRECA::Event objects vector, events_local).
		/// @param[in] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @return local (per MPI process) rate.
		/// @note This function also dumps surface coverage file (if activated by the user).
		
		if( papreca_config.getSurfaceCoverageFile( ).isActive( ) ){ papreca_config.calcSurfaceCoverage( ); } //ALWAYS calculate the surface coverage BEFORE calculating the sticking coefficients. Do that because the calcStickingCoeffs function zeroes the depsites and deptries and those variables are neccesary to calcualte the surface coverage.
		papreca_config.calcStickingCoeffs( );
		
		double rate_local = 0.0;
		for( const auto &event : events_local ){
			
			//For deposition events, scale event rate by sticking coefficient (the sticking coefficient of all events is initialized to 1 and the tweaked in this function).
			if( event->getType( ) == "DEPO" ){
				Deposition *depo = dynamic_cast< Deposition* >( event );	
				event->setRate( event->getRate( ) * depo->getDepoTemplate( )->getStickingCoeff( ) ); //For non-variable sticking coefficient events the getStickingCoeff function will return the constant sticking coeff. For variable sticking coefficients, the relevant value is obtained from the calcStickingCoeffs function of the PaprecaConfig class.
			}
			
			rate_local += event->getRate( );	
		}
		
		return rate_local;
		
	}
	
	void fillAndSortIndexedRatesVec( double *arr , const int &arr_size , DOUBLE2INTPAIR_VEC &rates_indexed ){
		
		/// Uses std::sort to sort a double array. Typically used to sort rates vector in ascending order. Sorting helps with faster event selection. It also ensures repeatability (if the same random number seed is used), because the event enumeration is maintained even if the domain decomposition numbering changes (i.e., if MPI decides to decompose the domain differently). The sorting is not performed on the array itself but on an auxiliary std::vector< std::pair< double , int > > containers that helps us sort while maintaining the index of events. This is necessary, since each time we stochastically select an event using the rates array, we have to link it (by index) to the correct PAPRECA::Event.
		/// @param[in] arr double array container.
		/// @param[in] arr_size size of arr.
		/// @param[in,out] rates_indexed std::vector< std::pair< double , int > > container used for ascending sorting.
		
		rates_indexed.reserve( arr_size );
		
		for( int i = 0; i < arr_size; ++i ){
			rates_indexed.push_back( std::make_pair( arr[i] , i ) );
		}
		


		std::sort( rates_indexed.begin( ) , rates_indexed.end( ) );
	}
	
	int selectProcessStochastically( double *arr , const int &arr_size , double &rnum , double &rates_sum ){
		
		/// Selects a process stochastically (through the classic N-FOLD way selection process). Can be used to select an MPI proc that fires an event or can be used to select the next kMC event inside the proc. The event selection is performed on a std::vector< std::pair< double , int > > container that holds the rate and the associated index of events.
		/// @param[in] arr array of rates.
		/// @param[in] arr_size size of arr.
		/// @param[in] rnum uniformly distributed (between 0 and 1) pseudorandom number (usually drawn on master proc).
		/// @param[in] rates_sum total sum of rates
		/// @return integer that corresponds to the ith position of the provided vector (i.e., index of the selected process) or -1 if the function could not select a process (indicates error).
		/// @see PAPRECA::selectAndExecuteEvent()
		/// @note See paper and equation 3 here for more information regarding the classic N-FOLD way and the event selection process: https://www.sciencedirect.com/science/article/pii/S0927025623004159
		
		DOUBLE2INTPAIR_VEC rates_indexed; //We will be sorting the rates vec so we need a vector of <double,int> pairs (i.e., rate and proc num). This helps us know the proc num AFTER sorting the processes by rate.
		//The rates vec is either a proc rates vec (index is the proc_num) or an event rates vec (index is the event num).
		fillAndSortIndexedRatesVec( arr , arr_size , rates_indexed );
		
		rnum *= rates_sum; //Scale rnum by the sum of rates (could introduce floating point errors).
		double rate_cur = 0.0;
		
		for( int i = 0; i < arr_size; ++i ){

			double rate = rates_indexed[i].first;
			int index = rates_indexed[i].second;
			
			if( rate > 0 ){ //Avoid processes (or procs) with zero rate.
			
				rate_cur += rate;
				if( rnum <= rate_cur ){
					return index;
					
				}
			}
			
			
			
		}
		
		allAbortWithMessage( MPI_COMM_WORLD , "No event was selected in function selectProcessStochastically in papreca.cpp" );
		return -1; //This means that some error occurred and no process was selected
		
	}

} //End of PAPRECA Namespace
