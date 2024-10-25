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
/// @brief Declarations for functions that execute events
///

#ifndef EVENT_EXECUTE_H
#define EVENT_EXECUTE_H

//System Headers
#include <vector>
#include <mpi.h>


//LAMMPS headers
#include "lammps.h"
/// \cond
#include "pointers.h"
/// \endcond

//KMC headers
#include "bond.h"
#include "event.h"
#include "event_list.h"
#include "event_select.h"
#include "sim_clock.h"
#include "papreca_config.h"
#include "papreca_error.h"
#include "utilities.h"

namespace PAPRECA{
	
	//Formation events
	void fillFormTransferDataArr( BondForm *bond_form , int *form_data );
	void deserializeFormTransferDataArr( int *form_data , int &bond_type , int &delete_atoms );
	void executeBondForm( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &KMC_loopid , double &time , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event );
	
	//Bond-breaking events
	void executeBondBreak( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &KMC_loopid , double &time , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event , ATOM2BONDS_MAP &atomID2bonds );
	
	//Deposition events
	void fillDepoDataTransfArr( double *depo_data , Deposition *depo );
	void deserializeDepoTransfData( double *depo_data , double *site_pos , double *rot_pos , double &rot_theta , double &insertion_vel );
	void executeDeposition( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &KMC_loopid , double &time , PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event );
	
	//Diffusion events
	void fillIntegerDiffDataTransfArray( int *diff_intdata , Diffusion *diff );
	void fillDoubleDiffDataTransfArray( double *diff_doubledata , Diffusion *diff );
	void deserializeIntegerDiffDataArr( int *diff_intdata , int &parent_type , int &is_displacive , int &diffused_type );
	void deserializeDoubleDiffDataArr( double *diff_doubledata , double *vac_pos , double &insertion_vel );
	void executeDiffusion( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &KMC_loopid , double &time , PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event );
	
	//Monoatomic desorption events
	void executeMonoatomicDesorption( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &KMC_loopid , double &time , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event );
	
	//General event execution functions
	void printStepInfo( PaprecaConfig &papreca_config , const int &KMC_loopid , const double &time , const double &film_height , const double &proc_rates_sum );
	void executeEvent( LAMMPS_NS::LAMMPS *lmp , int &KMC_loopid , double &time , PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs , const int &event_proc , const int &event_num , char *event_type , std::vector< Event* > &events_local , ATOM2BONDS_MAP &atomID2bonds );
	int selectAndExecuteEvent( LAMMPS_NS::LAMMPS *lmp , int &KMC_loopid , double &time , char *event_type , int &proc_id , int &nprocs , PaprecaConfig &papreca_config , std::vector< Event* > &events_local , ATOM2BONDS_MAP &atomID2bonds , double &film_height );

}//end of PAPRECA namespace 


#endif
