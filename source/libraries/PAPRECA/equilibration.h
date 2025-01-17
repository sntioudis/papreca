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
/// @brief Declarations for functions equilibrate the system (i.e., Molecular Dynamics runs, system minimization runs) and equilibration-related functions (e.g., functions that delete desorbed atoms).
///

#ifndef EQUILIBRATION_H
#define EQUILIBRATION_H

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
#include "papreca_config.h"
#include "lammps_wrappers.h"
#include "sim_clock.h"
#include "utilities.h"

namespace PAPRECA{

	//Delete desorbed atoms
	void fillDelidsLocalVec( LAMMPS_NS::LAMMPS *lmp , const double &desorb_cut , std::vector< LAMMPS_NS::tagint > &delids_local , ATOM2BONDS_MAP &atomID2bonds );
	bool delidsLocalVectorsAreEmpty( std::vector< LAMMPS_NS::tagint > &delids_local );
	void gatherAndTrimDelIdsOnDriverProc( const int &proc_id , const int &nprocs , std::vector< LAMMPS_NS::tagint > &delids_local , std::vector< LAMMPS_NS::tagint > &delids_global );
	int fillDelidsVec( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const double &desorb_cut , std::vector< LAMMPS_NS::tagint > &delids , ATOM2BONDS_MAP &atomID2bonds );
	void broadcastDelidsFromMasterProc( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , int &delids_num , std::vector< LAMMPS_NS::tagint > &delids );
		
	//Equilibration
	void equilibrateFluidAtoms( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , double &time , const unsigned long int &trajectory_duration);
	void equilibrate( LAMMPS_NS::LAMMPS *lmp , int &proc_id , const int &nprocs , double &time , PaprecaConfig &papreca_config , double &film_height , int &zero_rate , const int &KMC_loopid , ATOM2BONDS_MAP &atomID2bonds );
	
}//end of PAPRECA namespace 


#endif
