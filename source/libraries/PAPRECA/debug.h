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
///@brief Functions that can potentially be useful for debugging purposes. These function are not used in the main code. Also, some of those might need a modification to work for your specific debugging requirements.

#ifndef DEBUG_H
#define DEBUG_H

//System Headers
#include <unordered_map>
#include <unordered_set>
#include <iostream>

//PAPRECA headers
#include "bond.h"
#include "event.h"
#include "event_list.h"

//LAMMPS headers
#include "lammps.h"
/// \cond
#include "library.h"
#include "pointers.h"
#include "atom.h"
/// \endcond


namespace PAPRECA{
	

	
	void debugPrintBondMapPairs( ATOM2BONDS_MAP const &bonds_map , const int &proc_num );
	void debugPrintBasicAtomInfo( LAMMPS_NS::LAMMPS *lmp , const int &proc_id );
	void debugPrintNeighborLists( LAMMPS_NS::LAMMPS *lmp , const int &proc_id );
	void debugCheckBondsInNeibLists( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , ATOM2BONDS_MAP &atomID2bonds );
	void debugPrintBondsList( LAMMPS_NS::tagint *bonds_list , LAMMPS_NS::bigint &bonds_num , const int &proc_id );
	void debugPrintType2SigmaMap( INTPAIR2DOUBLE_MAP &types2sigma );
	void debugPrintEventInfo( Event *event , const int &proc_id );
	void debugCheckDeposition( );
	
} //End of PAPRECA namespace




#endif
