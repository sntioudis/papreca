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
///@brief Header for geometry calculation functions (e.g., film height, interference).

#ifndef GEOMETRY_CALC_H
#define GEOMETRY_H

//System Headers
#include <mpi.h>
#include <vector>

//LAMMPS headers
#include "lammps.h"
/// \cond
#include "domain.h"
/// \endcond

//kMC Headers
#include "papreca_error.h"
#include "papreca_config.h"

namespace PAPRECA{
	
	//Film Height calculation
	void calcLocalMassAndFillMassProfile( LAMMPS_NS::LAMMPS *lmp , double **mass_profiles , double &local_mass , const int &atom_type , double *atom_xyz , const double &atom_mass , const double &bin_width , const int &bins_num );
	double **initMassProfilesArr( const int &types_num , const int &bins_num );
	void deleteMassProfilesArr( double **mass_profiles , const int &bins_num  );
	void fillMassProfilesTotalArrFromMassProfilesLocal( const int &bins_num , const int &types_num , double **mass_profiles_total , double **mass_profiles_local );
	void getFilmHeightFromMassBinsMethod( PaprecaConfig &papreca_config , LAMMPS_NS::LAMMPS *lmp , const int &proc_id , double &film_height , double **mass_profiles_total , const double &local_mass , double *atom_mass , const int &bins_num , const int &types_num , const double &bin_width );
	void calcFilmHeight( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const int &KMC_loopid , PaprecaConfig &papreca_config , double &film_height );
	
	//Interference between atoms
	const bool atomsCollide( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , double *atom1_xyz , const int &atom1_type , double *atom2_xyz , const int &atom2_type );
	
}//end of namespace PAPRECA

#endif
