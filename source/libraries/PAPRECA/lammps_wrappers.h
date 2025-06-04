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
/// @brief Declarations for LAMMPS wrapper functions. 
///
/// Typically functions that calls lmp->input->one or lammps_command to execute a specific LAMMPS command with given inputs.

#ifndef LAMMPS_WRAPPERS_H
#define LAMMPS_WRAPPERS_H

//System Headers
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <array>
#include <cmath>

//LAMMPS headers
#include "lammps.h"
/// \cond
#include "atom.h"
#include "molecule.h"
#include "library.h"
#include "input.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
/// \endcond

//kMC Headers
#include "papreca_config.h"
#include "papreca_error.h"
#include "utilities.h"


namespace PAPRECA{
	class PaprecaConfig; //Forward declaration to resolve compiler errors due to circular dependencies
	
	//Initialize LAMMPS
	void initializeLMP( LAMMPS_NS::LAMMPS **lmp );
	void readLMPinput( const std::string &lmp_input , LAMMPS_NS::LAMMPS *lmp );
		
	//Setup integrators
	void setupNveLimIntegrator( LAMMPS_NS::LAMMPS *lmp ,  PaprecaConfig &papreca_config );
	void removeNveLimIntegrator( LAMMPS_NS::LAMMPS *lmp ,  PaprecaConfig &papreca_config );
	
	//Execute LAMMPS
	void runLammps( LAMMPS_NS::LAMMPS *lmp , const int &timesteps_num );
	void MPIBcastAndExecuteCommand( LAMMPS_NS::LAMMPS *lmp , std::string &command ); //This function gets a line command (std::string), casts it to all other procs, and executes the command
	
	//Period Box Operations
	void remap3DArrayInPeriodicBox( LAMMPS_NS::LAMMPS *lmp , double *arr ); //Receives a 3D array and remaps it inside the existing periodic box.
	
	//kMC operations
	void deleteAtoms( LAMMPS_NS::LAMMPS *lmp , LAMMPS_NS::tagint *atom_ids , const int &num_atoms , const std::string &delete_bonds , const std::string &delete_molecule );
	void deleteAtoms( LAMMPS_NS::LAMMPS *lmp , std::vector< LAMMPS_NS::tagint > &atom_ids , const std::string &delete_bonds , const std::string &delete_molecule ); //Overloaded function of the deleteAtoms function to work with vectors
	void deleteAtomsInBoxRegion( LAMMPS_NS::LAMMPS *lmp , double &boxxlo , double &boxxhi , double &boxylo , double &boxyhi , double &boxzlo , double &boxzhi , const std::string &delete_bonds , const std::string &delete_molecule );
	void createAtom( LAMMPS_NS::LAMMPS *lmp , const double atom_pos[3] , const int &atom_type );
	void deleteBond( LAMMPS_NS::LAMMPS *lmp , const LAMMPS_NS::tagint &atom1id , const LAMMPS_NS::tagint &atom2id , const bool special );
	void formBond( LAMMPS_NS::LAMMPS *lmp , const LAMMPS_NS::tagint &atom1id , const LAMMPS_NS::tagint &atom2id , const int &bond_type );
	void resetMobileAtomsGroups( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config );
	void insertMolecule( LAMMPS_NS::LAMMPS *lmp , const double site_pos[3] , const double rot_pos[3] , const double &rot_theta , const int &mol_id , const char *mol_name );
	void diffuseAtom( LAMMPS_NS::LAMMPS *lmp , const double vac_pos[3] , const LAMMPS_NS::tagint &parent_id , const int &parent_type , const int &is_displacive , const int &diffused_type );
	
	//Sigmas
	void initType2SigmaFromLammpsPairCoeffs( LAMMPS_NS::LAMMPS *lmp , INTPAIR2DOUBLE_MAP &type2sigma );
	
	//Neibs lists
	int getMaskedNeibIndex( int *neighbors , int &j );
	
	//Bond Lists
	void initAndGatherBondsList( LAMMPS_NS::LAMMPS *lmp , LAMMPS_NS::tagint **bonds_list , LAMMPS_NS::bigint &bonds_num );
	
	//Molecules
	static inline const int getMolIndexFromMolName( LAMMPS_NS::LAMMPS *lmp , std::string mol_name );
	void computeMolCenter( LAMMPS_NS::LAMMPS *lmp , std::string mol_name );
	
	//Files
	void dumpRestart( LAMMPS_NS::LAMMPS *lmp , const int &KMC_loopid , const int &dump_freq );
	
	//LAMMPS Maths wrappers
	double get3DSqrDistWithPBC( LAMMPS_NS::LAMMPS *lmp , const double *x1 , const double *x2 );
	
	
}//end of namespace PAPRECA


#endif
