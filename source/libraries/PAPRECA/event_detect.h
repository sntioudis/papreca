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
/// @brief Declarations for functions that detect predefined events
///

#ifndef EVENT_DETECT_H
#define EVENT_DETECT_H

//System Headers
#include <vector>
#include <mpi.h>


//LAMMPS headers
#include "lammps.h"
/// \cond
#include "pointers.h"
/// \endcond

//KMC headers
#include "event.h"
#include "event_list.h"
#include "bond.h"
#include "papreca_config.h"
#include "lammps_wrappers.h"
#include "geometry_calc.h"
#include "utilities.h"

namespace PAPRECA{

	//Diffusion events
	const bool feCandidateHas4PO4Neibs( PaprecaConfig &papreca_config , PredefinedDiffusionHop *diff_template , LAMMPS_NS::tagint *atom_ids , int *atom_types , int *neighbors , int &neighbors_num , ATOM2BONDS_MAP &atomID2bonds );
	void getDiffPointCandidateCoords( LAMMPS_NS::LAMMPS *lmp  , PaprecaConfig &papreca_config , double *iatom_xyz , double *candidate_xyz , PredefinedDiffusionHop *diff_template );
	const bool candidateDiffHasCollisions( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int *neighbors , int &neighbors_num , double *candidate_xyz , const int &diffused_type , double *iatom_xyz , const int &iatom_type );
	void getDiffEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &iatom , int *neighbors , int &neighbors_num , std::vector< Event* > &events_local , ATOM2BONDS_MAP &atomID2bonds );
	
	//Deposition events
	const bool atomIsInDepoScanRange( PaprecaConfig &papreca_config , double *iatom_xyz , double &film_height );
	void getDepoPointCandidateCoords( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , double *iatom_xyz ,  double *candidate_xyz , PredefinedDeposition *depo_template );
	const bool depoCandidateIsBelowRejectionHeight( PaprecaConfig &papreca_config , double *candidate_xyz , const double &film_height );
	void getMolCoords( LAMMPS_NS::LAMMPS *lmp , double **mol_xyz , double **mol_dx , const int &mol_natoms , double *candidate_center );
	void initMolCoordsArr( double ***mol_xyz , const int &mol_natoms );
	void deleteMolCoordsArr( double **mol_xyz , const int &mol_natoms );
	bool atomHasCollisionWithMolAtoms( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , double *atom_xyz , const int &atom_type , const int &mol_natoms , double **mol_xyz , int *mol_atomtype );
	bool candidateDepoHasCollisions( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const int &nprocs , PaprecaConfig &papreca_config , int *neighbors , int neighbors_num , double *candidate_center , double *iatom_xyz , const int &iatom_type , PredefinedDeposition *depo_template );
	void getDepoEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs , const int &iatom , int *neighbors , int &neighbors_num , double &film_height , std::vector< Event* > &events_local );
	
	//Bond-Breaking and formation events
	const bool headAtomIsCatalyzed( PredefinedReaction *reaction_template , int *atom_types , int *neighbors , int &neighbors_num );
	void getBondBreakingEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &iatom , int *neighbors , int &neighbors_num , std::vector<Event*> &events_local , ATOM2BONDS_MAP &atomID2bonds );
	const bool atomsBelong2TheSameMol( const LAMMPS_NS::tagint &iatom_mol , const LAMMPS_NS::tagint &jneib_mol );
	const bool atomHasMaxBonds( PaprecaConfig &papreca_config , ATOM2BONDS_MAP &atomID2bonds , const LAMMPS_NS::tagint &atom_id , const int atom_type );
	bool bondBetweenAtomsExists( ATOM2BONDS_MAP &atomID2bonds , const LAMMPS_NS::tagint &atom1_id , const LAMMPS_NS::tagint &atom2_id );
	const bool atomCandidatesAreLone( const LAMMPS_NS::tagint atom1_id , const LAMMPS_NS::tagint atom2_id , ATOM2BONDS_MAP &atomID2bonds );
	const bool atomHasMaxBondTypes( PaprecaConfig &papreca_config , ATOM2BONDS_MAP &atomID2bonds , const LAMMPS_NS::tagint &atom_id , const int &atom_type , const int &bond_type );
	void getBondFormEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &iatom , int *neighbors , int &neighbors_num , std::vector<Event*> &events_local , ATOM2BONDS_MAP &atomID2bonds );
	
	//Monoatomic Desorption events
	void getMonoDesEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &iatom , std::vector< Event* > &events_local , ATOM2BONDS_MAP &atomID2bonds );
	
	//General Functions
	void loopAtomsAndIdentifyEvents( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , int &nprocs , const int &KMC_loopid , PaprecaConfig &papreca_config , std::vector<Event*> &events_local , ATOM2BONDS_MAP &atomID2bonds , double &film_height );
		

}//end of PAPRECA namespace 


#endif
