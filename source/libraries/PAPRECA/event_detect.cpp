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

#include "event_detect.h"

namespace PAPRECA{

	//Diffusion events
	const bool feCandidateHas4PO4Neibs( PaprecaConfig &papreca_config , PredefinedDiffusionHop *diff_template , LAMMPS_NS::tagint *atom_ids , int *atom_types , int *neighbors , int &neighbors_num , ATOM2BONDS_MAP &atomID2bonds ){
	
		/// Scans neighbors of a diffusion parent atom to check whether 4 distinct PO4 (phosphate) structures exist in the neighbors list. This function is only called when the custom PAPRECA::PredefinedDiffusionHop style: Fe_4PO4neib is active.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] diff_template Diffusion Hop template (PAPRECA::PredefinedDiffusionHop) as initialized by the user (in the PAPRECA input file).
		/// @param[in] atom_ids LAMMPS atom ids.
		/// @param[in] atom_types LAMMPS atom types.
		/// @param[in] neighbors LAMMPS array storing the IDs of neighbors of the current (parent) atom type.
		/// @param[in] neighbors_num total number of neighbors.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @return true/false depending on whether at least 4 PO4 (phosphate) structures exist or not in the neighborhood of the diffusion parent candidate atom.
		/// @see PAPRECA::getDiffEventsFromAtom()
		
		int po4_neibs_num = 0;
		
			
		for( int i = 0; i < neighbors_num; ++i ){
				
			int ineib = getMaskedNeibIndex( neighbors , i );
			int neib_type = atom_types[ineib];
			LAMMPS_NS::tagint neib_id = atom_ids[ineib];
				
			if( ( diff_template->getStyleAtomTypes( ) )[0] == neib_type ){ //the atom style types vector of custom Fe_4PO4neib diffusion events contains the P atom (we use its bond list to identify the PO4 neibs).
				
				BOND_VECTOR &bonds = atomID2bonds[ neib_id ]; //Get bonds of P neib
				if( bonds.size( ) == papreca_config.getMaxBondsFromSpecies( ( diff_template->getStyleAtomTypes( ) )[0] ) ){ 
					++po4_neibs_num;	
					if( po4_neibs_num >= 4 ){ return true; } //immediately return if you find at least 4 PO4 neibs
						
						
				}
				
			}
					
			
		}
		
		return false; //If you scanned all neighbors but there were no 4 po4 structures exit with false.
		

	}
	
	void getDiffPointCandidateCoords( LAMMPS_NS::LAMMPS *lmp  , PaprecaConfig &papreca_config , double *iatom_xyz , double *candidate_xyz , PredefinedDiffusionHop *diff_template ){
	
		/// Calculates the diffusion point (vacancy) coordinates for a given parent atom. Depending on the user settings, the diffusion point can be directly above the parent atom or at the surface of a sphere centered at the coordinates of the parent atom.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] iatom_xyz coordinates of parent atom.
		/// @param[in,out] candidate_xyz coordinates of diffusion point.
		/// @param[in] diff_template Diffusion template (PAPRECA::PredefinedDiffusionHop) as initialized by the user (in the PAPRECA input file).
		/// @see PAPRECA::getDiffEventsFromAtom()
		
		//For definition of random_vacancy see getDepoPointCandidateCoords function...
		
		//Get Diffusion distance
		double diff_dist = diff_template->getDiffusionDist( );
		
		if( papreca_config.diffVecsAreRandom( ) ){
			
			const double rnum1 = papreca_config.getUniformRanNum( );
			const double rnum2 = papreca_config.getUniformRanNum( );
			
			double theta = 2.0 * M_PI * rnum1; //Gives a number between 0 and 2PI
			double phi = 0.0;
			
			if( papreca_config.getRandomDiffVecsStyle( ) == "2D" ){
				phi = 0.5 * M_PI * rnum2; //Gives a number between 0 and pi/2. This means that the random diffvec can only be above the parent type.
			}else if( papreca_config.getRandomDiffVecsStyle( ) == "3D" ){
				phi = M_PI * rnum2;
			}else{
				allAbortWithMessage( MPI_COMM_WORLD , "Unkown random diffvecs style " + papreca_config.getRandomDiffVecsStyle( ) );
			}
			
			candidate_xyz[0] = iatom_xyz[0] + diff_dist * sin( phi ) * cos( theta );
			candidate_xyz[1] = iatom_xyz[1] + diff_dist * sin( phi ) * sin( theta );
			candidate_xyz[2] = iatom_xyz[2] + diff_dist * cos( phi );
				
		}else{
			
			//Same equations as random vectors case (if phi and theta are 0). If depo vectors are not random there is no need to make unnecessary multiplications and/or risk floating point errors.
			candidate_xyz[0] = iatom_xyz[0];
			candidate_xyz[1] = iatom_xyz[1];
			candidate_xyz[2] = iatom_xyz[2] + diff_dist;
			
		}

		remap3DArrayInPeriodicBox( lmp , candidate_xyz );
		
	}
	
	const bool candidateDiffHasCollisions( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int *neighbors , int &neighbors_num , double *candidate_xyz , const int &diffused_type , double *iatom_xyz , const int &iatom_type ){

		/// Checks for collisions between the diffusion atom (i.e., atom moving to the vacancy) and existing atoms in the simulation.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] neighbors IDs of neighbors of the parent atom.
		/// @param[in] neighbors_num number of neighbors of the parent atom.
		/// @param[in] candidate_xyz coordinates (x,y, and z) of vacancy (point towards which the atom diffuses).
		/// @param[in] diffused_type atom type of diffused atom.
		/// @param[in] iatom_xyz coordinates of the parent atom.
		/// @param[in] iatom_type atom type of parent atom.
		/// @return true or false if the diffused candidate atom has or does not have collisions with existing atoms, respectively.
		/// @see PAPRECA::getDiffEventsFromAtom()
		
		LAMMPS_NS::tagint *id = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "id" );
		int *type = (int *) lammps_extract_atom( lmp , "type" );
		double **pos = ( double **)lammps_extract_atom( lmp , "x" );
		
		
		//This check has to be done in a separate call, since iatom is not in the neib list of iatom
		if( atomsCollide( lmp , papreca_config , iatom_xyz , iatom_type , candidate_xyz , diffused_type ) ){ return true; }
		
		//Check for collisions between candidate_xyz and neib atoms of iatom
		for( int i = 0; i < neighbors_num; ++i ){ //Check for collisions between parent atom (iatom) neibs and 
			int ineib = getMaskedNeibIndex( neighbors , i ); //get Masked index from neib list
			if( atomsCollide( lmp , papreca_config , pos[ineib] , type[ineib] , candidate_xyz , diffused_type ) ){
				return true;
			}
		}
		
		return false;
	}
	
	void getDiffEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &iatom , int *neighbors , int &neighbors_num , std::vector< Event* > &events_local , ATOM2BONDS_MAP &atomID2bonds ){
		
		/// Checks if the current atom is candidate to diffusion events (a.k.a. PAPRECA::PredefinedDiffusionHop). Detected deposition events are inserted in the events_local vector of PAPRECA::Event objects (storing all events of the current MPI process).
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] iatom local index of current atom.
		/// @param[in] neighbors IDs of neighbors of the current atom.
		/// @param[in] neighbors_num number of neighbors of the current atom.
		/// @param[in,out] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::loopAtomsAndIdentifyEvents()
		
		
		if( !papreca_config.predefinedCatalogHasDiffusionHopEvents( ) ){ return; }
		
		LAMMPS_NS::tagint *atom_ids = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "id" ); //extract atom ids
		double **atom_xyz = ( double **)lammps_extract_atom( lmp , "x" );//extract atom positions
		int *atom_types = ( int *)lammps_extract_atom( lmp , "type" );//extract atom types
		
		const LAMMPS_NS::tagint iatom_id = atom_ids[iatom];
		const LAMMPS_NS::tagint iatom_type = atom_types[iatom];
		double *iatom_xyz = atom_xyz[iatom];
		
		
		//Dictate if your current atom is a candidate for diffusion events
		PredefinedDiffusionHop *diff_template = papreca_config.getDiffusionHopFromAtomType( iatom_type );
		
		if( diff_template ){ //diff_template will be NULL if the parent_type (iatom type) is not linked to a diffusion event
		
				
			double candidate_xyz[3];
			getDiffPointCandidateCoords( lmp , papreca_config , iatom_xyz , candidate_xyz , diff_template );
			
			const int diffused_type = diff_template->getDiffusedAtomType( );
			
			if( diff_template->getCustomStyle( ) == "Fe_4PO4neib" && !feCandidateHas4PO4Neibs( papreca_config , diff_template , atom_ids , atom_types , neighbors , neighbors_num , atomID2bonds ) ){ return; }
			
			if( !candidateDiffHasCollisions( lmp , papreca_config , neighbors , neighbors_num , candidate_xyz , diffused_type , iatom_xyz , iatom_type ) ){
				
				const double rate = diff_template->getRate( );
				const int is_displacive = diff_template->isDisplacive( );
							
				Diffusion *diff = new Diffusion( rate , candidate_xyz , iatom_id , iatom_type , is_displacive , diffused_type , diff_template );
				events_local.push_back( diff );
						
			}
				
		
		}
		
		
	}
	
	//Deposition events
	const bool atomIsInDepoScanRange( PaprecaConfig &papreca_config , double *iatom_xyz , double &film_height ){
		
		/// Checks if the parent atom is within the acceptable scan range for deposition events. The acceptable scan ranges can be modified through the PAPRECA input file.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] iatom_xyz coordinates of candidate atom
		/// @param[in] film_height height in current PAPRECA step.
		/// @return true if the parent atom is within the acceptable scan range or false otherwise.
		/// @see PAPRECA::getDepoEventsFromAtom()
		
		// If getHeightDepoScan( ) returns the default value (i.e., -1), then immediately return true as the user has not selected a deposition scan range.
		return( ( ( papreca_config.getHeightDepoScan( ) == -1 ) || ( ( iatom_xyz[2] <= film_height + papreca_config.getHeightDepoScan( )  ) && ( iatom_xyz[2] >= film_height - papreca_config.getHeightDepoScan( ) ) ) ) ? true : false );
	}
	
	void getDepoPointCandidateCoords( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , double *iatom_xyz ,  double *candidate_xyz , PredefinedDeposition *depo_template ){
		
		/// Fills candidate_xyz array with the coordinates of the deposition candidate. The deposition candidate coordinates DO NOT coincide with the parent atom coordinates. Depending on the user input settings in PAPRECA input, the deposition candidate coordinates can either be directly above the parent atom or on the surface of the upper hemisphere centered at the parent atom coordinates.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] iatom_xyz coordinates of parent atom.
		/// @param[in,out] candidate_xyz coordinates of candidate point for the deposition event.
		/// @param[in] depo_template Deposition template (PAPRECA::PredefinedDeposition) as initialized by the user (in the PAPRECA input file).
		/// @note In the current version deposition candidates are always placed above the candidate atom. This means that you cannot use this version to place a molecule underneath the parent atom. Hence, you cannot get a film growing towards -z.
		/** @note The user is given the chance to produce a point at a random position on the surface of a sphere and at a distance equal to depo_offset. depo_offset is defined as the offset distance for the deposition event (using the create_Deposition command in the PAPRECA input file).
		//Here, we use the random number generator (on the parent event proc) to generate 2 random numbers. Those 2 random numbers are used to generate a random vector on a sphere with a radius of depo_offset around the parent atom.*/
		/// @see PAPRECA::getDepoEventsFromAtom()
		
		const double depo_offset = depo_template->getDepoOffset( );
		double *mol_center = depo_template->getCenter( );
		
		if( papreca_config.depoVecsAreRandom( ) ){
			
			const double rnum1 = papreca_config.getUniformRanNum( );
			const double rnum2 = papreca_config.getUniformRanNum( );
			
			double theta = 2.0 * M_PI * rnum1; //Gives a number between 0 and 2PI.
			double phi = 0.5 * M_PI * rnum2; //Gives a number between 0 and pi/2 (for phi between pi/2 and pi negative z coords would be produced, which is unwanted for the case of thin film growth (i.e., deposition on top of current film)
			
			//1) Account for the mol_center, 2) offset candidate_xyz on a sphere around
			candidate_xyz[0] = iatom_xyz[0] + mol_center[0] + depo_offset * sin( phi ) * cos( theta );
			candidate_xyz[1] = iatom_xyz[1] + mol_center[1] + depo_offset * sin( phi ) * sin( theta );
			candidate_xyz[2] = iatom_xyz[2] + mol_center[2] + depo_offset * cos( phi );
		
		}else{
			//For non-random deposition vectors the candidate sits directly above the parent atom. 
			candidate_xyz[0] = iatom_xyz[0] + mol_center[0]; //Same equations as random vectors case (if phi and theta are 0). If depo vectors are not random there is no need to make unnecessary multiplications and/or risk floating point errors.
			candidate_xyz[1] = iatom_xyz[1] + mol_center[1];
			candidate_xyz[2] = iatom_xyz[2] + mol_center[2] + depo_offset;
		}

		remap3DArrayInPeriodicBox( lmp , candidate_xyz ); //Remap is necessary because if the candidate center lies outside of the simulation box, then create_atoms (LAMMPS function) will NOT create an atom AND WILL ALSO NOT THROW A WARNING!
		
	}
	
	const bool depoCandidateIsBelowRejectionHeight( PaprecaConfig &papreca_config , double *candidate_xyz , const double &film_height ){
		
		/// Checks if the deposition candidate point is below the rejection height (as set in the PAPRECA input).
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] candidate_xyz array storing the coordinates of the candidate deposition point.
		/// @param[in] film_height film height at the current PAPRECA step.
		/// @return true or false if the candidate point is below or above the rejection height, respectively.
		/// @see PAPRECA::getDepoEventsFromAtom()
		
		//Immediately return yes if getHeightDepoReject returns the default value (i.e., -1) meaning that the user did not set a reject deposition height.
		return( ( ( papreca_config.getHeightDepoReject( ) == -1 ) || ( candidate_xyz[2] <= papreca_config.getHeightDepoReject( ) + film_height ) ) ? true : false );
		
	}
	
	void getMolCoords( LAMMPS_NS::LAMMPS *lmp , double **mol_xyz , double **mol_dx , const int &mol_natoms , double *candidate_center ){
		
		/// Fills mol_xyz with the coordinates of the candidate molecule for deposition.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in,out] mol_xyz Temporary array storing the candidate (for deposition) molecule coordinates. Container mol_xyz is only used to check for possible collisions between the inserted molecule and system molecules.
		/// @param[in] mol_dx array storing the x,y, and z distances of each molecule atom from the molecule center. See molecule.h and molecule.cpp files in the LAMMPS source directory for more information.
		/// @param[in] mol_natoms total number of molecule atoms.
		/// @param[in] candidate_center coordinates of the candidate deposition point.
		/// @see PAPRECA::getDepoEventsFromAtom()
		
		for( int i = 0; i < mol_natoms; ++i ){
			
			for( int j = 0; j < 3; ++j ){
				
				mol_xyz[i][j] = candidate_center[j] + mol_dx[i][j];
			}
			//Remap coords into periodic box. This helps with the calculation of distances between atoms and also avoids error with the create_atoms function (the function will not create new atoms if they are inserted outside of the periodic box).
			//Careful, remapping has to be done after all 3 assignments per mol_xyz[i], otherwise you are sending garbage into the remapper and you can get stuck in infinite loops
			remap3DArrayInPeriodicBox( lmp , mol_xyz[i] );
				
		}
		
	}
	
	void initMolCoordsArr( double ***mol_xyz , const int &mol_natoms ){
		
		/// Initializes the molecule coordinates array used to check for collisions between the deposited molecule and system atoms.
		/// @param[in,out] mol_xyz array storing the x,y, and z distances of each molecule atom from the molecule center. See molecule.h and molecule.cpp files in the LAMMPS source directory for more information.
		/// @param[in] mol_natoms total number of molecule atoms.
		/// @note Molecule coordinates arrays initialized with this function have to be manually deleted using the deleteMolcoordsArr() function.
		/// @see PAPRECA::getDepoEventsFromAtom(), PAPRECA::deleteMolCoordsArr()

		*mol_xyz = new double*[mol_natoms];
		
		for( int i = 0; i < mol_natoms; ++i ){
		
			(*mol_xyz)[i] = new double[3];
		
		}

	}
	
	void deleteMolCoordsArr( double **mol_xyz , const int &mol_natoms ){
		
		/// Deletes the molecule coordinates array used to check for collisions between the deposited molecule and system atoms.
		/// @param[in,out] mol_xyz array storing the x,y, and z distances of each molecule atom from the molecule center. See molecule.h and molecule.cpp files in the LAMMPS source directory for more information.
		/// @param[in] mol_natoms total number of molecule atoms.
		/// @see PAPRECA::getDepoEventsFromAtom(), PAPRECA::initMolCoordsArr()
		
		for( int i = 0; i < mol_natoms; ++i ){
		
			delete[ ] mol_xyz[i];
		
		}
		
		delete[ ] mol_xyz;

	}
	
	bool atomHasCollisionWithMolAtoms( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , double *atom_xyz , const int &atom_type , const int &mol_natoms , double **mol_xyz , int *mol_atomtype ){

		/// Checks if the parent event atom has collisions with any of the inserted molecule atoms.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] atom_xyz coordinates of the parent atom.
		/// @param[in] atom_type atom type of parent atom.
		/// @param[in] mol_natoms total number of atoms in the inserted molecule.
		/// @param[in] mol_xyz a mol_natomsx3 array storing the coordinates of all inserted molecule atoms.
		/// @param[in] mol_atomtype array storing the atom types of all inserted molecule atoms.
		/// @return true or false if atom has collisions with molecule atoms or not, respectively.
		/// @see PAPRECA::getDepoEventsFromAtom(), PAPRECA::atomsCollide(), PAPRECA::candidateDepoHasCollisions()

		for( int j = 0; j < mol_natoms; ++j ){
			
			if( atomsCollide( lmp , papreca_config , mol_xyz[j] , mol_atomtype[j] , atom_xyz , atom_type  ) ){
				return true;	
			}
			
		}
		
		return false;
		
	}
	
	
	bool candidateDepoHasCollisions( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const int &nprocs , PaprecaConfig &papreca_config , int *neighbors , int neighbors_num , double *candidate_center , double *iatom_xyz , const int &iatom_type , PredefinedDeposition *depo_template ){
		
		/// Checks if the inserted molecule atoms have collisions with 1) the parent atom, 2) all neighbors of the parent event atom.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] neighbors IDs of all neighbors of the parent event atom.
		/// @param[in] neighbors_num number of atoms.
		/// @param[in] candidate_center coordinates of the center of mass of the candidate inserted molecule.
		/// @param[in] iatom_xyz coordinates of parent atom.
		/// @param[in] iatom_type atom type of parent atom.
		/// @param[in] depo_template Deposition template (PAPRECA::PredefinedDeposition) as initialized by the user (in the PAPRECA input file).
		/// @return true or false if candidate deposition has collisions or not, respectively.
		/// @see PAPRECA::getDepoEventsFromAtom(), PAPRECA::atomHasCollisionWithMolAtoms(), PAPRECA::atomsCollide()
		/// @note This function assumes that any potential collision between the inserted molecule and existing atoms in the system can be detected using the parent atom neighbors. For very big molecules there is a chance that a molecule atom collides with existing atoms not included in the parent atom neighbor list. Hence, this function might require modifications in the future.
		
		LAMMPS_NS::tagint *id = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "id" );
		int *type = (int *) lammps_extract_atom( lmp , "type" );
		double **pos = ( double **)lammps_extract_atom( lmp , "x" );
		
		//get mol information from mol name
		double **mol_dx = depo_template->getCoords( );
		int *mol_atomtype = depo_template->getAtomTypes( );
		const int mol_natoms = depo_template->getAtomsNum( );
		
		//get candidate mol coordinates
		double **mol_xyz = NULL;
		initMolCoordsArr( &mol_xyz , mol_natoms );
		getMolCoords( lmp , mol_xyz , mol_dx , mol_natoms , candidate_center );

		//Checking for collisions between the current atom (iatom) and the mol atoms has to be done in a separate function call(because the iatom coordinates are not in the iatom neighbor list).
		if( atomHasCollisionWithMolAtoms( lmp , papreca_config , iatom_xyz , iatom_type , mol_natoms , mol_xyz , mol_atomtype ) ){
			deleteMolCoordsArr( mol_xyz , mol_natoms );
			return true;
			
		}
		
		//Check for collisions will all the neighbors of the parent atom!
		for( int i = 0; i < neighbors_num; ++i ){ 
			int ineib = getMaskedNeibIndex( neighbors , i ); //get Masked index from neib list
			if( atomHasCollisionWithMolAtoms( lmp , papreca_config , pos[ineib] , type[ineib] , mol_natoms , mol_xyz , mol_atomtype ) ){
				deleteMolCoordsArr( mol_xyz , mol_natoms );
				return true;
			}
			
		}
		
		//If there are no collisions remember to delete the MolCoords array. We use that array just to check for collisions. The LAMMPS object stores the mol template information
		deleteMolCoordsArr( mol_xyz , mol_natoms );
		return false;
		
	}
	
	

	void getDepoEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs , const int &iatom , int *neighbors , int &neighbors_num , double &film_height , std::vector< Event* > &events_local ){
		
		/// Checks if a system atom is parent to a deposition event (a.k.a. PAPRECA::PredefinedDeposition). Detected deposition events are inserted in the events_local vector of PAPRECA::Event objects (storing all events of the current MPI process).
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] iatom local index of current atom.
		/// @param[in] neighbors array containing the IDs of the neighbors to the current atom.
		/// @param[in] neighbors_num total number of neighbors
		/// @param[in] film_height film height at the current PAPRECA step.
		/// @param[in,out] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @see PAPRECA::loopAtomsAndIdentifyEvents()
		
		
		if( !papreca_config.predefinedCatalogHasDepositionEvents( ) ){ return; }
		
		LAMMPS_NS::tagint *atom_ids = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "id" ); //extract atom ids
		double **atom_xyz = ( double **)lammps_extract_atom( lmp , "x" );//extract atom positions
		int *atom_types = ( int *)lammps_extract_atom( lmp , "type" );//extract atom types
		
		const LAMMPS_NS::tagint iatom_id = atom_ids[iatom];
		const int iatom_type = atom_types[iatom];
		double *iatom_xyz = atom_xyz[iatom];
		
		PredefinedDeposition *depo_template = papreca_config.getDepositionFromParentAtomType( iatom_type );
		if( depo_template ){ //depo_template==NULL if the parent atom type (iatom_type) is not linked to a deposition
		
			if( atomIsInDepoScanRange( papreca_config , iatom_xyz , film_height ) ){ 
				
				double candidate_center[3];
				getDepoPointCandidateCoords( lmp , papreca_config , iatom_xyz , candidate_center , depo_template );
				
				if( depoCandidateIsBelowRejectionHeight( papreca_config , candidate_center , film_height ) ){ //reject depo candidates above a certain point
				
						if( depo_template->hasVariableStickingCoeff( ) || papreca_config.getSurfaceCoverageFile( ).isActive( ) ){ depo_template->incrementDepositionTries( ); }//No need to reset in the beginning. Variables are reset within the calcVariableStickingCoeff member function of PredefinedDeposition, immediately after the calculation of the sticking coefficient
						
						if( !candidateDepoHasCollisions( lmp , proc_id , nprocs , papreca_config , neighbors , neighbors_num , candidate_center , iatom_xyz , iatom_type , depo_template ) ){
							
							double rot_pos[3] = { 0.0 , 0.0 , 1.0 }; //In this version we don't rotate the molecule at all, so just define a rotation axis and set theta to zero!
							
							Deposition *depo = new Deposition( depo_template->getRate( ) , candidate_center , rot_pos , 0.0 , 0 , depo_template->getAdsorbateName( ) , depo_template );
							events_local.push_back( depo );
							if( depo_template->hasVariableStickingCoeff( ) || papreca_config.getSurfaceCoverageFile( ).isActive( ) ){ depo_template->incrementDepositionSites( ); }
							
						}
					
					
				}
				
			}
			
		}
		
	}

	
	//Bond-breaking and formation events
	const bool headAtomIsCatalyzed( PredefinedReaction *reaction_template , int *atom_types , int *neighbors , int &neighbors_num ){
		
		/// Scans neighbors of parent predefined reaction candidate to check if catalyzing types exist in the neighborhood. The catalyzing types are provided in the PAPRECA input file.
		/// @param[in] reaction_template Reaction template (PAPRECA::PredefinedReaction or PAPRECA::PredefinedBondForm) as initialized by the user (in the PAPRECA input file).
		/// @param[in] atom_types number of LAMMPS atom types.
		/// @param[in] neighbors array storing the neighbor IDs of the parent PAPRECA::PredefinedReaction candidate atom.
		/// @param[in] neighbors_num number of neighbors of the parent candidate atom.
		/// @return true/false depending on whether a catalyzing type exists or not in the neighborhood of the reaction parent candidate atom.
		/// @note: This function can be called by either the getBondBreakingEventsFromAtom() or the getBondFormEventsFromAtom(). This means that the reaction template can be either a PAPRECA::PredefinedReaction (parent class) or a PAPRECA::PredefinedBondForm (derived class) object. If this function is called with the derived class object there is still no need to cast it to PAPRECA::PredefinedBondForm since we only use the parent class (PAPRECA::PredefinedReaction) functions here.
		/// @see PAPRECA::getBondBreakingEventsFromAtom(), PAPRECA::getBondFormEventsFromAtom()
		
		if( reaction_template->getCatalyzingTypes( ).empty( ) ){ return true; } //Immediately return true if there are no catalyzing types to avoid unnecessary scans over neibs
		
		//Scan through all the neighbors of iatom and see if there is at least one catalyzing type
		for( int i = 0; i < neighbors_num; ++i ){
				
			int ineib = getMaskedNeibIndex( neighbors , i ); //get Masked index from neib list
			int neib_type = atom_types[ineib];
			if( elementIsInVector( reaction_template->getCatalyzingTypes( ) , neib_type ) ){
				return true;
			}
				
				
		}
		
		return false;
			
				
	}
	
	void getBondBreakingEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &iatom , int *neighbors , int &neighbors_num , std::vector<Event*> &events_local , ATOM2BONDS_MAP &atomID2bonds ){

		/// Scans all atoms in the current MPI process for PAPRECA::BondBreak (a.k.a. PAPRECA::PredefinedReaction) events. Discovered events are inserted in the events_local vector of PAPRECA::Event objects (storing all events of the current MPI process).
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] iatom local (LAMMPS, per MPI process) index of candidate atom.
		/// @param[in] neighbors array storing the IDs of LAMMPS neighbors of the candidate atom.
		/// @param[in] neighbors_num total number of neighbors of the candidate atom.
		/// @param[in,out] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::loopAtomsAndIdentifyEvents()
		
		
		if( !papreca_config.predefinedCatalogHasBondBreakEvents( ) ){ return; }
		
		LAMMPS_NS::tagint *atom_ids = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "id" ); //extract atom ids
		int *atom_types = (int *) lammps_extract_atom( lmp , "type" ); //extract atom types
		LAMMPS_NS::tagint *molids = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "molecule" ); //Extract molecule
		
		const LAMMPS_NS::tagint iatom_id = atom_ids[iatom];
		const LAMMPS_NS::tagint &molid = molids[iatom];
		
		BOND_VECTOR &bonds = atomID2bonds[ iatom_id ];
		
		for ( const auto &bond : bonds ){
			
			if ( bond.parentAtomIsHead( ) ){ //This means that our current atom is head atom to that specific bond. We do that to prevent breaking events from getting "defined" twice, which would lead to incorrect KMC total events frequency (and incorrect time advancement).
				
				const int bond_type = bond.getBondType( );
				PredefinedReaction *break_template = papreca_config.getReactionFromBondType( bond_type ); //The event_list.h function will return NULL if the bond does not participate in breaking events
				if( break_template ){ //So, if break_template==NULL you will not enter this part
					if( headAtomIsCatalyzed( break_template , atom_types , neighbors , neighbors_num ) ){
						
						const double rate = break_template->getRate( );
						BondBreak *bond_break = new BondBreak( rate , iatom_id , bond.getBondAtom( ) , bond.getBondType( ) , break_template );
						events_local.push_back( bond_break ); //Polymorphism allows pushing back of children of Event class.
						//But, we definitely need pointers to correctly manage the memory of the general events_local container, since Event children can have different sizes and this would create slicing issues
						
					}
					
				}
				
			}
			
			
		}
		

	}
	
	const bool atomsBelong2TheSameMol( const LAMMPS_NS::tagint &iatom_mol , const LAMMPS_NS::tagint &jneib_mol ){

		/// Uses LAMMPS molecule IDs to decide if two atoms pbelong to the same mol.
		/// @param[in] iatom_mol ID of the first atom.
		/// @param[in] jneib_mol ID of the second atom (which is typically one of the neighbors of iatom).
		/// @return true if atoms belong to the same molecule or false if the atoms do not belong to the same molecule.
		
		return ( iatom_mol == jneib_mol ? true : false );

	}
	
	
	const bool atomHasMaxBonds( PaprecaConfig &papreca_config , ATOM2BONDS_MAP &atomID2bonds , const LAMMPS_NS::tagint &atom_id , const int atom_type ){

		/// Determines the total number of bonds for an atom by retrieving its PAPRECA::Bond objects vector. Then compares it to the user defined species max bonds.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @param[in] atom_id ID of atom.
		/// @param[in] atom_type type of current atom.
		/// @return true if the size of PAPRECA::Bond objects vector is greater or equal to the user defined species max bonds.

		
		return( ( atomID2bonds[atom_id].size( ) >= papreca_config.getMaxBondsFromSpecies( atom_type ) ) ? true : false );
		
	}
	
	bool bondBetweenAtomsExists( ATOM2BONDS_MAP &atomID2bonds , const LAMMPS_NS::tagint &atom1_id , const LAMMPS_NS::tagint &atom2_id ){ //Just search on the bond list of one of the atoms to see if you can find the other one

		/// The PAPRECA::Bond objects vector is retrieved from the atomID2bonds container for two atoms to determine if the two atoms are already bonded.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @param[in] atom1_id ID of the first atom.
		/// @param[in] atom2_id ID of the second atom.
		/// @return true or false depending on whether there is an existing bond between the two atoms or not, respectively.
		
		BOND_VECTOR &bonds = atomID2bonds[ atom1_id ];
		for( const auto &bond : bonds ){
		
			const LAMMPS_NS::tagint bond_atom = bond.getBondAtom( );
			
			if ( bond_atom == atom2_id ){
			
				return true;
			
			}
		
		
		} 

		return false;

	}
	
	const bool atomCandidatesAreLone( const LAMMPS_NS::tagint atom1_id , const LAMMPS_NS::tagint atom2_id , ATOM2BONDS_MAP &atomID2bonds ){
		
		/// Checks if two atoms are lone (i.e., if they do not have any bonds). This is done by retrieving the PAPRECA:Bond objects vector for each atom from the atomID2bonds container.
		/// @param[in] atom1_id ID of the first atom.
		/// @param[in] atom2_id ID of the second atom.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @return true if both atoms are lone (not bonded to any other atom) or false if one or both atoms contain atoms (i.e., if the size of their PAPRECA::Bond objects vector is non-zero.
		

		
		return( ( atomID2bonds[atom1_id].size( ) == 0 && atomID2bonds[atom2_id].size( ) == 0 ) ? true : false );

	}
	
	
	const bool atomHasMaxBondTypes( PaprecaConfig &papreca_config , ATOM2BONDS_MAP &atomID2bonds , const LAMMPS_NS::tagint &atom_id , const int &atom_type , const int &bond_type ){
		
		/// Checks if the number of bonds of a specific bond type for the current atom is greater or equal than the user defined max bond types of species. This is done by retrieving the PAPRECA::Bond objects vector of the current atom from the atomID2bonds map.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @param[in] atom_id id of current atom.
		/// @param[in] atom_type of current atom.
		/// @param[in] bond_type to be checked
		/// @return true if the number of bonds of a specific bond type for the current atom is greater or equal than the user defined max bond types of species. False otherwise.
		
		const int bonds_max = papreca_config.getMaxBondTypesOfSpecies( atom_type , bond_type );
		
		if(  bonds_max != -1 ){ //If the mapping exists search. Otherwise, return false immediately
			int bonds_cur = 0;
			BOND_VECTOR &bonds = atomID2bonds[atom_id];
			for( const auto &bond : bonds ){
				if( bond.getBondType( ) == bond_type ){ ++bonds_cur; }
				if( bonds_cur >= bonds_max ){ return true; }
			}
			
		}
		
		return false; //This can happen if we've scanned all bonds and did not exceed the bond limit or if the MaxBondTypes mapping returned -1 (which means that no mapping was se int the PAPRECA input file.
	}
	
	void getBondFormEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &iatom , int *neighbors , int &neighbors_num , std::vector<Event*> &events_local , ATOM2BONDS_MAP &atomID2bonds ){

		/// Scans a half neighbor list and detect bond formation (a.k.a. PAPRECA::PredefinedBondForm). Detected events are inserted in events_local vector, which is a container of PAPRECA::Event object and stores the events detected on the current MPI process.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] iatom local index of current atom.
		/// @param[in] neighbors IDs of neighbors of the current atom.
		/// @param[in] neighbors_num number of neighbors of the current atom.
		/// @param[in,out] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::loopAtomsAndIdentifyEvents()
		
		if( !papreca_config.predefinedCatalogHasBondFormEvents( ) ){ return; }
		
		//Get Lammps pointers
		LAMMPS_NS::tagint *atom_ids = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "id" ); //extract atom ids
		double **atom_xyz = ( double **)lammps_extract_atom( lmp , "x" );//extract atom positions
		int *atom_types = ( int *)lammps_extract_atom( lmp , "type" );//extract atom typess
		LAMMPS_NS::tagint *atom_mol = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "molecule" ); //Extract molecule of specific atom
		
		//Define atoms properties based on iatom
		const LAMMPS_NS::tagint iatom_id = atom_ids[iatom];
		
		LAMMPS_NS::tagint iatom_mol;
		if( atom_mol != NULL ){ iatom_mol = atom_mol[iatom]; } //Ensure that atom_mol exists before obtaining specific imol value (non-molecular systems do not have molecule ids and this can cause segmentation faults!)
		
		double *iatom_xyz = atom_xyz[iatom];
		const int iatom_type = atom_types[iatom];
			
		
		
		for( int j = 0; j < neighbors_num; ++j ){ //Scan all neibs of iatom (parent atom) on the half list
		
			const int jneib = getMaskedNeibIndex( neighbors , j ); //Get masked index
			
			//Get neib properties
			const LAMMPS_NS::tagint jneib_id = atom_ids[jneib];
			const LAMMPS_NS::tagint jneib_mol = atom_mol[jneib];
			double *jneib_xyz = atom_xyz[jneib];
			const int jneib_type = atom_types[jneib];
			
		
			
			INT_PAIR type_pair( iatom_type , jneib_type );
			PredefinedBondForm *form_template = papreca_config.getBondFormFromAtomTypesPair( type_pair ); //The inverse pair is already in the catalog so no need to search twice
				
			if( form_template ){ //This indicates that the bond is formable (i.e., form_template is not NULL).
				
				if( atom_mol != NULL && !form_template->isSameMol( ) && atomsBelong2TheSameMol( iatom_mol , jneib_mol ) ){ continue; } //This is only for formation events from templates with same_mol=false;
				
				if ( !atomHasMaxBonds( papreca_config , atomID2bonds , iatom_id ,iatom_type ) && !atomHasMaxBonds( papreca_config , atomID2bonds , jneib_id , jneib_type ) && !bondBetweenAtomsExists( atomID2bonds , iatom_id , jneib_id ) ){ //Those conditions are now checked sequencially based on computational cost (i.e., as you move further in the nested loop in costs more to check if the condition is true).			
				
					const int bond_type = form_template->getBondType( );
					if( atomHasMaxBondTypes( papreca_config , atomID2bonds , iatom_id , iatom_type , bond_type ) || atomHasMaxBondTypes( papreca_config , atomID2bonds , jneib_id , jneib_type , bond_type ) ){ continue; }
					if( form_template->isLone( ) && !atomCandidatesAreLone( iatom_id , jneib_id , atomID2bonds ) ){ continue; } //This is only for formation events involving lone candidates.
							
					const double bond_sqr_dist = form_template->getBondDistSqr( );
													
					double pair_sqr_dist = get3DSqrDistWithPBC( lmp  , iatom_xyz , jneib_xyz ); //Run Custom minimum image distance function with periodicity along the x-, and y-directions (but now on the y-direction).
					if( pair_sqr_dist <= bond_sqr_dist ){ //If pair distance smaller than bonding distance add forming event to local events table
						const double rate = form_template->getRate( );
						BondForm *bond_form = new BondForm( rate , iatom_id , jneib_id , bond_type , form_template );
						events_local.push_back( bond_form );					
					}
				}
			}
		}
		
	}
	
	//Monoatomic Desorption events
	void getMonoDesEventsFromAtom( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , const int &iatom , std::vector< Event* > &events_local , ATOM2BONDS_MAP &atomID2bonds ){
		
		/// Checks if the parent atom is candidate to monoatomic desorption events (a.k.a. PAPRECA::PredefinedMonoatomicDesorption). Detected events are inserted in events_local vector, which is a container of PAPRECA::Event object and stores the events detected on the current MPI process.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] iatom local index of current atom.
		/// @param[in,out] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::loopAtomsAndIdentifyEvents()
		
		if( !papreca_config.predefinedCatalogHasMonoDesEvents( ) ){ return; }
		
		LAMMPS_NS::tagint *atom_ids = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "id" ); //extract atom ids
		int *atom_types = ( int *)lammps_extract_atom( lmp , "type" );//extract atom types
		
		const LAMMPS_NS::tagint iatom_id = atom_ids[iatom];
		BOND_VECTOR &bonds = atomID2bonds[iatom_id]; //retrieve unordered set from map, if map already exists
		
		const LAMMPS_NS::tagint iatom_type = atom_types[iatom];
		
		//Dictate if your current atom is a candidate for diffusion events
		PredefinedMonoatomicDesorption *monodes_template = papreca_config.getMonoatomicDesorptionFromAtomType( iatom_type );
		
		if( monodes_template ){
				
				if( bonds.empty( ) ){ //Currently, only LONE (non-bonded) single atom desorption is supported. This line checks if the current atom has bonds with other atoms.
					const double rate = monodes_template->getRate( );
					const int parent_type = monodes_template->getParentAtomType( );
					MonoatomicDesorption *monodes = new MonoatomicDesorption( rate , iatom_id , parent_type , monodes_template );
					events_local.push_back( monodes );
				}
			
		}
		
		
	}
	
	//General functions
	void  loopAtomsAndIdentifyEvents( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , int &nprocs , const int &KMC_loopid , PaprecaConfig &papreca_config , std::vector<Event*> &events_local , ATOM2BONDS_MAP &atomID2bonds , double &film_height ){
	
		/// 1) Calculates film height (if that is requested by the user). 2) Discovers events and inserts them in the PAPRECA::Event objects vector (storing all events detected on the current MPI process). PAPRECA::BondBreak, PAPRECA::Deposition, PAPRECA::MonoatomicDesorption, and PAPRECA::Diffusion are discovered through a full neighbors list. PAPRECA::BondForm are detected using a half-neighbors list. This is done because bondbreak, deposition, monoatomic desorption, and diffusion events perform collision checks which require full neighbor lists (because we need to know all the neighbors of each atom). On the other hand, it is more efficient to check for PAPRECA:BondForm events using a half neighbors list that include each neighbor pair once.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] KMC_loopid current PAPRECA step.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in,out] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @param[in,out] film_height film height at current PAPRECA step.
		/// @see PAPRECA::calcFilmHeight(), PAPRECA::getBondBreakingEventsFromAtom(), PAPRECA::getDepoEventsFromAtom(), PAPRECA::getDiffEventsFromAtom(), PAPRECA::getMonoDesEventsFromAtom(), PAPRECA::getBondFormEventsFromAtom()
		/// @note The user is advised to consult the LAMMPS documentation (https://docs.lammps.org/) for more information about neighbors lists as well as details related to the lammps_neighlist_num_elements andlammps_neighlist_element_neighbors functions used to retrieve the neighbor lists containers.
		
		calcFilmHeight( lmp , proc_id , KMC_loopid ,  papreca_config , film_height );
		
		int neiblist_id = lammps_find_fix_neighlist( lmp , "papreca" , 1 ); //Get neighbors list with ID 1 (full list as in the papreca fix)
		if( neiblist_id == -1 ){ allAbortWithMessage( MPI_COMM_WORLD , "Lammps could not find full neib list from fix papreca (fix papreca all papreca) command. Please ensure that the fix papreca command is present in your LAMMPS input file." ); }
		int iatom = -1 , neighbors_num = -1 , *neighbors = NULL;
		
		int atoms_num = lammps_neighlist_num_elements( lmp , neiblist_id ); //Find number of atoms in the "zero" neighbor list.
		//Loop over  full list
		for ( int i = 0; i < atoms_num; ++i ){
			//Get neibs list and iatom index
			lammps_neighlist_element_neighbors( lmp , neiblist_id , i , &iatom , &neighbors_num , &neighbors ); //get local atom index (iatom), number of neighbors of iatom, and indexes of iatom neighbors
			//Get events
			getBondBreakingEventsFromAtom( lmp , papreca_config , iatom , neighbors , neighbors_num , events_local , atomID2bonds );
			getDepoEventsFromAtom( lmp  , papreca_config , proc_id , nprocs , iatom , neighbors , neighbors_num  , film_height , events_local );
			getDiffEventsFromAtom( lmp , papreca_config , iatom , neighbors , neighbors_num , events_local , atomID2bonds );
			getMonoDesEventsFromAtom( lmp , papreca_config , iatom , events_local , atomID2bonds );
			
		}
		
		//Reset neib list variables
		neiblist_id = -1;
		iatom = -1;
		neighbors_num = -1;
		neighbors = NULL;
		
		//Now loop half neiblist atoms to get the bond formation events
		neiblist_id = lammps_find_fix_neighlist( lmp , "papreca" , 2 ); //Get neighbors list with ID 2 (half list as in the papreca fix)
		if( neiblist_id == -1 ){ allAbortWithMessage( MPI_COMM_WORLD , "Lammps could not find full neib list from fix papreca (fix papreca all papreca) command. Please ensure that the fix papreca command is present in your LAMMPS input file." ); }
		
		atoms_num = lammps_neighlist_num_elements( lmp , neiblist_id );
		for( int i = 0; i < atoms_num; ++i ){
			lammps_neighlist_element_neighbors( lmp , neiblist_id , i , &iatom , &neighbors_num , &neighbors ); //get local atom index (iatom), number of neighbors of iatom, and indexes of iatom neighbors
			getBondFormEventsFromAtom( lmp , papreca_config , iatom , neighbors , neighbors_num , events_local , atomID2bonds );
		}
		
	}

} //End of PAPRECA Namespace
