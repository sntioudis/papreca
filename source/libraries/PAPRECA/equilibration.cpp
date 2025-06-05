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
/// @brief Definitions for equilibration.h.

#include "equilibration.h"

namespace PAPRECA{
	
	//Delete Desorbed atoms
	void fillDelidsLocalVec( LAMMPS_NS::LAMMPS *lmp , const double &desorb_cut , std::vector< LAMMPS_NS::tagint > &delids_local , ATOM2BONDS_MAP &atomID2bonds ){
		
		/// Called by deleteDesorbedAtoms() and only when the delete_desorbed algorithm is set to gather_local. The function compares the z-coordinate of an atom. If the atom z-value is higher than desorb_cut, the atom ID is marked for deletion (i.e., inserted in the delids_local container) and it is later deleted along with its bonded atoms (retrieved from atomID2bonds map).
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] desorb_cut cutoff distance for atom deletion. Atoms whose z-coordinate is equal to or greater than desorb_cut are marked for deletion.
		/// @param[in,out] delids_local vector of collected atom IDs (delids_local will be different on each MPI process).
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::deleteDesorbedAtoms(), PAPRECA::Bond::recursiveCollectBondedAtoms()
		/// @note When deleting an atom, all atoms bonded with the deleted atom have to be deleted as well. Otherwise, LAMMPS will likely throw a "Missing Bond Atoms in proc %d" error. See lammps documentation (https://docs.lammps.org/) for more information.
		/// @note This function is not to be confused with fillDelidsVec(). fillDelidsVec() is only called when the user sets the delete_desorbed algorithm to gather_all, while this function (i.e., fillDelidsLocalVec() ) is called if the user set the delete_desorbed algorithm to gather_local.
		
		double **atom_xyz = ( double **)lammps_extract_atom( lmp , "x" );//extract atom positions
		const int natoms = *( ( int *)lammps_extract_global( lmp , "nlocal" ) );
		LAMMPS_NS::tagint *id = ( int *)lammps_extract_atom( lmp , "id" );
		TAGINT_SET delids_set; //This is to refrain from collecting the same id twice on the same proc (we use it in the recursive function that collects bonded atoms)
		
		for( int i = 0; i < natoms; ++i ){
			
			if( atom_xyz[i][2] >= desorb_cut ){
				
				if( !elementIsInUnorderedSet( delids_set , id[i] ) ){ //Avoid collecting atoms twice
				
					delids_set.insert( id[i] );
					delids_local.push_back( id[i] );
					PAPRECA::Bond::recursiveCollectBondedAtoms( id[i] , delids_local , delids_set , atomID2bonds ); //Collect all bonded atoms of inserted atom.
				}	
				
			}	
		}
		
		
	}
	
	bool delidsLocalVectorsAreEmpty( std::vector< LAMMPS_NS::tagint > &delids_local ){
		
		/// Communicates information between MPI processes to check if every local (i.e., on every MPI process) delids_local container is empty. Effectively, this means that the z-coordinates of all atoms in the simulation are lower than desorb_cut and that no deletions have to be performed.
		/// @param[in] delids_local vector of collected atom IDs.
		/// @return true if the delids_local vector is empty on all MPI processes. False, otherwise.
		/// @see PAPRECA::deleteDesorbedAtoms()
		
		const int num_delids_local = delids_local.size( );
		int num_delids_global;
		
		MPI_Allreduce( &num_delids_local , &num_delids_global , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );
		
		return num_delids_global == 0 ? true : false;
		
	}

	void gatherAndTrimDelIdsOnDriverProc( const int &proc_id , const int &nprocs , std::vector< LAMMPS_NS::tagint > &delids_local , std::vector< LAMMPS_NS::tagint > &delids_global ){
		
		/// Called by deleteDesorbedAtoms() and only when the delete_desorbed algorithm is set to gather_local. The gather_local algorithm compares the z-coordinates of atoms with desorb_cut locally (i.e., individually for each MPI process). After all comparisons are performed, this function (i.e., gatherAndTrimDelIdsOnDriverProc() ) communicates information among the MPI processes. Then, the function fills delids_global with the IDs of all atoms marked for deletion and makes sure that no duplicate deletions are performed.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs number of MPI processes.
		/// @param[in] delids_local contains the IDs of atoms marked for deletion on a specific MPI process (i.e., the delids_local container is different on each MPI process).
		/// @param[in,out] delids_global vector containing the IDs of all atoms marked for deletion (contains the same data on each MPI process as casted from the master process).
		/// @see PAPRECA::deleteDesorbedAtoms(), PAPRECA::fillDelidsLocalVec()
		
		/// Different MPI processes can return a different number of delids so we need to calculate the counts/displacements first and then use Gatherv
		std::vector< LAMMPS_NS::tagint > delids_gathered;
		int *recv_counts = new int[nprocs];
		int *displ = new int[nprocs];
		int num_gathered = 0;
		const int delids_local_size = delids_local.size( );
		
		
		MPI_Gather( &delids_local_size , 1 , MPI_INT , recv_counts, 1 , MPI_INT , 0 , MPI_COMM_WORLD ); //Get information regarding the sizes of the respective vectors (on master proc id 0 ).
		
		//Calculate displacements and recv_counts (needed for a call to Gatherv).
		if ( proc_id == 0 ){ 
			for( int i = 0; i < nprocs; ++i ){
				displ[i] = num_gathered;
				num_gathered += recv_counts[i];
			}
			delids_gathered.resize( num_gathered );
		}
		
		MPI_Gatherv( delids_local.data( ) , delids_local_size , MPI_INT , delids_gathered.data() , recv_counts , displ , MPI_INT , 0 , MPI_COMM_WORLD ); //Now deilds_gathered on proc 0 contains all the del ids.
		
		//we now have to trim the delids_gathered array to avoid duplicate ids .Gathered from 2 separate procs (e.g., atom on parent proc and same atom on ghost proc)
		int num_trim = 0;
		if( proc_id == 0 ){ 
			TAGINT_SET delids_trim_set; //We will use this to trim duplicates
			delids_trim_set.reserve( delids_gathered.size( ) );
			delids_global.reserve( delids_gathered.size( ) );
			
			for( const auto &id : delids_gathered ){
				
				if( !elementIsInUnorderedSet( delids_trim_set , id ) ){
					delids_trim_set.insert( id );
					delids_global.push_back( id );
				}
				
			}
			
			num_trim = delids_global.size( );
		}
		
		//Now proc 0 will bcast the number of trimed ids to allow the delids_global vector to resize correctly in all other procs
		MPI_Bcast( &num_trim , 1 , MPI_INT , 0 , MPI_COMM_WORLD ); 
		
		
		if( proc_id != 0 ){ delids_global.resize( num_trim ); } //resize delids global vector to number of trimmed to avoid segmentation faults or hangs when Bcasting
		
		//Finally, we Bcast the delids_global from proc 0 to all other procs
		MPI_Bcast( delids_global.data( ) , num_trim , MPI_INT , 0 , MPI_COMM_WORLD );


		delete [ ]recv_counts;
		delete [ ]displ;
		
		
	}

	int fillDelidsVec( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const double &desorb_cut , std::vector< LAMMPS_NS::tagint > &delids , ATOM2BONDS_MAP &atomID2bonds ){
		
		/// Called by deleteDesorbedAtoms() and only when the delete_desorbed algorithm is set to gather_all. The function compares the z-coordinate of an atom. If the atom z-coordinate is higher than desorb_cut, the atom ID is marked for deletion (i.e., inserted in the delids container) and it is deleted along with its bonded atoms (retrieved from atomID2bonds map). Here, a gather operation collects data from all atoms on the master MPI process (i.e., proc_id==0) before comparing the z-coordinates of atoms with desorb_cut.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] desorb_cut cutoff distance for atom deletion. Atoms whose z-coordinate is equal to or greater than desorb_cut are marked for deletion.
		/// @param[in,out] delids vector of collected atom IDs on the master MPI process.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::deleteDesorbedAtoms(), PAPRECA::broadcastDelidsFromMasterProc(), PAPRECA::Bond::recursiveCollectBondedAtoms()
		/// @note When deleting an atom, all atoms bonded with the deleted atom have to be deleted as well. Otherwise, LAMMPS will likely throw a "Missing Bond Atoms in proc %d" error. See lammps documentation (https://docs.lammps.org/) for more information.
		/// @note This function is not to be confused with fillDelidsLocalVec(). fillDelidsLocalVec() is only called when the user sets the delete_desorbed algorithm to gather_local, while this function (i.e., fillDelidsVec() ) is called if the user set the delete_desorbed algorithm to gather_all.
		
		int delids_num = 0;
		
		//Only the master proc performs this calculation and then broadcasts delids to other procs
		int natoms = *( int *)lammps_extract_global( lmp , "natoms" );
		LAMMPS_NS::tagint *atom_id = new LAMMPS_NS::tagint[natoms];
		int *atom_type = new int[natoms];
		double *atom_xyz = new double[3 * natoms];
			
		//Extract positions, ids, and types. Only proc 0 performs calculations. In later versions we should gather atoms ONLY on proc 0 (master proc).
		//gather concat is more expensive but necessary (since atom ids in our case might NOT be consecutive: see library.cpp and library.h from LAMMPS source files).
		lammps_gather_atoms_concat(lmp,(char *) "id" , 0 , 1 , atom_id );
		lammps_gather_atoms_concat(lmp,(char *) "type" , 0 , 1 , atom_type );
		lammps_gather_atoms_concat(lmp,(char *) "x" , 1 , 3 , atom_xyz );
		
		if( proc_id == 0 ){
		
			TAGINT_SET delids_set; //This is to refrain from collecting the same id twice on the same proc (we use it in the recursive function that collects bonded atoms)
			
			for( int i = 0; i < natoms; ++i ){
				
				if( atom_xyz[3 * i+2] >= desorb_cut ){ //3*(i+2) because gathered xyz data are stored in triplets (i.e., atom i: 3(i+0)->x, 3(i+1)->y, 3(i+2)->z
					
					if( !elementIsInUnorderedSet( delids_set , atom_id[i] ) ){ //Avoid collecting atoms twice
					
						delids_set.insert( atom_id[i] );
						delids.push_back( atom_id[i] );
						PAPRECA::Bond::recursiveCollectBondedAtoms( atom_id[i] , delids , delids_set , atomID2bonds ); //Collect all bonded atoms of inserted atom.
					}			
				}
			}
			delids_num = delids.size( );
		}
		
		//Delete malloc'ed arrays
		delete [ ] atom_id;
		delete [ ] atom_type;
		delete [ ] atom_xyz;
		
		MPI_Bcast( &delids_num, 1 , MPI_INT , 0 , MPI_COMM_WORLD ); //Broadcast delids_num so every proc knows that there are no delids to communicate
		
		return delids_num;
		
	}
	
	void broadcastDelidsFromMasterProc( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , int &delids_num , std::vector< LAMMPS_NS::tagint > &delids ){

		/// Broadcasts the delids vector (containing the IDs of atoms marked for deletion) from the master MPI process to all other processes. This function is only called when the delete_desorbed option is active and the gather_all algorithm is selected.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] delids_num number of atom IDs marked for deletion.
		/// @param[in,out] delids vector of collected atom IDs initialized by the master MPI process and broadcasted to all other MPI processes.
		/// @see PAPRECA::deleteDesorbedAtoms(), PAPRECA::fillDelidsVec(), PAPRECA::Bond::recursiveCollectBondedAtoms()
		
		
		//At this stage the delids vector on the master proc has the correct size (and ids). The delids vector on all other procs are empty.
		if( proc_id != 0 ){
		
			delids.resize( delids_num ); //For any other proc resize to delids_num (otherwise brodcast will fail).
		
		}

		//Now we can safely broadcast the delids data to all other procs
		MPI_Bcast( delids.data( ) , delids_num , MPI_INT , 0 , MPI_COMM_WORLD );
	}

	void deleteDesorbedAtoms( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &proc_id , const int &nprocs , double &film_height , ATOM2BONDS_MAP &atomID2bonds ){	
		
		/// This function is always called but performs computations only if the user has set a desorption height cutoff in the PAPRECA input file. The present function compares the z-coordinate of each atom with the desorption height cutoff. If any z-coordinate value is greater than or equal to the desorption height cutoff, the associated atom (along with its bonded atoms) is deleted. Currently, the user can select between two different algorithms: 1) gather_local (see fillDelidsLocalVec() function description/notes), and 2) gather_all (see fillDelidsVec() function description/notes). A comparison of the performance between the two algorithms is not currently available. However, as a quick note, it can be mentioned that gather_all is expected to be more memory intensive, since it calls the LAMMPS function lammps_gather_atoms_concat() to gather the coordinates of non-consecutive IDs on the master proc.
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] film_height height at current PAPRECA step.
		/// @param[in,out] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::fillDelidsLocalVec(), PAPRECA::fillDelidsVec()
		/// @note The user is advised to refer to the LAMMPS developer documentation (https://docs.lammps.org/) to understand how lammps_gather_atoms_concant works (called by gather_all deletion algorithm).
		
		if( papreca_config.getDesorptionHeight( ) == -1 ){ return; } //Immediately exit this function if the desorption height is not set by the user (i.e., if the desorption height in papreca config is equal to the default value (-1). No need to delete atoms in that case.
		double desorb_cut = film_height + papreca_config.getDesorptionHeight( ); //Units consistent with units in LAMMPS input. This is the distance above which we consider atoms to be desorbed (scaled by current film height).
			
		//Regardless of the trimming method Update bondslist and atomic positions
		atomID2bonds.clear( );
		runLammps( lmp , 0 ); //invoke run 0 to update neibors lists before gathering and deleting atoms
		//For molecular systems bond sort id is enabled by default. Hence, to use atomIDd2bonds maps we need to update our atomID2bonds maps. For non molecular systems IT MIGHT BE OK TO NOT UPDATE atomID2bonds but we do it anyway for safety.
		PAPRECA::Bond::initAtomID2BondsMap( lmp , proc_id , atomID2bonds );
		
		
		if( papreca_config.getDesorptionStyle( ) == "gather_local" ){ //Gather local means that we go through all atoms on all procs to find desorbed atoms. Then we gather all trim ids in the master proc and process (to avoid duplicate deletion ids).
			
			std::vector< LAMMPS_NS::tagint > delids_local;
			std::vector< LAMMPS_NS::tagint > delids_global;
			
			
			fillDelidsLocalVec( lmp , desorb_cut , delids_local , atomID2bonds );
			if( delidsLocalVectorsAreEmpty( delids_local ) ){
				return;	
			}else{
				
				gatherAndTrimDelIdsOnDriverProc( proc_id , nprocs , delids_local , delids_global );
				if( delids_global.size( ) <= papreca_config.getDesorbDelMax( ) ){ //Only delete atoms if the number of delids is smaller than the permitted (by the used) maximum number of atoms that can be deleted at once.
					deleteAtoms( lmp , delids_global , "no" , "no" );
					resetMobileAtomsGroups( lmp , papreca_config );
				}
			}
			
		}else if( papreca_config.getDesorptionStyle( ) == "gather_all" ){ //Gather all means we immediately gather all atoms in the master proc and process there. This option requires less inter-processor communication BUT probably necessitates more RAM.
			
			std::vector< LAMMPS_NS::tagint > delids;
			int delids_num = fillDelidsVec( lmp , proc_id , desorb_cut , delids , atomID2bonds );
			
			//Only perform those steps if there is at least one delid to delete
			if( delids_num != 0 ){
				broadcastDelidsFromMasterProc( lmp , proc_id , delids_num , delids );
				if( delids_num <= papreca_config.getDesorbDelMax( ) ){
					deleteAtoms( lmp , delids , "no" , "no" );
					resetMobileAtomsGroups( lmp , papreca_config );
				}
			}
		}else if( papreca_config.getDesorptionStyle( ) == "LAMMPS_region" ){
			
			if( desorb_cut < lmp->domain->boxhi[2] ){
				deleteAtomsInBoxRegion( lmp , lmp->domain->boxlo[0] , lmp->domain->boxhi[0] , lmp->domain->boxlo[1] , lmp->domain->boxhi[1] , desorb_cut , lmp->domain->boxhi[2] , "yes" , "no" );
				resetMobileAtomsGroups( lmp , papreca_config );
			}
		
		}else if( !papreca_config.getDesorptionStyle( ).empty( ) ){
			allAbortWithMessage( MPI_COMM_WORLD , "Desorbed atoms style is not an acceptable style (deleteDesorbedAtoms function in papreca.cpp)." );
		}

	}

	//Equilibration
	void equilibrateFluidAtoms( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , double &time , const unsigned long int &trajectory_duration ){
		
		/// Performs a LAMMPS simulation on the fluid atom types (as defined in the PAPRECA input). Then, updates the simulation clock by timestep*trajectory_duration (as defined by the user in the LAMMPS and PAPRECA inputs). Additionally, might perform minimizations before/after the LAMMPS trajectory (if an appropriate LAMMPS minimization command is defined by the user).
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in,out] time current time.
		/// @see PAPRECA::runLammps()
		resetMobileAtomsGroups( lmp , papreca_config ); //Reset mobile atom groups (i.e., add/remove atoms from the fluid group so you can be ready to run LAMMPS.
		
		//Minimization before traj
		if( !papreca_config.getMinimize1( ).empty( ) ){ lmp->input->one( papreca_config.getMinimize1( ).c_str( ) ); } //Only call the minimize functions IF a minimize LAMMPS command is defined! otherwise you will get a runtime error in LAMMPS
		
		//Set up nve limited groups if required
		if( papreca_config.nveLimGroupsAreActive( ) && !papreca_config.nveLimGroupIsEmpty( ) ){ setupNveLimIntegrator( lmp , papreca_config ); }
		
		//Run trajectory
		runLammps( lmp , trajectory_duration );
		
		if( papreca_config.nveLimGroupsAreActive( ) && !papreca_config.nveLimGroupIsEmpty( ) ){
			removeNveLimIntegrator( lmp , papreca_config );
			papreca_config.updateNveLimGroup( ); //Increments nve limited steps for limited atoms and removes atoms from relevant group if required
		}
		
		//Minimization (after trajectory)
		if( !papreca_config.getMinimize2( ).empty( ) ){ lmp->input->one( papreca_config.getMinimize2( ).c_str( ) ); }
		
		advanceSimClockFromLAMMPS( papreca_config , time );
		
	}
	
	void equilibrate( LAMMPS_NS::LAMMPS *lmp , int &proc_id , const int &nprocs , double &time , PaprecaConfig &papreca_config , double &film_height , int &zero_rate , const int &KMC_loopid , ATOM2BONDS_MAP &atomID2bonds ){

		/// This function performs a LAMMPS run on the current system configuration, every KMC_per_MD (as set by the user in the PAPRECA input file). Then, it deletes atoms whose z-coordinate is equal to or greater than the desorption height cutoff (defined in the PAPRECA input file).
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in,out] time current PAPRECA simulation time.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] film_height current PAPRECA simulation height.
		/// @param[in] zero_rate 0 if the total event rate at the current step is zero, or 1 otherwise.
		/// @param[in] KMC_loopid current PAPRECA simulation step.
		/// @param[in,out] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::equilibrateFluidAtoms(), PAPRECA::deleteDesorbedAtoms()
		/// @note The function also calculates the execution times during the LAMMPS (MD) step (if the executionTimes file has been activated in the PAPRECA input file).
		
			
		if( KMC_loopid % papreca_config.getKMCperLongMD( ) == 0 ){
		
			papreca_config.setMDTimeStamp4ExecTimeFile( KMC_loopid );
			equilibrateFluidAtoms( lmp , papreca_config , time , papreca_config.getLongTrajDuration( ) );
			papreca_config.calcMDTime4ExecTimeFile( nprocs , KMC_loopid );
			
			deleteDesorbedAtoms( lmp , papreca_config , proc_id , nprocs , film_height , atomID2bonds );
		
		}else if( KMC_loopid % papreca_config.getKMCperMD( ) == 0 || zero_rate ){
		
			papreca_config.setMDTimeStamp4ExecTimeFile( KMC_loopid );
			equilibrateFluidAtoms( lmp , papreca_config , time , papreca_config.getTrajDuration( ) );
			papreca_config.calcMDTime4ExecTimeFile( nprocs , KMC_loopid );
			
			deleteDesorbedAtoms( lmp , papreca_config , proc_id , nprocs , film_height , atomID2bonds );
		
		}
		
		

	}
	
} //End of PAPRECA Namespace
