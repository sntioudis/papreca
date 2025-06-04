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

#include "event_execute.h"

namespace PAPRECA{
	
	//Formation events
	void fillFormTransferDataArr( BondForm *bond_form , int *form_data ){
		
		/// Serializes (prepares) data for calls to MPI function related to the execution of PAPRECA::BondForm events.
		/// @param[in] bond_form PAPRECA::BondForm event to be executed.
		/// @param[in,out] form_data 2-element array of int data.
		/// @note form_data[0] stores the bond type. form_data[0] tells us if we need to delete atoms after bond form or not (1 delete/0 don't delete).
		/// @see PAPRECA::deserializeFormTransferDataArr(), PAPRECA::executeBondForm()
		
		form_data[0] = bond_form->getBondType( );
		
		PredefinedBondForm *form_template = bond_form->getFormTemplate( );
		form_data[1] = form_template->isDeleteAtoms( );
	}
	
	void deserializeFormTransferDataArr( int *form_data , int &bond_type , int &delete_atoms ){
		
		/// Deserializes (post-processes) data after calls to MPI function related to the execution of PAPRECA::BondForm events.
		/// @param[in] form_data 2-element array of serialized (in fillFormTransferDataArr() function) data.
		/// @param[in,out] bond_type type of bond.
		/// @param[in,out] delete_atoms 0 if atoms are not deleted after the execution of the bond formation event, or 1 if atoms are deleted.
		/// @see PAPRECA::fillFormTransferDataArr(), PAPRECA::executeBondForm()
		
		bond_type = form_data[0];
		delete_atoms = form_data[1];
		
	}
	
	void executeBondForm( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &KMC_loopid , double &time , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event ){
		
		/// Executes PAPRECA::BondForm event. This is done by 1) communicating information from the MPI process that detected the event to all other MPI processes, and 2) calling formBond() and potentially deleteAtoms() after the appropriate data have been communicated.
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] event_proc MPI process calling the function (i.e., proc that detected this event).
		/// @param[in] selected_event PAPRECA::Event selected (to be executed). Casted to the correct type (i.e., PAPRECA::BondForm) before execution.
		/// @see PAPRECA::formBond(), PAPRECA::deleteAtoms()
		
		LAMMPS_NS::tagint atom_ids[2];
		atom_ids[0] = -1;
		atom_ids[1] = -2;
		int form_data[2];
		int bond_type = -3;
		int delete_atoms = -4;
		
		if( proc_id == event_proc ){ //retrieve information in event_proc !EXTRA CAUTION HERE: The selected_event pointer is NULL for all other procs EXCEPT for the event proc
			
			BondForm *bond_form = dynamic_cast<BondForm*>( selected_event ); //Cast as BondForm to access member variables of bond form
			
			atom_ids[0] = bond_form->getAtom1ID( );
			atom_ids[1] = bond_form->getAtom2ID( );
			fillFormTransferDataArr( bond_form , form_data );
			printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~EVENTS INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  \n Executing bond formation event from proc %d, BOND_TYPE=%d , ATOM1_ID=%d , ATOM2_ID=%d \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n" , proc_id , form_data[0] , atom_ids[0] , atom_ids[1] );
		}
		
		
		//Here we use 2 separate MPI_Bcasts for safety, since tagint adotps either sizeof( int_32 ) or sizeof( int_64 ) during compilation and casting everything as MPI_INT types might cause runtime errors or undefined behavior.
		MPI_Bcast( atom_ids , 2 , MPI_LMP_TAGINT , event_proc , MPI_COMM_WORLD ); //Using 2 different casts because there is a difference between LMP_TAGINT and MPI_INT types.
		MPI_Bcast( form_data , 2 , MPI_INT , event_proc , MPI_COMM_WORLD );
		deserializeFormTransferDataArr( form_data , bond_type , delete_atoms );

		formBond( lmp , atom_ids[0] , atom_ids[1] , bond_type ); //Now we can safely call this on all procs, since all procs know the important event details (i.e., atom1id, atom2id, bond_type )
		if( proc_id == 0 ){ papreca_config.getLogFile( ).appendBondForm( KMC_loopid , time , atom_ids[0] , atom_ids[1] , bond_type ); }
		
		if( delete_atoms ){ //This means that we have to deal with the bond formation even between 2 lone oxygens
		
			deleteAtoms( lmp , atom_ids , 2 , "no" , "no" );
		
		}
		
	}
	
	//Bond-breaking events
	void executeBondBreak( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &KMC_loopid , double &time , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event , ATOM2BONDS_MAP &atomID2bonds ){
		
		/// Executes PAPRECA::Bondbreak event. This is done by 1) communicating information from the MPI process that detected the event to all other MPI processes, and 2) calling deleteBond() after the appropriate data have been communicated.
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] event_proc MPI process calling the function (i.e., proc that detected this event).
		/// @param[in] selected_event selected (for execution) PAPRECA::Event. Casted to the correct type (i.e., PAPRECA::BondBreak) before execution.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::deleteBond()
		
		LAMMPS_NS::tagint atom_ids[2];
		atom_ids[0] = -1;
		atom_ids[1] = -2;//Initialize those values at something non realistic for all procs
		int bond_type = -3;
		
		
		if( proc_id == event_proc ){ //retrieve information in event_proc !EXTRA CAUTION HERE: The selected_event pointer is NULL for all other procs EXCEPT for event proc
			
			BondBreak *bond_break = dynamic_cast<BondBreak*>( selected_event ); //Cast as BondBreak to access member variables of bond break
			
			atom_ids[0] = bond_break->getAtom1ID( );
			atom_ids[1] = bond_break->getAtom2ID( );
			bond_type = bond_break->getBondType( );
			printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~EVENTS INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n Executing bond break event from proc %d, bond_type=%d , atom1_id = %d , atom2_id = %d \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n" , proc_id , bond_type , atom_ids[0] , atom_ids[1] );
		}
		
		MPI_Bcast( atom_ids , 2 , MPI_LMP_TAGINT , event_proc , MPI_COMM_WORLD ); //Proc selection was done on proc event_proc so now we need to communicate the event_proc to all procs
		MPI_Bcast( &bond_type , 1 , MPI_INT , event_proc , MPI_COMM_WORLD ); //Cast bond type separately to avoid data type runtime errors
		
		//Breakbond is part of lammps_wrappers
		deleteBond( lmp , atom_ids[0] , atom_ids[1] , 1 ); //Now we can safely call this on all procs, since all procs know the important event details (i.e., atom1id, atom2id ). Delete special if you are using fix_shake and/or you want to recompute the pairwise lists.
		
		//Configure internal nve/limit integrator if selected by user in the PAPRECA input file
		
		if( papreca_config.getNveLimSteps( ) != - 1 ){ //Meaning that the nve limit option was set in the input file
				papreca_config.insertEventAtomIDs2NveLimGroup( {atom_ids[0] , atom_ids[1]} ); //Initialize a TAGINT_VEC from the two communicated tagints for current breaking event
		}
		
		

		if( proc_id == 0 ){ papreca_config.getLogFile( ).appendBondBreak( KMC_loopid , time , atom_ids[0] , atom_ids[1] , bond_type ); }
		
	}


	void fillDepoDataTransfArr( double *depo_data , Deposition *depo ){
		
		/// Serializes (prepares) data for calls to MPI function related to the execution of PAPRECA::Deposition events.
		/// @param[in,out] depo_data 8-element array of double data for transfer.
		/// @param[in] depo PAPRECA::Deposition event to be executed.
		/// @note depo_data[0], depo_data[1], depo_data[2] are the deposition site coordinates (x,y, and z). depo_data[3], depo_data[4], depo_data[5] are the coordinates of the center of rotation of the molecule. depo_data[6] is the angle of rotation. depo_data[7] is the insertion velocity.
		/// @see PAPRECA::deserializeDepoTransfData(), PAPRECA::executeDeposition()
		
		for( int i = 0; i < 3; ++i ){
				
			depo_data[i] = depo->getSitePos( )[i];
				
		}
			
		for( int i = 3; i < 6; ++i ){
				
			depo_data[i] = depo->getRotPos( )[i-3];
				
		}
			
		depo_data[6] = depo->getRotTheta( );
		depo_data[7] = depo->getDepoTemplate( )->getInsertionVel( );

	}
	
	void deserializeDepoTransfData( double *depo_data , double *site_pos , double *rot_pos , double &rot_theta , double &insertion_vel ){
		
		/// Deserializes (post-processes) data after calls to MPI function related to the execution of PAPRECA::Deposition events.
		/// @param[in] depo_data 8-element array of serialized (in fillDepoDataTransfArr() function) data.
		/// @param[in,out] site_pos array containing the coordinates of the center-of-mass of the inserted molecule/atom.
		/// @param[in,out] rot_pos array containing the coordinates of the center-of-rotation of the inserted molecule.
		/// @param[in,out] rot_theta angle of rotation of inserted molecule
		/// @param[in,out] insertion_vel velocity of inserted molecule.
		/// @see PAPRECA::fillDepoDataTransfArr(), PAPRECA::executeDeposition()
		
		for( int i = 0; i < 3; ++i ){
				
			site_pos[i] = depo_data[i]; //Get the pos array from the first three values of depo_data
				
		}
			
		for( int i = 3; i < 6; ++i ){
				
			rot_pos[i-3] = depo_data[i]; //Get the rotation from positions 3-6
				
		}
		
		rot_theta = depo_data[6]; //Get the rotation angle theta from position 7 in the depo_data matrix
		insertion_vel = depo_data[7]; //Get the insertion velocity from the last position
		
		
	}
	
	void executeDeposition( LAMMPS_NS::LAMMPS *lmp , int &KMC_loopid , double &time , PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event ){
		
		/// Executes PAPRECA::Deposition event. This is done by 1) communicating information from the MPI process that detected the event to all other MPI processes, and 2) calling insertMolecule() after the appropriate data have been communicated.
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] event_proc MPI process calling the function (i.e., proc that detected this event).
		/// @param[in] selected_event PAPRECA::Event selected (to be executed). Casted to the correct type (i.e., PAPRECA::Deposition) before execution.
		/// @see PAPRECA::insertMolecule(), PAPRECA::resetMobileAtomsGroups()
		
		double depo_data[8];
		char mol_name[1024] = "NONE";
		int n_mol_name = -1;
		
		double site_pos[3];
		double rot_pos[3];
		double rot_theta;
		double insertion_vel;
		
		if( proc_id == event_proc ){ //Go to depo event proc
		
			Deposition *depo = dynamic_cast<Deposition*>( selected_event ); //Cast as deposition to access member variables of deposition
			fillDepoDataTransfArr( depo_data , depo );
			
			strcpy( mol_name , depo->getMolName( ).c_str( ) ); //Get adsorbate name
			n_mol_name = strlen( mol_name ) + 1; //Get string size of adsorbate name
			
			printf( " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~EVENTS INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n Executing deposition event from proc %d, MOL_NAME=%s center_pos=(%f,%f,%f) \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n" , proc_id , mol_name , depo_data[0] , depo_data[1] , depo_data[2] );
			
			
			
		}
		
		//Now cast the doubles
		MPI_Bcast( depo_data , 8 , MPI_DOUBLE , event_proc , MPI_COMM_WORLD );
		
		//Finally, BCast the char (first needs casting of n_mol_name to get char size)
		MPI_Bcast( &n_mol_name , 1 , MPI_INT , event_proc , MPI_COMM_WORLD );
		MPI_Bcast( mol_name , n_mol_name , MPI_CHAR , event_proc , MPI_COMM_WORLD );
		
		//Now deserialize the data on all procs (from the transferred vector) for use with the deposit function
		deserializeDepoTransfData( depo_data , site_pos , rot_pos , rot_theta , insertion_vel );
		
		//Now we are ready to call the insertMolecule function from the lammps_wrappers header, on all procs.
		insertMolecule( lmp , site_pos , rot_pos , rot_theta , 0 , mol_name );
		
		if( proc_id == 0 ){ papreca_config.getLogFile( ).appendDeposition( KMC_loopid , time , site_pos , rot_pos , rot_theta , insertion_vel , mol_name ); }

		if( insertion_vel != 0.0 ){
			lmp->input->one( "group new_mol subtract all fluid frozen" ); //This method is used to "select" all newly inserted atoms to a group. This allows to assign velocities. Since the newly inserted atom IS NOT IN ANY GROUP AT THIS MOMENT, you can subtract  fluid+frozen from all to get the new molecule atoms.
			std::string input_str = "velocity new_mol set NULL NULL " + std::to_string(insertion_vel) + " units box"; //NULL for x-y velocities means that we only the vertical velocities are set.
			lmp->input->one( input_str.c_str( ) );
			lmp->input->one( "group new_mol delete" ); //delete that group otherwise you will get a LAMMPS error in the next deposition (because of redefinition of new_mol group).
			resetMobileAtomsGroups( lmp , papreca_config ); //Reset atom groups immediately so you'll be ready for the next LAMMPS run.
		}
		
	}
	
	//Diffusion events
	void fillIntegerDiffDataTransfArray( int *diff_intdata , Diffusion *diff ){
		
		/// Serializes (prepares) integer data for calls to MPI function related to the execution of PAPRECA::Diffusion events.
		/// @param[in,out] diff_intdata 3-element array of integer data for transfer.
		/// @param[in] diff PAPRECA::Diffusion event to be executed.
		/// @note diff_intdata[0] contains the parent type. diff_intdata[1] is 0 if the diffusion type is displacive (see PAPRECA::Diffusion, and PAPRECA::PredefinedDiffusionHop definitions). diff_intdata[2] contains the diffused type.
		/// @see PAPRECA::deserializeIntegerDiffDataArr(), PAPRECA::executeDiffusion()
		
		diff_intdata[0] = diff->getParentType( );
		diff_intdata[1] = diff->isDisplacive( );
		diff_intdata[2] = diff->getDiffusedType( );
		
	}
	
	void fillDoubleDiffDataTransfArray( double *diff_doubledata , Diffusion *diff ){
		
		/// Serializes (prepares) double data for calls to MPI function related to the execution of PAPRECA::Diffusion events.
		/// @param[in,out] diff_doubledata 4-element array of double data for transfer.
		/// @param[in] diff PAPRECA::Diffusion event to be executed.
		/// @note diff_doubledata[0], diff_doubledata[1], diff_doubledata[2] contain the x,y, and z coordinates of the diffusion point (vacancy). diff_doubledata[3] contains the insertion velocity.
		/// @see PAPRECA::deserializeDoubleDiffDataArr(), PAPRECA::executeDiffusion()
		
		copyDoubleArray3D( diff_doubledata , diff->getVacancyPos( ) ); //We can use the copyDoubleArray3D (utilities.h/cpp) function here as it only copies the first 3 elements of the vector (default start/end arguments).
		diff_doubledata[3] = diff->getDiffTemplate( )->getInsertionVel( ); //The last double position contains the insertion velocity
		
	}
	
	void deserializeIntegerDiffDataArr( int *diff_intdata , int &parent_type , int &is_displacive , int &diffused_type ){
		
		/// Deserializes (post-processes) integer data after calls to MPI function related to the execution of PAPRECA::Diffusion events.
		/// @param[in] diff_intdata 3-element array of serialized (in fillIntegerDiffDataTransfArray() function) data.
		/// @param[in,out] parent_type atom type of parent atom.
		/// @param[in,out] is_displacive 0 or 1 depending on whether the diffusion event is diplacive or not.
		/// @param[in,out] diffused_type atom type of diffused atom.
		/// @see PAPRECA::fillIntegerDiffDataTransfArray(), PAPRECA::executeDiffusion()
		
		parent_type = diff_intdata[0];
		is_displacive = diff_intdata[1];
		diffused_type = diff_intdata[2];
		
	}
	
	void deserializeDoubleDiffDataArr( double *diff_doubledata , double *vac_pos , double &insertion_vel ){

		/// Deserializes (post-processes) double data after calls to MPI function related to the execution of PAPRECA::Diffusion events.
		/// @param[in] diff_doubledata 4-element array of serialized (in fillDoubleDiffDataTransfArray() function) data.
		/// @param[in,out] vac_pos coordinates of vacancy.
		/// @param[in,out] insertion_vel velocity of diffused atom.
		/// @see PAPRECA::fillDoubleDiffDataTransfArray(), PAPRECA::executeDiffusion()
		
		copyDoubleArray3D( vac_pos , diff_doubledata );
		insertion_vel = diff_doubledata[3];
	}

	void executeDiffusion( LAMMPS_NS::LAMMPS *lmp , int &KMC_loopid , double &time , PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event ){
		
		/// Executes PAPRECA::Diffusion event. This is done by 1) communicating information from the MPI process that detected the event to all other MPI processes, and 2) calling diffuseAtom() after the appropriate data have been communicated.
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] event_proc MPI process calling the function (i.e., proc that detected this event).
		/// @param[in] selected_event PAPRECA::Event selected (to be executed). Casted to the correct type (i.e., PAPRECA::Diffusion) before execution.
		/// @see PAPRECA::diffuseAtom(), PAPRECA::resetMobileAtomsGroups()
		
		double vac_pos[3] , insertion_vel;
		LAMMPS_NS::tagint parent_id = -1;
		int diff_intdata[3]; //This array stores the parent_type, is_displacive, and diffused_type in positions 0, 1, and 2 respectively.
		double diff_doubledata[4]; //This array stores the vac_pos in positions 1-3 and the insertion velocity in the last position
		
		int parent_type , is_displacive , diffused_type;
		
		if( proc_id == event_proc ){
			
			Diffusion *diff = dynamic_cast<Diffusion*>( selected_event ); //Cast as diffusion to access member variables of diffusion
			fillIntegerDiffDataTransfArray( diff_intdata , diff );
			fillDoubleDiffDataTransfArray( diff_doubledata , diff );
			parent_id = diff->getParentId( ); //This is a tagint, hence, it is communicated separately from the other ints
			
			printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~EVENTS INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n Executing diffusion event from proc %d, parent_id=%d , vac_pos=(%f,%f,%f) \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n" , proc_id , parent_id , diff_doubledata[0] , diff_doubledata[1] , diff_doubledata[2] );
		}
		
		//BCast VacPos (vacancy position) (positions 1-3) and insertion velocity (last position).
		MPI_Bcast( diff_doubledata , 4 , MPI_DOUBLE , event_proc , MPI_COMM_WORLD );
		//Bcast tagint (parent_id);
		MPI_Bcast( &parent_id , 1 , MPI_LMP_TAGINT , event_proc , MPI_COMM_WORLD );
		//BCast ints (parent_type , is_displacive , diffused_type );
		MPI_Bcast( diff_intdata , 3 , MPI_INT , event_proc , MPI_COMM_WORLD );
		
		//Deserialize the BCasted data
		deserializeIntegerDiffDataArr( diff_intdata , parent_type , is_displacive , diffused_type );
		deserializeDoubleDiffDataArr( diff_doubledata , vac_pos , insertion_vel );
		
		//Now safely call the relevant lammps_wrappers function
		diffuseAtom( lmp , vac_pos , parent_id , parent_type , is_displacive , diffused_type );
		if( proc_id == 0 ){ papreca_config.getLogFile( ).appendDiffusion( KMC_loopid , time , vac_pos , parent_id , parent_type , insertion_vel , is_displacive , diffused_type ); }
		
		if( insertion_vel != 0.0 ){
			lmp->input->one( "group new_atom subtract all fluid frozen" ); //Same as deposition insertion velocities. Probably an overkill to select a single atom using a subtract group. Can be made faster/better in future versions.
			std::string input_str = "velocity new_atom set NULL NULL " + std::to_string(insertion_vel) + " units box";
			lmp->input->one( input_str.c_str( ) );
			lmp->input->one( "group new_atom delete" );
			resetMobileAtomsGroups( lmp , papreca_config );
		}
		
	}
	
	//Monoatomic desorption events
	void executeMonoatomicDesorption( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , int &KMC_loopid , double &time , const int &proc_id , const int &nprocs , const int &event_proc , Event *selected_event ){
		
		/// Executes PAPRECA::MonoatomicDesorption event. This is done by 1) communicating information from the MPI process that detected the event to all other MPI processes, and 2) calling deleteAtoms() after the appropriate data have been communicated.
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] event_proc MPI process calling the function (i.e., proc that detected this event).
		/// @param[in] selected_event PAPRECA::Event selected (to be executed). Casted to the correct type (i.e., PAPRECA::MonoatomicDesorption) before execution.
		/// @see PAPRECA::deleteAtoms(), PAPRECA::resetMobileAtomsGroups()
		
		
		LAMMPS_NS::tagint atom_ids[1];
		atom_ids[0] = -1;
		int parent_type = -2;
		
		if( proc_id == event_proc ){ //retrieve information in event_proc !EXTRA CAUTION HERE: The selected_event pointer is NULL for all other procs EXCEPT for the event proc
		
			MonoatomicDesorption *monodes = dynamic_cast<MonoatomicDesorption*>( selected_event ); //Cast as BondForm to access member variables of bond form
			atom_ids[0] = monodes->getParentId( );
			parent_type = monodes->getParentType( );
			
			printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~EVENTS INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  \n Executing monoatomic desorption event from proc %d, PARENT_TYPE=%d , ATOM_ID=%d \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n" , proc_id , parent_type , atom_ids[0] );
		}
		
		MPI_Bcast( atom_ids , 1 , MPI_LMP_TAGINT , event_proc , MPI_COMM_WORLD );
		MPI_Bcast( &parent_type , 1 , MPI_INT , event_proc , MPI_COMM_WORLD );
		
		//In this version all detected MonoAtomicDesorption events involve lone (non-bonded) atoms. This means that we can now simply delete the atomID to perform a "desorption" event.
		deleteAtoms( lmp , atom_ids , 1 , "no" , "no" );
		if( proc_id == 0 ){ papreca_config.getLogFile( ).appendMonoatomicDesorption( KMC_loopid , time , atom_ids[0] , parent_type ); }
		
	}
	
	//General event execution functions
	
	void printStepInfo( PaprecaConfig &papreca_config , const int &KMC_loopid , const double &time , const double &film_height , const double &proc_rates_sum ){

		/// Prints essential PAPRECA information to the screen/terminal.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] KMC_loopid number of current PAPRECA step.
		/// @param[in] time current time.
		/// @param[in] film_height current height.
		/// @param[in] proc_rates_sum total (gathered/sum from all MPI processes) event rate.
		/// @note The printed information is NOT appended to the papreca.log file. papreca.log file stores information in a more simplified format (see export_files.h and export_files.cpp).
		
		printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~KMC INFO~~~~~~~~~~~~~~~~~~~~~~~~~~ \n" );
		printf( "This is KMC/MD step #%d \n" , KMC_loopid );
		printf( "The current time is %E seconds \n" , time );
		if( !papreca_config.getHeightMethod( ).empty( ) ){ printf( "The current Height is %f (Angstroms) \n" , film_height ); }
		printf( "The total rate on this step is %E hz \n" , proc_rates_sum );
		printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n" );
		
	}

	void executeEvent( LAMMPS_NS::LAMMPS *lmp , int &KMC_loopid , double &time , PaprecaConfig &papreca_config , const int &proc_id , const int &nprocs , const int &event_proc , const int &event_num , char *event_type , std::vector< Event* > &events_local , ATOM2BONDS_MAP &atomID2bonds ){
		
		/// Casts event type from parent MPI process to all other processes and then executes event. Currently, only PAPRECA::BondForm, PAPRECA::BondBreak, PAPRECA::Deposition, PAPRECA::Diffusion, PAPRECA::MonoatomicDesorption events are supported.
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] event_proc ID of proc that discovered the event.
		/// @param[in] event_num local (on the MPI process that discovered the event) event index (in the PAPRECA::Event objects vector).
		/// @param[in,out] event_type type of selected event.
		/// @param[in] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see lammps_wrappers.h and lammps_wrappers.cpp
		/// @note Information should be communicated among the processes before the execution of an event. This happens because to execute an event we need to call lammps_command or lmp->input->one. lammps_command and lmp->input->one functions have to be called simultaneously by all MPI processes (otherwise the program will hang).
		/// @note This function will almost certainly require a few changes if the user decides to modify the existing classes of events (bond form/break, deposition, diffusion) or add a new class of events. If changes are made the user will have to find a way to communicate data between procs (similar solutions as the solutions above can be used, of course).

		Event *selected_event = NULL; //Initialize this to NULL for all procs
		int n_event_type = -1;

		if( proc_id == event_proc ){ //Go to event proc
			selected_event = events_local[event_num]; //Point at the selected event on the event_proc
			strcpy( event_type , selected_event->getType( ).c_str( ) );//get event type on event_proc
			n_event_type = strlen( event_type ) + 1;
		}
		
		//Now the event_proc knows the type of event and needs to broadcast it to all other procs.
		MPI_Bcast( &n_event_type , 1 , MPI_INT , event_proc , MPI_COMM_WORLD ); //Cast the char length first (this is to enable the next MPI_Bcast of the actual string).
		MPI_Bcast( event_type , n_event_type , MPI_CHAR , event_proc , MPI_COMM_WORLD ); //Then cast the type of event.
		
		//Careful when you use the functions below. We initialized selected event as NULL. So, at this point ONLY the event_proc points to something that is not NULL!!!
		//We pass selected event in the function so careful!!
		
		if( !strcmp( event_type , "RXN-FORM" ) ){ //call proper function for execution depending on the event type. strcmp compares 2 strings and gives 0 if they are equal.
			executeBondForm( lmp , papreca_config , KMC_loopid , time , proc_id , nprocs , event_proc , selected_event );	
		}else if( !strcmp( event_type , "RXN-BREAK" ) ){
			executeBondBreak( lmp , papreca_config , KMC_loopid , time , proc_id , nprocs , event_proc , selected_event , atomID2bonds );
		}else if( !strcmp( event_type , "DEPO" ) ){
			executeDeposition( lmp , KMC_loopid , time , papreca_config , proc_id , nprocs , event_proc , selected_event );
		}else if( !strcmp( event_type , "DIFF" ) ){
			executeDiffusion( lmp , KMC_loopid , time , papreca_config , proc_id , nprocs , event_proc , selected_event );
		}else if( !strcmp( event_type , "MONO-DES" ) ){
			executeMonoatomicDesorption( lmp , papreca_config , KMC_loopid , time , proc_id , nprocs , event_proc , selected_event );
		}else{
			std::string unknown_type( event_type );
			allAbortWithMessage( MPI_COMM_WORLD , "Unknown event type " + unknown_type + " in executeEvent function in papreca.cpp." );
		}
		
	}
	
	int selectAndExecuteEvent( LAMMPS_NS::LAMMPS *lmp , int &KMC_loopid , double &time , char *event_type , int &proc_id , int &nprocs , PaprecaConfig &papreca_config , std::vector< Event* > &events_local , ATOM2BONDS_MAP &atomID2bonds , double &film_height ){
		
		/// Gets total event rate by gathering all local rates (i.e., sum of rates on a single MPI process). Then, selects an event MPI process using the N-FOLD way. Afterwards, an event is chosen from the selected MPI process and executed on all procs.
		/// @param[in,out] lmp pointer to LAMMPS object.
		/// @param[in] KMC_loopid current PAPRECA step.
		/// @param[in,out] time current time.
		/// @param[in,out] event_type type of selected event.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] nprocs total number of MPI processes.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @param[in] film_height current height.
		/// @note Currently, we discover events on all MPI processes. However, KMC events are not executed in parallel as only one event from one MPI processes is fired at a time. We plan to introduce parallel event execution in subsequent versions of PAPRECA. Additional code will have to be written to prevent errors in neighboring events (e.g., execution of 2 deposition events that overlap, breaking of the same bond twice, etc.). Of course, executing events in parallel is expected to elevate the scalability and boost the efficiency of the code even further.
		/// @note See this paper for more information regarding the classic N-FOLD way and the selection of events: https://www.sciencedirect.com/science/article/pii/S0927025623004159
		
		strcpy( event_type , "NONE" );//Starts with NONE and returned as NONE ONLY and ONLY if on event is selected (i.e., if the rate is zero). In any other case this variable will hold the event type.
		int zero_rate = 0; //Starts with zero and becomes 1 if the rate is zero. Then used to exit function prematurely and avoid segmentation faults.
		
		double rate_local = getLocalRate( events_local , papreca_config ); //That is the total rate of a specific proc
		double *proc_rates = new double[nprocs];
		MPI_Gather( &rate_local , 1 , MPI_DOUBLE , proc_rates , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); //Gather all values from all procs in master proc (0).
		
		//Get Event proc
		int event_proc = -1; //Initialize this to -1 to understand if something failed. Also, if all rates are 0, this will stay at -1 and you'll know that at the current step there are no detected events
		
		if( proc_id == 0 ){ //Select the event_proc on proc 0
			double proc_rates_sum = doubleArrSum( proc_rates , nprocs ); //Get cumulative rate from all procs to be used in stochastic event selection
			if( proc_rates_sum <= 0.0 ){
				//Signal all other procs that total rate is zero and exit function immediately to avoid segmentation faults! (i.e., due to division by zero in dt calculation and in the selectProcessStochastically function.
				zero_rate = 1;
			}else{
				advanceSimClockFromKMC( papreca_config , proc_rates_sum , time ); //Advance clock but only when the total rate is non-zero! Otherwise you will divide by 0 and your time will become inf.
				double rnum = papreca_config.getUniformRanNum( ); //Draw random number on proc 0
				event_proc = selectProcessStochastically( proc_rates , nprocs , rnum , proc_rates_sum );
			}
			
			printStepInfo( papreca_config , KMC_loopid , time , film_height , proc_rates_sum );
		}
		
		
		delete [] proc_rates;
		MPI_Bcast( &zero_rate , 1 , MPI_INT , 0 , MPI_COMM_WORLD ); //Now all procs know if the total rate is zero and skip the following statement to exit function prematurely
		if( zero_rate ){ return zero_rate; } //Immediately exit if the rate is zero on all procs. No need to go through event selection or execute an event in this case.
		
			
		MPI_Bcast( &event_proc , 1 , MPI_INT , 0 , MPI_COMM_WORLD ); //Proc selection was done on proc 0 so now we need to communicate the event_proc to all procs
			
		/*Draw rnum on driver proc and cast it to other procs to select event. Refrain from drawing random numbers on any other proc rather than the driver proc because this will lead to non repeatable results. Explanation: The same atoms might end up on different domains. This doesn't mean that the decomposition is different. It just means that the same decomposition might assign different proc numbers to the same domains. Now, the driver proc always selects random numbers to advance the simulation clock and select the event proc. This means that the random number sequence on the driver proc is on a different stage (i.e., will give you a different rnum) on the driver proc compared to all other procs. Because the proc numbering changes, we won't always get the same number of events executed on proc 0. Hence, we will probably get different random numbers and start having different results between different runs, even if these runs were initialized from the same random seed.*/
		double rnum;
		if( proc_id == 0 ){ rnum = papreca_config.getUniformRanNum( ); }
		MPI_Bcast( &rnum , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );

		//Get kMC process on Event Proc
		int event_num = -1;
		if( proc_id == event_proc ){ //Now go to the event proc to select a kMC event
			double *event_rates = new double[events_local.size( )];
			Event::fillRatesArr( event_rates , events_local );
			double event_rates_sum = doubleArrSum( event_rates , events_local.size( ) );
			event_num = selectProcessStochastically( event_rates , events_local.size( ) , rnum , event_rates_sum );
			delete [] event_rates;
		}
			
		MPI_Bcast( &event_num , 1 , MPI_INT , event_proc , MPI_COMM_WORLD ); //Proc selection was done on proc event_proc so now we need to communicate the event_proc to all procs
		//The above Bcast might be unnecessary, since we will only access event_num INSIDE event_proc and event_num does indeed have the current value inside the event_proc. Could possibly omit, but kept here for consistency.
		
		if( proc_id == 0 ){
			//Immediately abort with an error if no event proc could  be selected.
			if( event_proc == -1 ){ allAbortWithMessage( MPI_COMM_WORLD , "Could not select event proc in selectAndExecuteEvent function (papreca.cpp)." ); }
			
		}
		
		executeEvent( lmp , KMC_loopid , time , papreca_config , proc_id , nprocs , event_proc , event_num , event_type , events_local , atomID2bonds );
		
		//Because time is advanced on the master proc, the time value has to be BCasted to all other procs now, before exiting (if the rate is zero you don't have to bcast and the function will exit on the previous return
		MPI_Bcast( &time , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
		return zero_rate;//; //This tells you if the rate was zero on this step. Even though proc0 calculates the total rate, the information is communicated to all other procs (through MPI_Bcast).
						  //Hence, if the rate is zero (i.e., no detected events), then zero_rate would be 1 in all procs and returned properly to the main function.
		
		
	}
	
} //End of PAPRECA Namespace
