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
/// @brief Definitions for LAMMPS wrapper functions. 

#include "lammps_wrappers.h"

namespace PAPRECA{

	///Initialize LAMMPS
	void initializeLMP( LAMMPS_NS::LAMMPS **lmp ){

		/// Instantiates a LAMMPS object
		/// @param[in,out] lmp pointer to LAMMPS object.
		
		*lmp = new LAMMPS_NS::LAMMPS( 0 , NULL , MPI_COMM_WORLD );
	
	}
	
	void readLMPinput( const std::string &lmp_input , LAMMPS_NS::LAMMPS *lmp ){
	
		/// Reads LAMMPS input and initializes LAMMPS parameters.
		/// @param[in] lmp_input name of LAMMPS input file (command-line argument).
		/// @param[in,out] lmp pointer to LAMMPS object.
		
		lmp->input->file( lmp_input.c_str( ) ); //LAMMPS constructor only uses c style strings so 
	}
	
	//Setup integrators
	void setupNveLimIntegrator( LAMMPS_NS::LAMMPS *lmp ,  PaprecaConfig &papreca_config ){
		
		/// Configures an nve/lim integrator for the Nve limited atoms (IDs included in papreca_config objects)
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] papreca_config configuration variable including basic simulation information. Used to retrieved the IDs of fluid atoms.
		/// @note This function uses this fix: https://docs.lammps.org/fix_nve_limit.html
		
		std::string input_str = "fix nve_limited_integration nve_limited nve/limit ";
		input_str += std::to_string( papreca_config.getNveLimDist( ) );
		lmp->input->one( input_str.c_str( ) );
			
		
	}
	
	void removeNveLimIntegrator( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config ){
		
		/// Removes nve/limit integrator previously set by setupNveLimIntegrator( ).
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] papreca_config configuration variable including basic simulation information. Used to retrieved the IDs of fluid atoms.
		
		lmp->input->one( "unfix nve_limited_integration");
		
	}
	
	
	//Execute LAMMPS
	void runLammps( LAMMPS_NS::LAMMPS *lmp , const int &timesteps_num ){
		
		/// Runs LAMMPS trajectory for given (previously initialized) LAMMPS object.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] timesteps_num number of timesteps to execute.
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/run.html
		
		if( timesteps_num < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Attempted to run trajectory with " + std::to_string( timesteps_num ) + " steps." ); }
		std::string input_str = "run " + std::to_string( timesteps_num ); 
		lmp->input->one( input_str.c_str( ) ); //convert c++ string to c++ string because LAMMPS can only operate with those
		
	}
	
	//Period Box Operations
	void remap3DArrayInPeriodicBox( LAMMPS_NS::LAMMPS *lmp , double *arr ){ 
	
		/// Receives an 1D array of 3 elements (i.e., x,y, and z coordinates) and remaps it inside the existing periodic box (as set in the provided LAMMPS object).
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] arr array of coordinates.
		/// @note See domain.h header for more information about the domain->remap function.
		
		lmp->domain->remap( arr );
	}	
	
	
	//kMC operations
	void resetMobileAtomsGroups4NveLimIntegration( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config ){
		
		/// Creates two mobile atom groups: one used for nve/limit integration and the "fluid" which contains all the mobile atoms minus those included in the nve/limit group
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] papreca_config configuration variable including basic simulation information. Used to retrieved the IDs of fluid atoms.
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/group.html.
		
		lmp->input->one( "group fluid clear" );
		lmp->input->one( "group nve_limited clear" );
		
		//Reset nve_limited group
		std::string input_str1 = "group nve_limited id ";
		input_str1 += papreca_config.exportNveLimIDs2String( );
		lmp->input->one( input_str1.c_str( ) ); //Now that the nve_limited group is defined, we need to subtract the limited atoms from the fluid group
		
		//Reset fluid group by subtracting nve_limited atoms from the fluid group
		std::string input_str2 = "group fluid_temp type ";
		for( auto &type : papreca_config.getFluidAtomTypes( ) ){ input_str2 += std::to_string( type ) + " "; }
		lmp->input->one( input_str2.c_str( ) );
		lmp->input->one( "group fluid subtract fluid_temp nve_limited");
		lmp->input->one( "group fluid_temp delete");
		
	}
	
	void resetMobileAtomsGroups( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config ){
		
		/// Clears the fluid group (containing the fluid atom types, as defined in the PAPRECA input file) and redefines it to. This ensures that all atoms of fluid types are included in the fluid group.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] papreca_config configuration variable including basic simulation information. Used to retrieved the IDs of fluid atoms.
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/group.html.
		/// @note This operation is necessary each time you add/remove atoms. This ensures that the correct group of atoms will keep moving in the simulation.
	
		lmp->input->one( "group fluid clear" );
		std::string input_str1 = "group fluid type ";
		for( auto &type : papreca_config.getFluidAtomTypes( ) ){ input_str1 += std::to_string( type ) + " "; }
		lmp->input->one( input_str1.c_str( ) );
		
		
	}
	
	void deleteAtoms( LAMMPS_NS::LAMMPS *lmp , LAMMPS_NS::tagint *atom_ids , const int &num_atoms , const std::string &delete_bonds , const std::string &delete_molecule ){

		/// Receives an array of atom IDs (atom_ids) and deletes all atoms in the simulation with the corresponding IDs.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] atom_ids array of IDs of atoms to be deleted.
		/// @param[in] num_atoms number of atoms to be deleted (length of atom_ids vector).
		/// @param[in] delete_bonds if "yes", then the LAMMPS delete_atoms command is called with the "bond yes" keyword.
		/// @param[in] delete_molecule if "yes", then the LAMMPS delete atoms command is called with the "mol yes" option.
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/delete_atoms.html.
		/// @note This function uses a dummy group named "deletion" to delete atoms. You should NOT use a group with the same name in your LAMMPS input file.
		
		if( delete_bonds != "yes" && delete_bonds != "no" ){ allAbortWithMessage( MPI_COMM_WORLD , "Unknown delete_bonds option: " + delete_bonds + " for deleteAtoms function in lammps_wrappers.cpp." ); }
		if( delete_molecule != "yes" && delete_molecule != "no" ){ allAbortWithMessage( MPI_COMM_WORLD , "Unknown delete_molecule option: " + delete_molecule + " for deleteAtoms function in lammps_wrappers.cpp." ); }
		
		//Insert all atoms to be deleted in the same group
		std::string input_str = "group deletion id ";
		for( int i = 0; i < num_atoms; ++i ){ input_str += std::to_string( atom_ids[i] ) + " "; }
		lmp->input->one( input_str.c_str( ) );
		
		if( ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "molecule" ) == NULL ){
			input_str = "delete_atoms group deletion"; //For non-molecular systems bond yes and mol yes options are not supported. Hence, they are discarded from the input command
		}else{
			input_str = "delete_atoms group deletion bond " + delete_bonds + " mol " + delete_molecule; //Delete desired atoms from deletion group and with selected bond and mol options
		}
		
		lmp->input->one( input_str.c_str( ) );
		
		input_str = "group deletion delete";
		lmp->input->one( input_str.c_str( ) );

	}
	
	void deleteAtoms( LAMMPS_NS::LAMMPS *lmp , std::vector< LAMMPS_NS::tagint > &atom_ids , const std::string &delete_bonds , const std::string &delete_molecule ){
		
		/// Receives a vector of atom IDs and deletes all atoms in the simulation with the corresponding IDs.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] atom_ids std::vector of IDs of atoms to be deleted.
		/// @param[in] delete_bonds if "yes", then the LAMMPS delete_atoms command is called with the "bond yes" keyword.
		/// @param[in] delete_molecule if "yes", then the LAMMPS delete atoms command is called with the "mol yes" option.
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/delete_atoms.html.
		/// @note This function uses a dummy group named "deletion" to delete atoms. You should NOT use a group with the same name in your LAMMPS input file.
		
		if( delete_bonds != "yes" && delete_bonds != "no" ){ allAbortWithMessage( MPI_COMM_WORLD , "Unknown delete_bonds option: " + delete_bonds + " for deleteAtoms function in lammps_wrappers.cpp." ); }
		if( delete_molecule != "yes" && delete_molecule != "no" ){ allAbortWithMessage( MPI_COMM_WORLD , "Unknown delete_molecule option: " + delete_molecule + " for deleteAtoms function in lammps_wrappers.cpp." ); }
		
		//Insert all atoms to be deleted in the same group
		std::string input_str = "group deletion id ";
		for( const auto &id :atom_ids ){
			
			input_str += std::to_string( id ) + " ";
			
		}
		lmp->input->one( input_str.c_str( ) );
		
		input_str = "delete_atoms group deletion bond " + delete_bonds + " mol " + delete_molecule; //Delete desired atoms from deletion group and with selected bond and mol options
		lmp->input->one( input_str.c_str( ) );
		
		input_str = "group deletion delete";
		lmp->input->one( input_str.c_str( ) );	
		
	}
	
	
	void deleteAtomsInBoxRegion( LAMMPS_NS::LAMMPS *lmp , double &boxxlo , double &boxxhi , double &boxylo , double &boxyhi , double &boxzlo , double &boxzhi , const std::string &delete_bonds , const std::string &delete_molecule ){
		
		/// Deletes all atoms in a block region defined by the input region bounds.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] boxxlo low deletion region x-boundary.
		/// @param[in] boxxhi high deletion region x-boundary.
		/// @param[in] boxylo low deletion region y-boundary.
		/// @param[in] boxyhi high deletion region y-boundary.
		/// @param[in] boxzlo low deletion region z-boundary.
		/// @param[in] boxzhi high deletion region z-boundary.
		/// @param[in] delete_bonds if "yes" then the delete_atoms LAMMPS command is called with the "bond yes" keyword. The delete_atoms command is called to delete atoms in the deletion region.
		/// @param[in] delete_molecule if "yes" then the delete_atoms LAMMPS command is called with the "mol yes" keyword. The delete_atoms command is called to delete atoms in the deletion region.
		/// @note This function is a wrapper of these LAMMPS commands: https://docs.lammps.org/region.html, https://docs.lammps.org/group.html, https://docs.lammps.org/delete_atoms.html
		/// @note This function uses a dummy group named "del_atoms" to delete atoms. You should NOT use a group with the same name in your LAMMPS input file.
		/// @note This function uses a dummy region named "del_region". You should NOT use a region with the same name in your LAMMPS input file.
		
		if( delete_bonds != "yes" && delete_bonds != "no" ){ allAbortWithMessage( MPI_COMM_WORLD , "Unknown delete_bonds option: " + delete_bonds + " for deleteAtomsInBoxRegion function in lammps_wrappers.cpp." ); }
		if( delete_molecule != "yes" && delete_molecule != "no" ){ allAbortWithMessage( MPI_COMM_WORLD , "Unknown delete_molecule option: " + delete_molecule + " for deleteAtomsInBoxRegion function in lammps_wrappers.cpp." ); }
		
		std::string input_str = "region del_region block " + std::to_string( boxxlo ) + " " + std::to_string( boxxhi ) + " " + std::to_string( boxylo ) + " " + std::to_string( boxyhi ) + " " + std::to_string( boxzlo ) + " " + std::to_string( boxzhi ) + " units box"; 
		lmp->input->one( input_str.c_str( ) ); //Create region with given box length
		
		lmp->input->one( "group del_atoms region del_region" ); //Define atoms to delete from region
		input_str = "delete_atoms group del_atoms bond " + delete_bonds + " mol " + delete_molecule + " ";
		lmp->input->one( input_str.c_str( ) ); //Delete atoms in the region
		
		lmp->input->one( "region del_region delete" ); //Delete region for later use of the "del_region" keyword.
	}
	
	void createAtom( LAMMPS_NS::LAMMPS *lmp , const double atom_pos[3] , const int &atom_type ){
		
		/// Inserts an atom in the simulation.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] atom_pos coordinates (x,y, and z) of insertion point.
		/// @param[in] atom_type type of insert atom
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/create_atoms.html
		/// @note CAREFUL: Currently, the inserted atom has zero velocity and charges!
		/// @note CAREFUL: LAMMPS will NOT create the atom if it lies exactly at or outside the periodic box! %PAPRECA uses other functions to remap coordinates inside the periodic box (e.g., PAPRECA::remap3DArrayInPeriodicBox).
		
		std::string input_str = "create_atoms " + std::to_string( atom_type ) + " single " + std::to_string( atom_pos[0] )+ " " + std::to_string( atom_pos[1] )+ " " + std::to_string( atom_pos[2] ) + " units box";
		lmp->input->one( input_str.c_str( ) );
		
		
	}
	
	void deleteBond( LAMMPS_NS::LAMMPS *lmp , const LAMMPS_NS::tagint &atom1id , const LAMMPS_NS::tagint &atom2id , const bool special ){
	
		/// Deletes an existing bond between atom1 and atom2.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] atom1id ID of the first atom.
		/// @param[in] atom2id ID of the second atom.
		/// @param[in] special if true, then the delete_bonds LAMMPS command is called with the "special" keyword.
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/delete_bonds.html
		/// @note The delete_bonds command in LAMMPS deletes ALL bonds of a specific type between 2 atoms
		/// Hence, to delete one specific bond, we temporarily move the bond to a dummy bond type group and then invoke the delete_bonds command
		/// Such method is safe since ONLY ONE BOND CAN BE FORMED BETWEEN A SET OF ATOMS (e.g., one P atom and one O atom should only have one bond between each other.
		/// If you want 2 specific atom IDS to form more that one bond between (e.g., a double/triple bond), consider using one bond type instead of 3 separate bond types.
		/// @note If an bond does not exist between the 2 atoms, LAMMPS will not throw an error! This happens because the bond deletion group will be empty by the time the delete_bonds command is called.
		/// @note This function uses a dummy group named bonddel to delete atoms. Do not define another group with the same name in your LAMMPS input file.
		
		std::string input_str = "group bonddel id " + std::to_string( atom1id ) + " " + std::to_string( atom2id ); //place atoms on the same bonddel group (via concentrating strings)
		lmp->input->one( input_str.c_str( ) );
		
		lmp->input->one( "set group bonddel bond 1" ); //Place all bonds between the two atoms in group 1 (dummy group)
		
		if( special ){ //The special command recomputes the 1-2, 1-3, 1-4 lists. It is NECCESARY  to use special yes when deploying a fix_shake and deleting bonds (you will probably get a segmentation fault if not used).
			input_str = "delete_bonds bonddel bond 1 remove special";//Remove all bonds of type 1 between atoms of atom1id and atom2id AND RECOMPUTE 1-2, 1-3, and 1-4 special lists.
		}else{
			input_str = "delete_bonds bonddel bond 1 remove";//Remove all bonds of type 1 between atoms of atom1id and atom2id (this is of course, one atom);
		}
		
		lmp->input->one( input_str.c_str( ) );
		lmp->input->one( "group bonddel delete" ); //Delete the bonddel group to avoid errors in later LAMMPS runs (i.e., doubly defined groups).
		
	}
	
	void formBond( LAMMPS_NS::LAMMPS *lmp , const LAMMPS_NS::tagint &atom1id , const LAMMPS_NS::tagint &atom2id , const int &bond_type ){
		
		/// Bonds atom1 and atom2 with a specific bond type (as defined in the LAMMPS input file).
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] atom1id ID of the first atom.
		/// @param[in] atom2id ID of the second atom.
		/// @param[in] bond_type type of bond to be formed.
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/create_bonds.html
		
		std::string input_str = "create_bonds single/bond " + std::to_string( bond_type ) + " " + std::to_string( atom1id ) + " " + std::to_string( atom2id );
		lmp->input->one( input_str.c_str( ) );

	}
	
	
	void insertMolecule( LAMMPS_NS::LAMMPS *lmp , const double site_pos[3] , const double rot_pos[3] , const double &rot_theta , const int &type_offset , const char *mol_name ){
		
		/// Inserts a molecule in the current LAMMPS system.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] site_pos array of insertion site coordinates. LAMMPS will insert a new molecule such that its center-of-mass coincides with site_pos.
		/// @param[in] rot_pos array of coordinates of the center of rotation.
		/// @param[in] rot_theta angle of molecule rotation.
		/// @param[in] type_offset add type_offset to the atom types of the molecule (as defined in the LAMMPS molecule template file). WARNING: For almost all cases you should use a value of 0 to ensure that the inserted atom types are identical to the atom types as defined in the LAMMPS molecule input file.
		/// @param[in] mol_name name of molecule (as defined in the LAMMPS input command. e.g., if you used this command: "molecule mmmTCP ./TCP.xyz", your mol_name should be "mmmTCP" ).
		/// @note This function is a wrapper of this LAMMPS command: https://docs.lammps.org/create_atoms.html
		/// @note the create_atoms command is called with the "units box keyword".
		
		//The first number you see in the command (0) corresponds to the offset in atom types. Please leave this to zero, because otherwise you are going to be adding some value to the actual atom type and will be messing data up.
		//The number you see below: 99999, is a random seed. We fix the random seed to a specific number because it essentially does not play any role.
		//Because, when you define all rot_theta and rot_pos the random seed is not used to orientate the molecule... The desired molecule orientation IS DEFINED BY THE TEMPLATE FILE!
		//In case you want to modify the orientation of the molecule on-the-fly (i.e., during the KMC run you can easily change the command line below.

		std::string input_str = "create_atoms " + std::to_string( type_offset ) + " single " + std::to_string( site_pos[0] ) + " " + std::to_string( site_pos[1] ) + " " + std::to_string( site_pos[2] ) + " mol " + mol_name + " 99999 rotate " + std::to_string( rot_theta ) + " " +
		std::to_string( rot_pos[0] ) + " " + std::to_string( rot_pos[1] ) + " " + std::to_string( rot_pos[2] ) + " units box"; //Units box are necessary because if you define a lattice your positions are going to be multiples of the lattice spacing

		lmp->input->one( input_str.c_str( ) );
		
	}
	
	void diffuseAtom( LAMMPS_NS::LAMMPS *lmp , const double vac_pos[3] , const LAMMPS_NS::tagint &parent_id , const int &parent_type , const int &is_displacive , const int &diffused_type ){
			
			/// Executes a diffusion operation (i.e., moves if the diffusion is displacive, or creates an atom if the diffusion is non-displacive).
			/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
			/// @param[in] vac_pos coordinates (x,y, and z) of vacancy.
			/// @param[in] parent_id ID of parent candidate atom.
			/// @param[in] parent_type atom type of parent candidate atom.
			/// @param[in] is_displacive 1 if the diffusion is displacive (i.e., if the parent atom moves), or 0 if the diffusion is non-displacive (i.e., if a new atom is created on the vacancy position).
			/// @param[in] diffused_type type of diffused atom. Can be the same as parent type or can be set to a different type if you wish to change the atom type after performing a diffusion event.
			/// @see createAtom(), deleteAtoms()
			/// @note If diffusion is displacive, the original atom moves. Instead of deleting bonds to move an atom it would be easier to delete the diffusing atom and spawn it again.
			/// This approximation introduces an error: it assumes that the charge of the atom becomes zero during diffusion.
			/// However, if you use any charge equilibration scheme this shouldn't be an issue because the atom is going to obtain the correct charge during equilibration.
			/// More sophisticated solutions would include: Actual moving of atom or even transferring (using set functions of charges).
			/// @note If diffusion is non-displacive, we create a new atom in the vacancy position instead of moving the parent atom.
			/// Again, we assume that the charge of the atom is 0 to begin with.
			
			if( is_displacive == 0 ){ //Now we simply create an atom at the vacancy pos
				createAtom( lmp , vac_pos , diffused_type );
				
			}else if( is_displacive == 1 ){//Displacive diffusion
				
				LAMMPS_NS::tagint *ids = new LAMMPS_NS::tagint[1];
				ids[0] = parent_id;
				
				deleteAtoms( lmp , ids , 1 , "yes" , "no" );
				createAtom( lmp , vac_pos , diffused_type );
				
				delete [ ] ids;
	
			}
			
	}
	
	//Sigmas
	void initType2SigmaFromLammpsPairCoeffs( LAMMPS_NS::LAMMPS *lmp , INTPAIR2DOUBLE_MAP &type2sigma ){
	
		/// Retrieves the pairstyle sigmas from LAMMPS to initialize the type sigmas in PAPRECA.
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in,out] type2sigma PAPRECA::INTPAIR2DOUBLE_MAP mapping a pair of atom types to their corresponding sigma.
		/// @note WARNING: You need to have a sigma pairstyle to use this function (i.e., sigma LJ or equivalent.).
		/// If you don't have such pair style you would have to MANUALLY add those parameters. The input file helps you setup/manage those options

		runLammps( lmp , 0 ); //Always run 0 before initializing sigmas from LAMMPS. This ensures that you retrieve the crossterms effectively (pair_modify mixes are computed on runtime)
		
		int types_num = *( int *)lammps_extract_global( lmp , const_cast<char*>( "ntypes" ) ); //Obtain number of types to allocate/retrieve sigma array (containing sigma pairstyle coeffs).
		int dim;
		double **sigma_temp = ( double **)lmp->force->pair->extract( const_cast<char*>( "sigma" ) , dim ); //retrieve sigma temp from LAMMPS
		
		//Certain sigma values have to be modified. There are 0 values inside this array (to avoid multiple storing).
		//However, CAREFUL because the values are not exactly equal to zero (due to rounding errors).
		//To check if something is close to zero we can obtain the numeric limits for doubles.
		for( int i = 1; i < types_num + 1; ++i ){
			
			for( int j =1; j < types_num + 1; ++j ){
				
				INT_PAIR pair( i , j ); //create int pair for pair types
				INT_PAIR pair_inverse( j , i );
				
				if( sigma_temp[i][j] < std::numeric_limits< double >::epsilon( ) ){ //COMPARE TO NUMERIC LIMITS EPSILON FOR DOUBLE AND NOT TO ZERO! OTHERWISE YOUR MAP WILL END UP HAVING 0 VALUES!
					type2sigma[pair] = sigma_temp[j][i]; //If the sigma_temp value is zero you can find the corresponding pair value using reverse i,j indexes
					type2sigma[pair_inverse] = sigma_temp[j][i];
				}else{
					type2sigma[pair] = sigma_temp[i][j];
					type2sigma[pair_inverse] = sigma_temp[i][j];
				}
				
			}
			
			
		}
	
	
	}
	
	//Neib Lists
	int getMaskedNeibIndex( int *neighbors , int &j ){
		
		/// Neighbor indexes in neighbor lists are masked with extra bits. This function unmasks the index and returns a valid local index to avoid segmentation faults.
		/// @param[in] neighbors array of IDs of neighbors.
		/// @param[in] j neighbor index.
		/// @return unmasked (valid) neighbor index.
		/// @note See here:https://matsci.org/t/requesting-a-full-neighbor-list-in-c-code-using-lammps-as-a-library/48294/4
		
		
		int jneib = neighbors[j]; 
		jneib &= NEIGHMASK;
		
		return jneib;
	}
	
	//Bond Lists
	void initAndGatherBondsList( LAMMPS_NS::LAMMPS *lmp , LAMMPS_NS::tagint **bonds_list , LAMMPS_NS::bigint &bonds_num ){

		/// Calls lammps_gather_bonds to fill array bonds_list. The container bonds_list contains the IDs of bonds of the LAMMPS object (lmp).
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in,out] bonds_list array containing the IDs of bonds.
		/// @param[in,out] bonds_num number of bonds in the bonds_list array.
		/// @see PAPRECA::Bond::fillAtomID2BondsContainer(), PAPRECA::Bond::addBond2BondVector(), PAPRECA::Bond::atomIDIsMapped(), PAPRECA::Bond::initAtomID2BondsMap()
		/// @note See LAMMPS documentation (https://docs.lammps.org/) for more information about lammps_gather_bonds.
		
		bonds_num = *( LAMMPS_NS::tagint *)lammps_extract_global( lmp , "nbonds" );
		*bonds_list = new LAMMPS_NS::tagint[ 3 * bonds_num ];
		lammps_gather_bonds( lmp , *bonds_list );

	}
	
	//Molecules
	const int getMolIndexFromMolName( LAMMPS_NS::LAMMPS *lmp , std::string mol_name ){
		
		/// Receives a molecule name and returns its molecule index. This allows one to access a LAMMPS molecule by using: lmp->atom->molecules[imol] (where imol is the molecule index).
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] mol_name valid molecule names (as defined in the LAMMPS input file. e.g., if you used this command: "molecule mmmTCP ./TCP.xyz", then mol_name should be "mmmTCP").
		/// @return valid molecule index. Aborts if index is -1 (i.e., if mol_name was not defined in the LAMMPS input file).
		/// @note This function wraps around the find_molecule LAMMPS function in atom.cpp header.
		
		const int index = lmp->atom->find_molecule( mol_name.c_str( ) );
		if( index == -1 ){
			allAbortWithMessage( MPI_COMM_WORLD , "Could not find mol_name " + mol_name + ". This typically happens due to a mismatch in the molecule names in LAMMPS input file and the PAPRECA input file (e.g., when defining a deposition event)." );
		}
		
		return index;
	}
	void computeMolCenter( LAMMPS_NS::LAMMPS *lmp , std::string mol_name ){
		
		
		/// Invokes the LAMMPS function compute_center for the mol_name.
		/// @param[in,out] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] mol_name valid molecule names (as defined in the LAMMPS input file. e.g., if you used this command: "molecule mmmTCP ./TCP.xyz", then mol_name should be "mmmTCP").
		/// @note This function is a wrapper around the compute_center of molecule.h.
		/// @note THis function can be helpful if you want to get the molecule center before the start of the run (i.e., to populate kMC events etc).
		/// @note This function HAS to be used to initialize the events list for deposition events. If not used you will get a segmentation fault as the molecule center will not be calculated.
		
		const int imol = getMolIndexFromMolName( lmp , mol_name );
		lmp->atom->molecules[imol]->compute_center( );
		
	
	}
	
	//Files
	void dumpRestart( LAMMPS_NS::LAMMPS *lmp , const int &KMC_loopid , const int &dump_freq ){
		
		/// Dumps a LAMMPS simulation restart file (can be used to restart a PAPRECA simulation).
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] KMC_loopid current PAPRECA step number.
		/// @param[in] dump_freq dump a restart file every dump_freq PAPRECA step.
		/// @note Restarts are controlled from this wrapper and NOT BY LAMMPS. LAMMPS will choose to create restarts based on the MD step, but we want to create restarts based on the KMC step.
		/// @note The dumping frequency can be set in the PAPRECA input file.
		
		if( KMC_loopid % dump_freq == 0 ){	
			lmp->input->one( "restart 1 ./papreca.restart" );
			lmp->input->one( "run 1" );
			lmp->input->one( "restart 0" );
			
			
		}
		
		
	}
	
	//Maths LAMMPS Wrappers
	double get3DSqrDistWithPBC( LAMMPS_NS::LAMMPS *lmp , const double *x1 , const double *x2 ){
		
		/// Calculates and returns the distances between 2 points, while accounting for any Periodic Boundary Conditions (PBC) in the system.
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] x1 array of coordinates of atom1.
		/// @param[in] x2 array of coordinates of atom2.
		/// @return minimum image distance between point1 (with coordinates as in array x1) and point2 (with coordinates as in array x2).
		double dx = x1[0] - x2[0];
		double dy = x1[1] - x2[1];
		double dz = x1[2] - x2[2];
		
		lmp->domain->minimum_image( "get3DSqrDistWithPBC func in lammps_wrappers.cpp of PAPRECA namespace" , 489 , dx , dy , dz ); //Use minimum_image to account for PBCs in the system
	
		return dx * dx + dy * dy + dz * dz;
		
		
	}
	

}//end of namespace PAPRECA
