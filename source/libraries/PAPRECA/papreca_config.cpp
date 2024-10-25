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
/// @brief Definitions for the PAPRECA::PaprecaConfig class.
#include "papreca_config.h"

namespace PAPRECA{
	
	//Constructors/Destructors
	PaprecaConfig::PaprecaConfig( ){ }
	PaprecaConfig::~PaprecaConfig( ){ delete rnum_gen; }
	
	//Essential parameters
	void PaprecaConfig::setKMCsteps( const unsigned long int &KMC_steps_in ){ KMC_steps = KMC_steps_in; }
	const unsigned long int &PaprecaConfig::getKMCsteps( ) const{ return KMC_steps; }
	void PaprecaConfig::setKMCperMD( const unsigned long int &KMC_per_MD_in ){ KMC_per_MD = KMC_per_MD_in; }
	const unsigned long int &PaprecaConfig::getKMCperMD( ) const{ return KMC_per_MD; }
	void PaprecaConfig::setTimeEnd( const double &time_end_in ){ time_end = time_end_in; }
	const double &PaprecaConfig::getTimeEnd( ){ return time_end; }
	
	//Random numbers
	void PaprecaConfig::initRanNumGenerator( LAMMPS_NS::LAMMPS *lmp , const int &seed ){ rnum_gen = new LAMMPS_NS::RanMars( lmp , seed ); }
	const double PaprecaConfig::getUniformRanNum( ) const{ return rnum_gen->uniform( ); } //Returns a pseudorandom number uniformly distributed between 0 and 1 (see RanMars.h header of LAMMPS).
	const bool PaprecaConfig::ranNumGeneratorIsInitialized( ) const{ return( rnum_gen == NULL ? false : true ); }
	
	//Atom groups
	void PaprecaConfig::setFluidAtomTypes( const std::vector< int > &fluid_atomtypes_in ){ fluid_atomtypes = fluid_atomtypes_in; } 
	const std::vector< int > &PaprecaConfig::getFluidAtomTypes( ) const{ return fluid_atomtypes; }
	void PaprecaConfig::setFrozenAtomTypes( const std::vector< int > &frozen_atomtypes_in ){ frozen_atomtypes = frozen_atomtypes_in; }
	const std::vector< int > &PaprecaConfig::getFrozenAtomTypes( ) const{ return frozen_atomtypes; }
	
	
	//PredefinedEvents
	PredefinedReaction *PaprecaConfig::getReactionFromBondType( const int &bond_type ){
		
		/// Returns a PAPRECA::PredefinedReaction object as defined in the PAPRECA:PredefinedEventsCatalog.
		/// param[in] bond_type type of bond.
		/// @return mapped PAPRECA::PredefinedReaction object.
		/// @note each bond can only be mapped to one PAPRECA::PredefinedReaction event.
		
		return( elementIsInUnorderedSet( predefined_catalog.bond_breaks_set , bond_type ) ? predefined_catalog.bond_breaks_map[bond_type] : NULL );
		
	}
	
	PredefinedBondForm *PaprecaConfig::getBondFormFromAtomTypesPair( const INT_PAIR &types_pair ){ 
	
		/// Returns a PAPRECA::PredefinedBondForm object as defined in the PAPRECA::PredefinedEventsCatalog.
		/// @param[in] types_pair std::pair( int , int ) of atom types.
		/// @return mapped PAPRECA::PredefinedBondForm object.
		/// @note each pair of types can only be mapped to one PAPRECA::PredefinedBondForm object.
		
		return( elementIsInUnorderedSet( predefined_catalog.bond_forms_set , types_pair ) ? predefined_catalog.bond_forms_map[types_pair] : NULL );
		
	}
	
	int PaprecaConfig::getMaxBondsFromSpecies( const int &atom_type ){ 
	
		/// Returns the maximum number of bonds of a specific atom type as defined in the PAPRECA::PredefinedEventsCatalog.
		/// @param[in] atom_type type of atom.
		/// @return maximum number of bonds for a specific atom type.
		
		return( mappingExists( predefined_catalog.bonds_max , atom_type ) ? predefined_catalog.bonds_max[atom_type] : std::numeric_limits<int>::max( ) ); 
		
	} //If no mapping is defined, return limits of int, then bonds will never be greater than limit of int (so no limit defined).
	
	int PaprecaConfig::getMaxBondTypesOfSpecies( const int &atom_type , const int &bond_type ){
		
		/// Returns the maximum number of bonds of a specific bond type, for a specific atom type.
		/// @param[in] atom_type type of atom.
		/// @param[in] bond_type type of bond to be checked.
		/// @return maximum number of bonds of a specific bond type for an atom type.
		
		if( mappingExists( predefined_catalog.bondtypes_max , atom_type ) ){
			INT2INT_MAP &species_maxbondtypes = predefined_catalog.bondtypes_max[atom_type];
			if( mappingExists( species_maxbondtypes , bond_type ) ){
				return species_maxbondtypes[bond_type];
			}
		}
		return -1; //In this case this means that the mapping does not exist (this is used in the atomHasMaxBondTypes function in main.cpp)
	}
	
	PredefinedDiffusionHop *PaprecaConfig::getDiffusionHopFromAtomType( const int &atom_type ){
		
		/// Returns a PAPRECA::PredefinedDiffusionHop object as defined in the PAPRECA::PredefinedEventsCatalog.
		/// @param[in] atom_type type of atom.
		/// @return mapped PAPRECA::PredefinedDiffusionHop object.
		/// @note each atom type can only be mapped to one PAPRECA::PredefinedDiffusionHop.
		
		return( elementIsInUnorderedSet( predefined_catalog.diffusions_set , atom_type ) ? predefined_catalog.diffusions_map[atom_type] : NULL );
		
	}
	
	PredefinedDeposition *PaprecaConfig::getDepositionFromParentAtomType( const int &atom_type ){
		
		/// Returns a PAPRECA::PredefinedDeposition object as defined in the PAPRECA::PredefinedEventsCatalog.
		/// @param[in] atom_type type of atom.
		/// @return mapped PAPRECA::PredefinedDeposition as defined in the PAPRECA::PredefinedEventsCatalog.
		/// @note each atom type can only be mapped to one PAPRECA::PredefinedDeposition.
		
		return( elementIsInUnorderedSet( predefined_catalog.depositions_set , atom_type ) ? predefined_catalog.depositions_map[atom_type] : NULL );
		
	}
	
	PredefinedMonoatomicDesorption *PaprecaConfig::getMonoatomicDesorptionFromAtomType( const int &atom_type ){
		
		/// Returns a PAPRECA::PredefinedMonoatomicDesorption object as defined in the PAPRECA::PredefinedEventsCatalog.
		/// @param[in] atom_type type of atom.
		/// @return mapped PAPRECA::PredefinedMonoatomicDesorption as defined in the PAPRECA::PredefinedEventsCatalog.
		/// @note each atom type can only be mapped to one PAPRECA::PredefinedMonoatomicDesorption.
		
		return( elementIsInUnorderedSet( predefined_catalog.monodes_set , atom_type ) ? predefined_catalog.monodes_map[atom_type] : NULL );
		
	}
	
	void PaprecaConfig::initPredefinedReaction( const int &atom1_type , const int &atom2_type , const int &bond_type , const double &rate , const std::vector< int > &catalyzing_types ){
		
		PredefinedReaction *reaction = NULL;
		
		if( catalyzing_types.size( ) == 0 ){
			reaction = new PredefinedReaction( atom1_type , atom2_type , bond_type , rate );
		}else{
			
			reaction = new PredefinedReaction( atom1_type , atom2_type , bond_type , rate , catalyzing_types );
		}

		predefined_catalog.bond_breaks_set.insert( bond_type );
		predefined_catalog.bond_breaks_map[bond_type] = reaction;
			
	}
	
	
	void PaprecaConfig::initPredefinedBondForm( const int &atom1_type , const int &atom2_type , const int &bond_type , const double &bond_dist , const int &delete_atoms , const int &lone_candidates , const bool &same_mol , const double &rate , const std::vector< int > &catalyzing_types ){
		
		PredefinedBondForm *bond_form = NULL;
	
		const double bond_dist_sqr = static_cast<double>( bond_dist * bond_dist );
		
		if( catalyzing_types.size( ) == 0 ){
			bond_form = new PredefinedBondForm( atom1_type , atom2_type , bond_type , rate , bond_dist_sqr , delete_atoms , lone_candidates , same_mol );
		}else{
			bond_form = new PredefinedBondForm( atom1_type , atom2_type , bond_type , rate , bond_dist_sqr , delete_atoms , lone_candidates , same_mol , catalyzing_types );
		}
		
		INT_PAIR type_pair( atom1_type , atom2_type );
		INT_PAIR type_pair_reverse( atom2_type , atom1_type );
		
		predefined_catalog.bond_forms_set.insert( type_pair );
		predefined_catalog.bond_forms_set.insert( type_pair_reverse );
		
		predefined_catalog.bond_forms_map[ type_pair ] = bond_form;
		predefined_catalog.bond_forms_map[ type_pair_reverse ] = bond_form;
		
	}
	
	void PaprecaConfig::initPredefinedDiffusionHop( const int &parent_type , const double &insertion_vel , const double &diff_dist , const bool &is_displacive , const int &diffused_type , const double &rate , const std::string &custom_style , const std::vector< int > &custom_atomtypes ){

		PredefinedDiffusionHop *diffusion = new PredefinedDiffusionHop( parent_type , insertion_vel , diff_dist,  rate , custom_style , custom_atomtypes , diffused_type , is_displacive );

		predefined_catalog.diffusions_set.insert( parent_type );
		predefined_catalog.diffusions_map[ parent_type ] = diffusion;
		
	}
	
	void PaprecaConfig::initPredefinedDeposition( LAMMPS_NS::LAMMPS *lmp , const int &parent_type , const double &depo_offset , const double &insertion_vel , const std::string &adsorbate_name , const double &rate , const bool &variable_sticking , const double &sticking_coeff ){
		
		PredefinedDeposition *depo = NULL;
		
		if( variable_sticking ){
			depo = new PredefinedDeposition( lmp , parent_type , rate , depo_offset , insertion_vel , adsorbate_name );
		}else{
			depo = new PredefinedDeposition( lmp , parent_type , rate , depo_offset , insertion_vel , adsorbate_name , sticking_coeff );
		}
		
		if( depo ){
			predefined_catalog.depositions_set.insert( parent_type );
			predefined_catalog.depositions_map[ parent_type ] = depo;
		}else{
			allAbortWithMessage( MPI_COMM_WORLD , "Improper deposition event initialization in papreca_config.cpp." );
		}
		
	}
	
	void PaprecaConfig::initPredefinedMonoatomicDesorption( const int &parent_type , const double &rate ){
		PredefinedMonoatomicDesorption *monodes = new PredefinedMonoatomicDesorption( parent_type , rate );
		
		predefined_catalog.monodes_set.insert( parent_type );
		predefined_catalog.monodes_map[parent_type] = monodes;
		
	}
	
	void PaprecaConfig::setSpeciesMaxBonds( const int &species , const int &bonds_max ){ predefined_catalog.bonds_max[ species ] = bonds_max; }
	void PaprecaConfig::setSpeciesMaxBondTypes( const int &species , const int &bond_type , const int &bonds_max ){ predefined_catalog.bondtypes_max[species][bond_type] = bonds_max; }
	
	void PaprecaConfig::calcStickingCoeffs( ){
	
		std::unordered_set< std::string > adsorbates;
		std::unordered_map< std::string , int > deposition_sites;
		std::unordered_map< std::string , int > deposition_tries;
		std::unordered_map< std::string , double > sticking_coeff;
		
		
		//Loop through all possible deposition events.
		for( auto it = predefined_catalog.depositions_map.begin( ); it != predefined_catalog.depositions_map.end( ); ++it ){
			
			PredefinedDeposition *depo = it->second;
			std::string adsorbate_name = depo->getAdsorbateName( );
			
			
			if( depo->hasVariableStickingCoeff( ) ){
				if( !elementIsInUnorderedSet( adsorbates , adsorbate_name ) ){
					adsorbates.insert( adsorbate_name );
					deposition_sites[adsorbate_name] = 0;
					deposition_tries[adsorbate_name] = 0;
					sticking_coeff[adsorbate_name] = -1.0; //We initialize this to 0, so there is no confusion between sticking=0 (i.e., when tries are 0) and "freshly-initialized" sticking coeffs.
					
				}
				
				int depdata_global[2] = { 0 , 0 };
				int depdata_local[2] = { depo->getDepositionSites( ) , depo->getDepositionTries( ) };
				MPI_Allreduce( depdata_local , depdata_global , 2 , MPI_INT , MPI_SUM , MPI_COMM_WORLD ); //Different processors have different events. Hence, they should have different deptries and depsites.
				
				//At this point, all procs have the same deptries/depsites FOR THAT SPECIFIC PREDEFINED DEPOSITION EVENT
				//However, the same adsorbate can have different sites and parent atoms. Hence, we collect all depsites and deptries corresponding to the same adsorbate name IN THE SAME CONTAINER
				deposition_sites[adsorbate_name] += depdata_global[0];
				deposition_tries[adsorbate_name] += depdata_global[1];
				
			}
			
		}
		
		//We now have to loop again through all deposition events and calculate the sticking coefficients for specific adsorbates.
		//No need to Allreduce data again here, since all processors should have the same deposition tries and deposition sites FOR THE SAME PREDEFINED EVENT.
		
		for( auto it = predefined_catalog.depositions_map.begin( ); it != predefined_catalog.depositions_map.end( ); ++it ){
			
			PredefinedDeposition *depo = it->second;
			std::string adsorbate_name = depo->getAdsorbateName( );
			
			
			if( depo->hasVariableStickingCoeff( ) ){
			
				if( sticking_coeff[adsorbate_name] == -1.0 ){ //This means that we haven't calculated the sticking coefficient for that specific adsorbate name yet.
					//Hence we calculate it here...
					if( deposition_tries[adsorbate_name] == 0 ){ //This is to prevent segmentation faults when no deposition tries are attempted at all
						sticking_coeff[adsorbate_name] = 0;
						
					}else{
						sticking_coeff[adsorbate_name] = static_cast< double >( deposition_sites[adsorbate_name] ) / static_cast< double >( deposition_tries[adsorbate_name] ); //Calculate the sticking coeff which is the number of free sites divided by the number of total sites searched
					}
					
				}
				
				//Regardless, we now set the sticking coefficient properly for that event.
				depo->setStickingCoeff( sticking_coeff[adsorbate_name] );
				
			}//No need to do anything if the deposition is not variable. If the deposition event has a constant sticking coeff, the relevant sticking coeff is already set properly.
		
			//After calculating the relevant rates make sure you reset the deposition tries and sites to prepare for the next run.
			depo->resetDepositionTriesAndSites( ); //Reset sites regardless of whether the deopsition event has a variable sticking or not (because you might use depsites/tries to calculate the surface coverage).
			//IMPORTANT: The calcStickingCoeff function has to be called AFTER the surface coverage function, otherwise the sites will be zeroed when you attempt to calculate the surface coverage.
		
		}
		
	}
	
	const bool PaprecaConfig::predefinedCatalogHasBondBreakEvents( ) const{ return( !predefined_catalog.bond_breaks_map.empty( ) ? true : false ); }
	const bool PaprecaConfig::predefinedCatalogHasBondFormEvents( ) const{ return( !predefined_catalog.bond_forms_map.empty( ) ? true : false ); }
	const bool PaprecaConfig::predefinedCatalogHasDiffusionHopEvents( ) const{ return( !predefined_catalog.diffusions_map.empty( ) ? true :  false ); }
	const bool PaprecaConfig::predefinedCatalogHasDepositionEvents( ) const{ return( !predefined_catalog.depositions_map.empty( ) ? true : false ); }
	const bool PaprecaConfig::predefinedCatalogHasMonoDesEvents( ) const{ return( !predefined_catalog.monodes_map.empty( ) ? true : false ); }		
	const bool PaprecaConfig::predefinedCatalogIsEmpty( ) const{
		
		return( ( !predefinedCatalogHasBondBreakEvents( ) && !predefinedCatalogHasBondFormEvents( ) && !predefinedCatalogHasDiffusionHopEvents( ) && !predefinedCatalogHasDepositionEvents( ) && !predefinedCatalogHasMonoDesEvents( ) ) ? true : false );
		
	}
	
	
	//Randon deposition/diffusion vectors
	void PaprecaConfig::setRandomDepoVecs( const bool &random_depovecs_in ){ random_depovecs = random_depovecs_in; }
	const bool &PaprecaConfig::depoVecsAreRandom( ) const{ return random_depovecs; }
	void PaprecaConfig::setRandomDiffVecs( const bool &random_diffvecs_in ){ random_diffvecs = random_diffvecs_in; }
	const bool &PaprecaConfig::diffVecsAreRandom( ) const{ return random_diffvecs; }
	void PaprecaConfig::setRandomDiffVecsStyle( const std::string &diffvecs_style_in ){ diffvecs_style = diffvecs_style_in; }
	const std::string &PaprecaConfig::getRandomDiffVecsStyle( ) const{ return diffvecs_style; }
	void PaprecaConfig::setDepoHeights( const double &height_deposcan_in , const double &height_deporeject_in ){
		
		height_deposcan = height_deposcan_in;
		height_deporeject = height_deporeject_in;
		
		
	}
	
	const double &PaprecaConfig::getHeightDepoScan( ) const{ return height_deposcan; }
	const double &PaprecaConfig::getHeightDepoReject( ) const{ return height_deporeject; }

	
	//Desorption settings
	void PaprecaConfig::setDesorptionHeight( const double &desorb_cut_in ){ desorb_cut = desorb_cut_in; }
	const double &PaprecaConfig::getDesorptionHeight( ) const{ return desorb_cut; }
	void PaprecaConfig::setDesorbDelMax( const int &desorb_delmax_in ){ desorb_delmax = desorb_delmax_in; }
	const int &PaprecaConfig::getDesorbDelMax( ) const{ return desorb_delmax; }
	void PaprecaConfig::setDesorptionStyle( const std::string &desorb_style_in ){ desorb_style = desorb_style_in; }
	const std::string &PaprecaConfig::getDesorptionStyle( ) const{ return desorb_style; }
	
	
	//Height calculation settings
	void PaprecaConfig::setHeightMethod( const std::string &height_method_in ){ height_method = height_method_in; }
	const std::string &PaprecaConfig::getHeightMethod( ) const{ return height_method; }
	void PaprecaConfig::setHeightPercentage( const double &height_percentage_in ){ height_percentage = height_percentage_in; }
	const double &PaprecaConfig::getHeightPercentage( ) const{ return height_percentage; }
	void PaprecaConfig::setBinWidth( const double &bin_width_in ){ bin_width = bin_width_in; }
	const double &PaprecaConfig::getBinWidth( ) const{ return bin_width; }
	
	
	//type2sigma
	void PaprecaConfig::initSigmasFromLammps( LAMMPS_NS::LAMMPS *lmp ){ initType2SigmaFromLammpsPairCoeffs( lmp , type2sigma ); }
	
	void PaprecaConfig::setSpeciesPair2Sigma( const int &species1 , const int &species2 , const double &sigma ){
		
		INT_PAIR pair( species1 , species2 );
		INT_PAIR pair_inverse( species1 , species2 );
		
		type2sigma[ pair ] = sigma;
		type2sigma[ pair_inverse ] = sigma;
		
		
	}
	
	void PaprecaConfig::setSigmaStyle( const std::string &sigmastyle_in ){ sigma_style = sigmastyle_in; }
	const std::string &PaprecaConfig::getSigmaStyle( )const{ return sigma_style; }
	void PaprecaConfig::setSigmaMix( const std::string &mixstyle_in ){ sigma_mix = mixstyle_in; }
	const std::string &PaprecaConfig::getSigmaMixStyle( )const{ return sigma_mix; }
	
	void PaprecaConfig::mixSigmas( LAMMPS_NS::LAMMPS *lmp ){

		if( !sigma_mix.empty( ) && sigma_mix != "no" ){
			
			int types_num = *( int *)lammps_extract_global( lmp , const_cast<char*>( "ntypes" ) ); //Obtain number of types to allocate/retrieve sigma array (containing sigma pairstyle coeffs).
			
			for( int i = 1; i < types_num + 1; ++i ){
				
				INT_PAIR iipair( i , i );
				
				if( !mappingExists( type2sigma , iipair ) || type2sigma[iipair] < std::numeric_limits< double >::epsilon( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Mixing attempted but diagonal term for atom type " + std::to_string( i ) + " was not set!" ); }
				double sigmaii = type2sigma[iipair];
				
				for( int j = 1; j < types_num + 1; ++j ){
					
					INT_PAIR jjpair( j , j );
					if( !mappingExists( type2sigma , jjpair ) || type2sigma[jjpair] < std::numeric_limits< double >::epsilon( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Mixing attempted but diagonal term for atom type " + std::to_string( i ) + " was not set!" ); }
					double sigmajj = type2sigma[jjpair];
					
					if( i != j ){
					
						INT_PAIR pair( i , j );	
						if( mappingExists( type2sigma , pair ) ){ warnAll( MPI_COMM_WORLD , "Mixing operation will overwrite manually set parameters for types " + std::to_string( i ) + " and " + std::to_string( j ) + ". \n You might want to review your sigmas_options!" ); }
						
						double sigmaij;
						if( sigma_mix == "geom" ){
							
							sigmaij = std::sqrt( sigmaii * sigmajj );
							
						}else if( sigma_mix == "arithm" ){
							sigmaij = 0.5 * ( sigmaii + sigmajj );
							
						}
						type2sigma[pair] = sigmaij;
						
					}
					
					
				}
				
			
			}
			
			
		}
		
		
	}
	
	const double PaprecaConfig::getSigmaFromAtomTypes( const int &atom1_type , const int &atom2_type ){
		
		INT_PAIR pair( atom1_type , atom2_type );
		double sigma = -1;
		
		if( mappingExists( type2sigma , pair ) ){
			
			sigma = type2sigma[pair];
			
		}else{
			
			allAbortWithMessage( MPI_COMM_WORLD , "Tried to return unmapped sigma values for species " + std::to_string( atom1_type ) + " and " + std::to_string( atom2_type ) + " in papreca_config.cpp." );
			
		}

		return sigma; //If this function returns -1, it indicates an error. In theory you should never reach this return, as the function will abort if the mapping does not exist.
		
	}
	
	const bool PaprecaConfig::type2SigmaMapIsEmpty( ) const{ return( type2sigma.empty( ) ? true : false ); }
	
	
	//Equilibrations - LAMMPS
	void PaprecaConfig::setMinimize1( const std::string &minimize1_in ){ minimize1 = minimize1_in; }
	const std::string &PaprecaConfig::getMinimize1( ) const{ return minimize1; }
	void PaprecaConfig::setMinimize2( const std::string &minimize2_in ){ minimize2 = minimize2_in; }
	const std::string &PaprecaConfig::getMinimize2( ) const{ return minimize2; }
	void PaprecaConfig::setTrajDuration( const int &traj_duration_in ){ traj_duration = traj_duration_in; }
	const int &PaprecaConfig::getTrajDuration( ) const{ return traj_duration; }
	void PaprecaConfig::setCtimeConvert( const double &c_convert_in ){ c_time_convert = c_convert_in; }
	const double &PaprecaConfig::getCtimeConvert( ){ return c_time_convert; }
	
	
	//Neighbor lists
	void PaprecaConfig::setNeibLists( const std::string &neiblist_half_in , const std::string &neiblist_full_in ){
			neiblist_half = neiblist_half_in;
			neiblist_full = neiblist_full_in;
	}
	
	const std::string &PaprecaConfig::getHalfNeibListName( ) const{ return neiblist_half; }
	const std::string &PaprecaConfig::getFullNeibListName( ) const{ return neiblist_full; }
	
	
	//ExportFiles
	Log &PaprecaConfig::getLogFile( ){ return log_file; }
	HeightVtime &PaprecaConfig::getHeightVtimeFile( ){ return heightVtime_file; }
	SurfaceCoverage &PaprecaConfig::getSurfaceCoverageFile( ){ return surfcoverage_file; }
	
	void PaprecaConfig::calcSurfaceCoverage( ){
		
		int depdata_global[2] = { 0 , 0 };
		int depdata_local[2] = { 0 , 0 };
		
		for( auto it = predefined_catalog.depositions_map.begin( ); it != predefined_catalog.depositions_map.end( ); ++it ){
			
			//Loop through all deposition events. Sum all deposition tries and sites of a single proc.
			PredefinedDeposition *depo = it->second;
			depdata_local[0] += depo->getDepositionSites( );
			depdata_local[1] += depo->getDepositionTries( );
			
		}
		
		//Different processors have different events. Hence, they should have different deptries and depsites.
		//Gather the different deposites and tries from ALL events
		
		MPI_Allreduce( depdata_local , depdata_global , 2 , MPI_INT , MPI_SUM , MPI_COMM_WORLD ); 
		
		if( depdata_global[1] == 0 ){ 
			surface_coverage = 0; //To avoid dividing by 0, set the surface coverage to zero if there are no detected deposition sites.
		}else{ 
			surface_coverage = static_cast< double >( 1.0 ) - static_cast< double >( depdata_global[0] ) / static_cast< double >( depdata_global[1] );  //depsites are the free sites for deposition. Hence, the surface coverage is 1 - sites/tries = 1 - sticking_coeff. This is the GLOBAL sticking coeff and not the per-event sticking coeff (as calculated by the calcStickingCoeff function).
		}

	}
	
	ElementalDistribution &PaprecaConfig::getElementalDistributionsFile( ){ return elementalDistribution_files; }
	ExecTime &PaprecaConfig::getExecTimeFile( ){ return execTime_file; }
	
	void PaprecaConfig::setupExportFiles( const int &proc_id ){
		
		if( proc_id == 0 ){ //Essential to open/write a/to the file USING ONE PROC ONLY! Opening a file with multiple procs at once can corrupt the file, lead to undefined behavior, or simply lead to doubly-appended lines.
			
			log_file.init( );
			if( heightVtime_file.isActive( ) ){ heightVtime_file.init( ); }
			if( surfcoverage_file.isActive( ) ){ surfcoverage_file.init( ); }
			if( execTime_file.isActive( ) ){ execTime_file.init( ); }
			
			//No need to init the ElementalDistributions file here as those are created during the run.
			
		}
		
	}
	
	void PaprecaConfig::setHybridStartTimeStamp4ExecTimeFile( const int &KMC_loopid ){
		
		if( execTime_file.isActive( ) && ( KMC_loopid % execTime_file.getPrintFreq( ) == 0 ) ){ execTime_file.setHybridStartTimeStamp( ); }
		
		
	}
	
	void PaprecaConfig::calcHybridAndKMCTimes4ExecTimeFile( const int &nprocs , const int &KMC_loopid ){
		
		if( execTime_file.isActive( ) && ( KMC_loopid % execTime_file.getPrintFreq( ) == 0 ) ){
			
			execTime_file.calcHybridTime( nprocs );
			execTime_file.calcKMCtime( nprocs );
			
		}
		
	}
	
	void PaprecaConfig::setMDTimeStamp4ExecTimeFile( const int &KMC_loopid ){
			
		if( execTime_file.isActive( ) && ( KMC_loopid % execTime_file.getPrintFreq( ) == 0 ) ){ execTime_file.setMDstartTimeStamp( ); }
		
	}
	
	void PaprecaConfig::calcMDTime4ExecTimeFile( const int &nprocs , const int &KMC_loopid ){
		
		if( execTime_file.isActive( ) && ( KMC_loopid % execTime_file.getPrintFreq( ) == 0 ) ){ execTime_file.calcMDtime( nprocs ); }
		
	}
	
	
	void PaprecaConfig::dumpElementalDistributionFile( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const int &KMC_loopid , double **mass_profiles_total , double *atom_mass , const int &bins_num , const int &types_num ){
		
		if( proc_id == 0 ){
			if( elementalDistribution_files.isActive( ) && ( KMC_loopid % elementalDistribution_files.getPrintFreq( ) == 0 ) ){
				elementalDistribution_files.init( KMC_loopid , types_num );
				elementalDistribution_files.append( lmp , mass_profiles_total , types_num , bins_num , bin_width , atom_mass );
				elementalDistribution_files.close( );
			}
		}
		
		
		
	}
	
	void PaprecaConfig::appendExportFiles( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const double &time , const char *event_type , const double &film_height , const int &KMC_loopid ){
	
		if( proc_id == 0 ){
			
			//Append to all files except for papreca.log (contains event specific information and is appended in the event_execute.cpp file)
			if( heightVtime_file.isActive( ) && ( KMC_loopid % heightVtime_file.getPrintFreq( ) == 0 ) ){ heightVtime_file.append( time , film_height ); }
			if( surfcoverage_file.isActive( ) && ( KMC_loopid % surfcoverage_file.getPrintFreq( ) == 0 ) ){ surfcoverage_file.append( time , surface_coverage ); } //surf_coverage here is a member variable of the papreca_config object. 
			if( execTime_file.isActive( ) && ( KMC_loopid % execTime_file.getPrintFreq( ) == 0 ) ){ 
				execTime_file.append( KMC_loopid , *( int *)lammps_extract_global( lmp , "natoms" ) );
			}
			
		}
		
	}
	
	void PaprecaConfig::closeExportFiles( const int &proc_id ){
		
		if( proc_id == 0 ){
			
			log_file.close( );
			if( heightVtime_file.isActive( ) ){ heightVtime_file.close( ); }
			if( surfcoverage_file.isActive( ) ){ surfcoverage_file.close( ); }
			if( execTime_file.isActive( ) ){ execTime_file.close( ); }
			
		}
		
		
	}
	
	void PaprecaConfig::setRestartDumpFreq( const int &restart_dumpfreq_in ){ restart_dumpfreq = restart_dumpfreq_in; }
	const int &PaprecaConfig::getRestartDumpFreq( ) const{ return restart_dumpfreq; }
	void PaprecaConfig::dumpLAMMPSRestart( LAMMPS_NS::LAMMPS *lmp , const int KMC_loopid ){
		
		dumpRestart( lmp , KMC_loopid , this->getRestartDumpFreq( ) );
		
		
	}
	
	
}//end of namespace PAPRECA
