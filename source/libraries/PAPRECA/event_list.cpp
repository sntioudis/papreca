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
/// @brief Definitions of PredefinedEvent classes.

#include "event_list.h"


namespace PAPRECA{
	
	//---------------------------------------PredefinedReaction class---------------------------------------
	//Constructors/Destructors
	PredefinedReaction::PredefinedReaction( const int &atom1_type_in , const int &atom2_type_in , const int &bond_type_in , const double &rate_in ) : atom1_type( atom1_type_in ) , atom2_type( atom2_type_in ) , bond_type( bond_type_in ) , rate( rate_in ){ }
	PredefinedReaction::PredefinedReaction( const int &atom1_type_in , const int &atom2_type_in , const int &bond_type_in , const double &rate_in , const std::vector< int > &catalyzing_types_in ) : atom1_type( atom1_type_in ) , atom2_type( atom2_type_in ) , bond_type( bond_type_in ) , rate( rate_in ) , catalyzing_types( catalyzing_types_in ){ }
	PredefinedReaction::~PredefinedReaction( ){ }
	
	//Functions
	const int &PredefinedReaction::getAtom1Type( ) const{ return atom1_type; }
	const int &PredefinedReaction::getAtom2Type( ) const{ return atom2_type; }
	const int &PredefinedReaction::getBondType( ) const{ return bond_type; }
	const double &PredefinedReaction::getRate( ) const{ return rate; }
	const std::vector< int > &PredefinedReaction::getCatalyzingTypes( ) const{ return catalyzing_types; }
	const bool &PredefinedReaction::isForm( ) const{ return is_form; }
	//------------------------------------END OF PredefinedReaction class------------------------------------
	
	

	//---------------------------------------PredefinedBondForm class---------------------------------------
	//Constructors/Destructors
	PredefinedBondForm::PredefinedBondForm( const int &atom1_type_in , const int &atom2_type_in , const int &bond_type_in , const double &rate_in , const double &bond_dist_sqr_in , const int &delete_atoms_in , const int &lone_candidates_in , const bool &same_mol_in ) : PredefinedReaction( atom1_type_in , atom2_type_in , bond_type_in , rate_in ) , bond_dist_sqr( bond_dist_sqr_in ) , delete_atoms( delete_atoms_in ) , lone_candidates( lone_candidates_in ) , same_mol( same_mol_in ){ is_form = true; }
	PredefinedBondForm::PredefinedBondForm( const int &atom1_type_in , const int &atom2_type_in , const int &bond_type_in , const double &rate_in , const double &bond_dist_sqr_in , const int &delete_atoms_in , const int &lone_candidates_in , const bool &same_mol_in , const std::vector< int > &catalyzing_types_in ) : PredefinedReaction( atom1_type_in , atom2_type_in , bond_type_in , rate_in , catalyzing_types_in ) , bond_dist_sqr( bond_dist_sqr_in ) , delete_atoms( delete_atoms_in ) , lone_candidates( lone_candidates_in ) , same_mol( same_mol_in ){ is_form = true; }
	PredefinedBondForm::~PredefinedBondForm( ){ }
	
	
	//Functions
	const double &PredefinedBondForm::getBondDistSqr( ) const{ return bond_dist_sqr; }
	const bool &PredefinedBondForm::isSameMol( ) const{ return same_mol; }
	const int &PredefinedBondForm::isDeleteAtoms( ) const{ return delete_atoms; }
	const int &PredefinedBondForm::isLone( ) const{ return lone_candidates; }
	//------------------------------------End of PredefinedBondForm class------------------------------------
	
	
	
	
	//---------------------------------------PredefinedDiffusionHop class---------------------------------------
	//Constructors/Destructors
	PredefinedDiffusionHop::PredefinedDiffusionHop( const int &parent_type_in , const double &insertion_vel_in , const double &diffusion_dist_in , const double &rate_in , const std::string &custom_style_in , const std::vector< int > &style_atomtypes_in ) : parent_type( parent_type_in ) , insertion_vel( insertion_vel_in ) , diffusion_dist( diffusion_dist_in ), rate( rate_in ) , custom_style( custom_style_in ) , style_atomtypes( style_atomtypes_in ) , diffused_type( parent_type_in ) , is_displacive( true ){ }
	PredefinedDiffusionHop::PredefinedDiffusionHop( const int &parent_type_in , const double &insertion_vel_in , const double &diffusion_dist_in , const double &rate_in , const std::string &custom_style_in , const std::vector< int > &style_atomtypes_in , const int &diffused_type_in , const bool &is_displacive_in ) : parent_type( parent_type_in ) , insertion_vel( insertion_vel_in ) , diffusion_dist( diffusion_dist_in ), rate( rate_in ) , custom_style( custom_style_in ) , style_atomtypes( style_atomtypes_in ) , diffused_type( diffused_type_in ) , is_displacive( is_displacive_in ){ }
	PredefinedDiffusionHop::~PredefinedDiffusionHop( ){ }
	
	//functions
	const int &PredefinedDiffusionHop::getParentAtomType( ) const{ return parent_type; }
	const double &PredefinedDiffusionHop::getInsertionVel( ) const{ return insertion_vel; }
	const int &PredefinedDiffusionHop::getDiffusedAtomType( ) const{ return diffused_type; }
	const bool &PredefinedDiffusionHop::isDisplacive( ) const{ return is_displacive; }
	const double &PredefinedDiffusionHop::getDiffusionDist( ) const{ return diffusion_dist; }
	const double &PredefinedDiffusionHop::getRate( ) const{ return rate; }
	const std::string &PredefinedDiffusionHop::getCustomStyle( ) const{ return custom_style; }
	const std::vector< int > &PredefinedDiffusionHop::getStyleAtomTypes( ) const{ return style_atomtypes; }
	//-------------------------------------End of PredefinedDiffusionHop class-------------------------------------
	
	
	
	//-------------------------------------------PredefinedDeposition-------------------------------------------
	//Constructors/Destructors
	PredefinedDeposition::PredefinedDeposition( LAMMPS_NS::LAMMPS *lmp , const int &parent_type_in , const double &rate_in , const double &depo_offset_in , const double &insertion_vel_in , const std::string &adsorbate_name_in ) : parent_type( parent_type_in ) , rate( rate_in )  , depo_offset( depo_offset_in ) , insertion_vel( insertion_vel_in ) , adsorbate_name( adsorbate_name_in ){ 
	
		int imol = lmp->atom->find_molecule( adsorbate_name.c_str( ) );
		coords = lmp->atom->molecules[imol]->dx; //See LAMMPS molecule.h and molecule.cpp for definitions
		atom_types = lmp->atom->molecules[imol]->type;
		atoms_num = lmp->atom->molecules[imol]->natoms;
		center = lmp->atom->molecules[imol]->center;
	
	
	}
	
	//Secondary PredefinedDeposition constructor for fixed sticking coefficient
	PredefinedDeposition::PredefinedDeposition( LAMMPS_NS::LAMMPS *lmp , const int &parent_type_in , const double &rate_in , const double &depo_offset_in , const double &insertion_vel_in , const std::string &adsorbate_name_in , const double &sticking_coeff_in ) : PredefinedDeposition( lmp , parent_type_in , rate_in , depo_offset_in , insertion_vel_in , adsorbate_name_in ){
		sticking_coeff = sticking_coeff_in;
		variable_sticking = false;
	};
	
	PredefinedDeposition::~PredefinedDeposition( ){ }
	
	//Functions
	const int &PredefinedDeposition::getParentType( ) const{ return parent_type; }
	const double &PredefinedDeposition::getRate( ) const{ return rate; }
	const double &PredefinedDeposition::getDepoOffset( ) const{ return depo_offset; }
	const double &PredefinedDeposition::getInsertionVel( ) const{ return insertion_vel; }
	const bool &PredefinedDeposition::hasVariableStickingCoeff( ) const{ return variable_sticking; }
	void PredefinedDeposition::setStickingCoeff( const double &sticking_coeff_in ){ sticking_coeff = sticking_coeff_in; }
	const double &PredefinedDeposition::getStickingCoeff( ) const{ return sticking_coeff; }
	void PredefinedDeposition::incrementDepositionTries( ){ ++deposition_tries; }
	const int &PredefinedDeposition::getDepositionTries( ) const{ return deposition_tries; }
	void PredefinedDeposition::resetDepositionTries( ){ deposition_tries = 0; }
	void PredefinedDeposition::incrementDepositionSites( ){ ++deposition_sites; }
	const int &PredefinedDeposition::getDepositionSites( ) const{ return deposition_sites; }
	void PredefinedDeposition::resetDepositionSites( ){ deposition_sites = 0; }
	void PredefinedDeposition::resetDepositionTriesAndSites( ){
		resetDepositionTries( );
		resetDepositionSites( );
	}
	const std::string &PredefinedDeposition::getAdsorbateName( ) const{ return adsorbate_name; }
	const int &PredefinedDeposition::getAtomsNum( ) const{ return atoms_num; }
	int *PredefinedDeposition::getAtomTypes( ){ return atom_types; }
	double *PredefinedDeposition::getCenter( ){ return center; }
	double **PredefinedDeposition::getCoords( ){ return coords; }
	//-------------------------------------End of PredefinedDeposition Class-------------------------------------
	
	//---------------------------------PredefinedMonoatomicDesorption Class--------------------------------
	//Constructors/Destructors
	PredefinedMonoatomicDesorption::PredefinedMonoatomicDesorption( const int &parent_type_in , const double &rate_in ) : parent_type( parent_type_in ) , rate( rate_in ){ }
	PredefinedMonoatomicDesorption::~PredefinedMonoatomicDesorption( ){ }
	
	//Functions
	const int &PredefinedMonoatomicDesorption::getParentAtomType( ) const{ return parent_type; }
	const double &PredefinedMonoatomicDesorption::getRate( ) const{ return rate; }
	//--------------------------------PredefinedMonoatomicDesorption Class--------------------------------
	
	
	//------------------------------------PredefinedEventsCatalog Class------------------------------------
	void PredefinedEventsCatalog::deletePredefinedDepositionsFromMap( ){
		
		std::vector< PredefinedDeposition* > depositions2del;
		
		for ( const auto &element : depositions_map ){ //Scan depositions to del. Two different atom types can be parent for deposition of the same adsorbate. Hence, need to gather deletion candidates and delete later.
													   //Comparison is easy and can be made through basic arithmetic pointer comparison
			PredefinedDeposition *deposition = element.second;
			bool deposition_is_del = false;
			
			
			for ( const auto &del_deposition : depositions2del ){
				
				if( del_deposition == deposition ){
					deposition_is_del = true;
					break;
				}
				
				
			}
			
			if ( !deposition_is_del ){
				
				depositions2del.push_back( deposition );
				
			}
			
					
			
		
		}
		
		
		for( const auto &del_deposition : depositions2del ){ //Loop through retrieved predefined deposition events and delete.
			
			delete del_deposition;
			
		}
		
		depositions_map.clear( );
		depositions_set.clear( );
		
		
		
	}
	
	void PredefinedEventsCatalog::deletePredefinedDiffusionsFromMap( ){

		for ( auto &element : diffusions_map ){
			
			
			if( element.second ){
				
				const int &diffusing_type = element.first;
				PredefinedDiffusionHop *diffusion = element.second;
				delete diffusion;
				diffusions_map[diffusing_type] = NULL;
			}
				
				
			
		
		}
		
		diffusions_map.clear( );
		diffusions_set.clear( );


	}
	
	void PredefinedEventsCatalog::deletePredefinedBondFormsFromMap( ){

		for ( auto &element : bond_forms_map ){
		
			if ( element.second ){ //This means that the mapped PredefinedReaction is not NULL (from previous object deletions)
			
				INT_PAIR type_pair = element.first;
				INT_PAIR type_pair_reverse = INT_PAIR( element.first.second , element.first.first ); //For reactions we map not only the actual pair but also the inverse pair.
															//Hence, we need to identify the inverse type pair and NULL the inverse pointer
															//to avoid deleting the same element twice
				PredefinedBondForm *bond_form = element.second;
				delete bond_form;
				
				bond_forms_map[ type_pair ] = NULL; //Map pairs to NULL to avoid double deletion! Do the same for the reverse pair for the same reasons
				bond_forms_map[ type_pair_reverse ] = NULL;

			}
				
		
		}
		
		bond_forms_map.clear( );
		bond_forms_set.clear( );



	}
	
	void PredefinedEventsCatalog::clearBondsMaxMap( ){ bonds_max.clear( ); }
	
	void PredefinedEventsCatalog::deletePredefinedBondBreaksFromMap( ){

		for ( const auto &element : bond_breaks_map ){
		
			if ( element.second ){
			
				const int &bond_type = element.first;
				PredefinedReaction *bond_break = element.second;
				delete bond_break;
				bond_breaks_map[bond_type] = NULL;
			
			}
		
		
		}
		
		bond_breaks_map.clear( );
		bond_breaks_set.clear( );


	}
	
	void PredefinedEventsCatalog::deletePredefinedMonoatomicDesorptionsFromMap( ){
		
		for( const auto &element : monodes_map ){
			
			if( element.second ){
				
				const int &parent_type = element.first;
				PredefinedMonoatomicDesorption *des = element.second;
				delete des;
				monodes_map[parent_type] = NULL;
				
			}
			
			
		}
		
		monodes_map.clear( );
		monodes_set.clear( );
		
		
	}
	
	
	
	PredefinedEventsCatalog::PredefinedEventsCatalog( ){ };
	PredefinedEventsCatalog::~PredefinedEventsCatalog( ){
		
		PredefinedEventsCatalog::deletePredefinedDepositionsFromMap( );
		PredefinedEventsCatalog::deletePredefinedDiffusionsFromMap( );
		PredefinedEventsCatalog::deletePredefinedBondFormsFromMap( );
		PredefinedEventsCatalog::clearBondsMaxMap( );
		PredefinedEventsCatalog::deletePredefinedBondBreaksFromMap( );
		PredefinedEventsCatalog::deletePredefinedMonoatomicDesorptionsFromMap( );

	}
	//-------------------------------End of PredefinedEventsCatalog Class-------------------------------


}//end of namespace PAPRECA
