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
/// @brief Definitions for PAPRECA:Event class and its children classes.

#include "event.h"

namespace PAPRECA{



	//---------------------------------------PARENT Event class---------------------------------------
	//Constructors/Destructors
	Event::Event( ): rate( 0.0 ), type( "NONE" ){ }
	Event::Event( const double &rate_in , const std::string &type_in ): rate( rate_in ), type( type_in ){ }
	Event::~Event( ){ };
	
	//Member functions
	void Event::assignRate( const double &rate_in ){ rate = rate_in; }
	void Event::assignType( const std::string &type_in ){ type = type_in; }
	const double &Event::getRate( )const{ return rate; }
	void Event::setRate( const double &rate_in ){ rate = rate_in; }
	const std::string &Event::getType( )const{ return type; }
	
	//Static functions
	void Event::fillRatesArr( double *event_rates , const std::vector< Event* > &events ){
		
		/// Receives a vector of pointers to event objects, and fills a (previously-initialized) double array with the rates of events.
		/// @param[in] event_rates previously-initialized array of doubles.
		/// @param[in] events vector to pointers of PAPRECA::Event objects (and/or objects of children of PAPRECA::Event).
		/// @see Event::fillRatesVec()
		/// @note This function assumes that the event_rates array was properly initialized (i.e., has as many elements as the number of events in the vector of pointers to PAPRECA::Event objects. Incorrectly initializing the array of doubles will lead to segmentation faults.
		
		
		for( int i = 0; i < events.size( ); ++i ){
			event_rates[i] = events[i]->getRate( );
		
		}
		
	}
	
	void Event::fillRatesVec( std::vector< double > &event_rates , const std::vector< Event* > &events ){
			
		/// Same as Event::fillRatesArr() but fills a vector of doubles instead of an array of doubles.
		/// @param[in] event_rates vector of doubles storing the rates of events.
		/// @param[in] events vector to pointers of PAPRECA::Event objects (and/or objects of children of PAPRECA::Event).
		/// @see Event::fillRatesArr()
		/// @note Even though a vector is used here, the present function can still lead to segmentation faults, since it uses C-style array access (i.e., event_rates[i]) instead of C++-style push_back functions. The use of C-style array access was mandatory to facilitate the communicator of data between MPI processes.
		
		for( int i = 0; i < events.size( ); ++i ){
			event_rates[i] = events[i]->getRate( );
				
		}
		
	}
	
	std::vector< double > Event::getRatesVec( const std::vector< Event* > &events ){
		
		/// Receives a vector of pointer to PAPRECA::Event objects and returns a vector of doubles of their corresponding rates.
		/// @param[in] events vector to pointers of PAPRECA::Event objects (and/or objects of children of PAPRECA::Event).
		/// @return vector of double rates of the input events.
		
		std::vector< double > rates;
		rates.reserve( events.size( ) );
		
		for( const auto &event : events ){
			
			rates.push_back( event->getRate( ) );
			
		}
		
		return rates;
		
	}
	
	double Event::getSumOfRates( const std::vector< Event* > &events ){ //This function uses the getRatesVec function and calculates the sum of all rates.
		std::vector< double > rates = Event::getRatesVec( events );
		return getSumOfVecElements( rates );
		
	}
	
	void  Event::deleteAndClearLocalEvents( LAMMPS_NS::LAMMPS *lmp , std::vector<Event*> &events_local ){
		
		/// Deletes all local events (different on each MPI process) and clears the events_local vector of PAPRECA::Event objects.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in,out] events_local vector containing all the PAPRECA::Event objects for a specific MPI process.
		
		
		for( const auto &event : events_local ){
			
			if( event->getType( ) == "RXN-BREAK" ){
				BondBreak *bond_break = dynamic_cast<BondBreak*>( event );
				delete bond_break;		
			}else if( event->getType( ) == "RXN-FORM" ){
				BondForm *bond_form = dynamic_cast<BondForm*>( event );
				delete bond_form;	
			}else if( event->getType( ) == "DEPO" ){
				Deposition *deposition = dynamic_cast<Deposition*>( event );
				delete deposition;			
			}else if( event->getType( ) == "DIFF" ){
				Diffusion *diffusion = dynamic_cast<Diffusion*>( event );
				delete diffusion;
			}else if( event->getType( ) == "MONO-DES" ){
				MonoatomicDesorption *monodes = dynamic_cast<MonoatomicDesorption*>( event );
				delete monodes;
			}
		
		}
		
		//Reset the vector and prepare for the next kMC step
		events_local.clear( ); //Clear the vector

	}
	
	//---------------------------------------END OF PARENT Event Class---------------------------------------
	
	
	
	
	//---------------------------------------------Reaction class---------------------------------------------
	//Child class constructor
	Reaction::Reaction( const double &rate_in , const LAMMPS_NS::tagint &atom1id_in , const LAMMPS_NS::tagint &atom2id_in , const int &bond_type_in ): Event::Event( rate_in , "RXN" ) , atom1id( atom1id_in ) , atom2id( atom2id_in ) , bond_type( bond_type_in ){ }
	Reaction::~Reaction( ){ }
	
	//Functions
	const LAMMPS_NS::tagint &Reaction::getAtom1ID( )const{ return atom1id; }
	const LAMMPS_NS::tagint &Reaction::getAtom2ID( )const{ return atom2id; }
	const int &Reaction::getBondType( )const{ return bond_type; }
	void Reaction::assignAtom1( const LAMMPS_NS::tagint &atom1id_in ){ atom1id = atom1id_in; }
	void Reaction::assignAtom2( const LAMMPS_NS::tagint &atom2id_in ){ atom2id = atom2id_in; }
	void Reaction::assignBondType( const int &bond_type_in ){ bond_type = bond_type_in; }
	
	void Reaction::initialize( const LAMMPS_NS::tagint &atom1id_in , const LAMMPS_NS::tagint &atom2id_in , const int &bond_type_in , const double &rate_in ){
	
		atom1id = atom1id_in;
		atom2id = atom2id_in;
		bond_type = bond_type_in;
		
		rate = rate_in;
		
	
	}
	
	void Reaction::resetEvent( ){
	
		atom1id = -1;
		atom2id = -2;
		bond_type = -3;
		
		rate = 0.0;
	
	}
	//--------------------------------------------END OF Reaction Class------------------------------------------
	
	
	
	
	//----------------------------------------------BondBreak class----------------------------------------------
	//Constructors/Destructors
	BondBreak::BondBreak( const double &rate_in , const LAMMPS_NS::tagint &atom1id_in , const LAMMPS_NS::tagint &atom2id_in , const int &bond_type_in , PredefinedReaction *break_template_in ): Reaction::Reaction( rate_in , atom1id_in , atom2id_in , bond_type_in ) , break_template( break_template_in ){
		type += "-BREAK"; //Update initialized reaction type from RXN to RXN-BREAK
	}
	BondBreak::~BondBreak( ){ }
	
	//Functions
	PredefinedReaction *BondBreak::getBreakTemplate( ){ return break_template; }
	//----------------------------------------END OF CHILD BondBreak Class----------------------------------------
	
	
	
	
	
	
	//-----------------------------------------------BondForm Class-----------------------------------------------
	//Constructors/Destructors
	BondForm::BondForm( const double &rate_in , const LAMMPS_NS::tagint &atom1id_in , const LAMMPS_NS::tagint &atom2id_in , const int &bond_type_in , PredefinedBondForm *form_template_in ): Reaction::Reaction( rate_in , atom1id_in , atom2id_in , bond_type_in ) , form_template( form_template_in ){
		type += "-FORM"; //Update initialized reaction type from RXN to RXN-FORM
	}
	BondForm::~BondForm( ){ }
	
	//Functions
	PredefinedBondForm *BondForm::getFormTemplate( ){ return form_template; }
	//----------------------------------------End of CHILD BondForm Class----------------------------------------
	
	
	
	
	
	//------------------------------------------CHILD Deposition Class------------------------------------------
	//Constructors/Destructors
	Deposition::Deposition( const double &rate_in , const double site_pos_in[3] , const double rot_pos_in[3] , const double &rot_theta_in , const int &mol_id_in , const std::string &mol_name_in , PredefinedDeposition *depo_template_in ) :  Event::Event( rate_in , "DEPO" ) , rot_theta( rot_theta_in ) , mol_id( mol_id_in ) , mol_name( mol_name_in ) , depo_template( depo_template_in ){
		copyDoubleArray3D( site_pos , site_pos_in );
		copyDoubleArray3D( rot_pos , rot_pos_in );
	}
	Deposition::~Deposition( ){ };
	
	//Functions
	double *Deposition::getSitePos( ){ return site_pos; };
	double *Deposition::getRotPos( ){ return rot_pos; };
	const double &Deposition::getRotTheta( )const{ return rot_theta; };
	const int &Deposition::getMolId( )const{ return mol_id; };
	const std::string &Deposition::getMolName( )const{ return mol_name; };
	PredefinedDeposition *Deposition::getDepoTemplate( ){ return depo_template; }
	//---------------------------------------END OF CHILD Deposition Class---------------------------------------
	
	
	
	
	//-------------------------------------------CHILD Diffusion CLASS-------------------------------------------
	//Constructors/Destructors
	Diffusion::Diffusion( const double &rate_in , const double vacancy_pos_in[3] , const LAMMPS_NS::tagint &parent_id_in , const int &parent_type_in , const int &is_displacive_in , const int &diffused_type_in , PredefinedDiffusionHop *diff_template_in ) : Event::Event( rate_in , "DIFF" ) , parent_id( parent_id_in ), parent_type( parent_type_in ) , is_displacive( is_displacive_in ) , diffused_type( diffused_type_in ) , diff_template( diff_template_in ){ 
		copyDoubleArray3D( vacancy_pos , vacancy_pos_in );
	};
	Diffusion::~Diffusion( ){ };
	
	//Functions
	double *Diffusion::getVacancyPos( ){ return vacancy_pos; };
	const LAMMPS_NS::tagint &Diffusion::getParentId( )const{ return parent_id; };
	const int &Diffusion::getParentType( )const{ return parent_type; };
	int Diffusion::isDisplacive( ){ return is_displacive; };
	const int &Diffusion::getDiffusedType( )const{ return diffused_type; };
	PredefinedDiffusionHop *Diffusion::getDiffTemplate( ){ return diff_template; }
	//-------------------------------------------END OF CHILD Diffusion CLASS-------------------------------------------
	
	//-----------------------------------------CHILD MonoatomicDesorption Class-----------------------------------------
	//Constructors/Destructors
	MonoatomicDesorption::MonoatomicDesorption( const double &rate_in , const LAMMPS_NS::tagint &parent_id_in , const int &parent_type_in , PredefinedMonoatomicDesorption *des_template_in ) : Event::Event( rate_in , "MONO-DES" ) , parent_id( parent_id_in ) , parent_type( parent_type_in ) , monodes_template( des_template_in ){ }
	MonoatomicDesorption::~MonoatomicDesorption( ){ }
	
	//Functions
	const int &MonoatomicDesorption::getParentId( ) const{ return parent_id; }
	const int &MonoatomicDesorption::getParentType( ) const{ return parent_type; }
	PredefinedMonoatomicDesorption *MonoatomicDesorption::getMonoDesTemplate( ){ return monodes_template; }
	//-------------------------------------END OF CHILD MonoatomicDesorption CLASS-------------------------------------
	
	
} //End of PAPRECA Namespace
