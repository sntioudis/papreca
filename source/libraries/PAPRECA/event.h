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
/// @brief Declarations for PAPRECA:Event class and its children classes.
///

#ifndef EVENT_H
#define EVENT_H

//System Headers
#include <array>
#include <vector>
#include <mpi.h>


//LAMMPS headers
#include "lammps.h"
/// \cond
#include "pointers.h"
/// \endcond

//KMC headers
#include "event_list.h"
#include "lammps_wrappers.h"
#include "utilities.h"

namespace PAPRECA{
	
	class Event{
	
		/// @class PAPRECA::Event
		/// @brief Parent class for PAPRECA events.
		///
		/// DISCLAIMER: PAPRECA::Event objects are not to be confused with PredefinedEventList objects. For more information about PredefinedEvent objects please refer to the relevant header and/or the documentation (event_list.h and event_list.cpp).
		/// PAPRECA::Event objects are the events discovered and (potentially) executed by the kMC algorithm. Every MPI proc can (and will) discover different PAPRECA::Event objects.
		/// The types of PAPRECA::Event supported by the initial version of this software are: PAPRECA::BondBreak, PAPRECA::BondForm, PAPRECA::Deposition (i.e., adsorption), PAPRECA::DiffusionHop, and PAPRECA::MonoatomicDesorption.
		/// It should be straightforward to extend the software to include new PAPRECA:Event types. A new child PAPRECA::Event class can be developed to enable the discovery of new classes of events.
		///	New event style implementations should be placed in this header.
	
		public:
		
		
			//Constructors/Destructors
			Event( );
			Event( const double &rate_in , const std::string &type_in );
			virtual ~Event( ); // Virtual destructor to ensure proper clean-up of derived classes
			
			//functions
			void assignRate( const double &rate_in );
			void assignType( const std::string &type_in );
			const double &getRate( ) const;
			void setRate( const double &rate_in );
			const std::string &getType( ) const;
			
			//Static functions
			static void fillRatesArr( double *event_rates , const std::vector< Event* > &events );
			static void fillRatesVec( std::vector< double > &event_rates , const std::vector< Event* > &events );
			static std::vector< double > getRatesVec( const std::vector< Event* > &events ); //This is a static function (i.e., is not bound to an instance of class Event). You can provide an events* vector here to get a (copy) rates vec.
			static double getSumOfRates( const std::vector< Event* > &events ); //This function uses the getRatesVec function and calculates the sum of all rates
			static void deleteAndClearLocalEvents( LAMMPS_NS::LAMMPS *lmp , std::vector<Event*> &events_local );
			
		protected:
			double rate;
			std::string type;
	
	};
	
	
	class Reaction : public Event{
		
		/// @class PAPRECA::Reaction
		/// @brief Child of PAPRECA:Event and parent of PAPRECA::BondBreak and PAPRECA::BondForm classes.
		
		public:
			//Child class constructor/destructor
			Reaction( const double &rate_in , const LAMMPS_NS::tagint &atom1id_in , const LAMMPS_NS::tagint &atom2id_in , const int &bond_type_in );
			~Reaction( );
			
			//Functions
			const LAMMPS_NS::tagint &getAtom1ID( ) const;
			const LAMMPS_NS::tagint &getAtom2ID( ) const;
			const int &getBondType( ) const;
			void initialize( const LAMMPS_NS::tagint &atom1id_in , const LAMMPS_NS::tagint &atom2id_in , const int &bond_type_in , const double &rate_in);
			void resetEvent( );
			void assignAtom1( const LAMMPS_NS::tagint &atom1id_in );
			void assignAtom2( const LAMMPS_NS::tagint &atom2id_in );
			void assignBondType( const int &bond_type_ind );
			
			
		protected:
			LAMMPS_NS::tagint atom1id = -1 , atom2id = -2;
			int bond_type = -3;
	
	};
	
	
	class BondBreak : public Reaction{
	
		/// @class PAPRECA::BondBreak
		/// @brief child of PAPRECA::Reaction dedicated to bond-breaking events.
		
		public:
			//Child class constructor/destructor
			BondBreak( const double &rate_in , const LAMMPS_NS::tagint &atom1id_in , const LAMMPS_NS::tagint &atom2id_in , const int &bond_type_in , PredefinedReaction *break_teamplate_in );
			~BondBreak( );
			
			//Functions
			PredefinedReaction *getBreakTemplate( );
			
			
			private:
				PredefinedReaction *break_template = NULL;
	
	};
	
	
	class BondForm : public Reaction{
	
		/// @class PAPRECA::BondForm
		/// @brief child of PAPRECA::Reaction dedicated to bond formation events.
		
		public:
			//Child class constructor/destructor
			BondForm( const double &rate_in , const LAMMPS_NS::tagint &atom1id_in , const LAMMPS_NS::tagint &atom2id_in , const int &bond_type_in , PredefinedBondForm *form_template_in );
			~BondForm( );
			
			//Functions
			PredefinedBondForm *getFormTemplate( );
			
			
			private:
				PredefinedBondForm *form_template = NULL;
	
	};
	
	
	class Deposition : public Event{
	
		/// @class PAPRECA::Deposition
		/// @brief Child of PAPRECA::Event dedicated to monoatomic or molecular adsorption.
		
		public:
		
			//Child class constructor/destructor
			Deposition( const double &rate_in , const double site_pos_in[3] , const double rot_pos_in[3] , const double &rot_theta_in , const int &mol_id_in , const std::string &mol_name_in , PredefinedDeposition *depo_template_in );
			~Deposition( );
			
			//Functions
			double *getSitePos( );
			double *getRotPos( );
			const double &getRotTheta( ) const;
			const int &getMolId( ) const;
			const std::string &getMolName( ) const;
			PredefinedDeposition *getDepoTemplate( );
		
		
		protected:
			double site_pos[3];
			double rot_pos[3]; //Required by create_atoms command with the mol option, defines the centre of rotation of the inserted molecule (angle in degrees).
			double rot_theta;
			int mol_id;
			std::string mol_name = "NONE";
			
			
		private:
			PredefinedDeposition *depo_template = NULL; //This is the associated predefined event to the discovered event. No need to new/delete anything here as that memory block is managed by the PredefinedDepositions class.
	
	};
	
	class Diffusion : public Event{
		
		/// @class PAPRECA::Diffusion
		/// @brief child of PAPRECA::Event class dedicated to diffusion events.
		
		public:
			//Child class constructor/destructor
			Diffusion( const double &rate_in , const double vacancy_pos_in[3] , const LAMMPS_NS::tagint &parent_id_in , const int &parent_type_in , const int &is_displacive_in , const int &diffused_type_in , PredefinedDiffusionHop *diff_template_in );
			~Diffusion( );
			
			//Functions
			double *getVacancyPos( );
			const LAMMPS_NS::tagint &getParentId( ) const;
			const int &getParentType( ) const;
			int isDisplacive( );
			const int &getDiffusedType( ) const;
			PredefinedDiffusionHop *getDiffTemplate( );
		
		protected:
			double vacancy_pos[3];
		
			LAMMPS_NS::tagint parent_id = -1;
			int parent_type = -2;	
			int is_displacive = 0;  ///< Displacive diffusion events "displace" the parent atom. This is defined as an integer (and not as bool, even though in realitiy, it is a bool) to facilitate the transfering of data using MPI.
			int diffused_type = -3;
		
		private:
			PredefinedDiffusionHop *diff_template = NULL;
		
		
	};
	
	class MonoatomicDesorption : public Event{
		
		/// @class PAPRECA::MonoatomicDesorption
		/// @brief child of PAPRECA::Event dedicated to monoatomic desorption events (i.e., where a single atom is ejected from the system).
		
		public:
			//Child class constructor/destructor
			MonoatomicDesorption( const double &rate_in , const LAMMPS_NS::tagint &parent_id_in , const int &parent_type_in , PredefinedMonoatomicDesorption *des_template_in );
			~MonoatomicDesorption( );
			
			//Functions
			const LAMMPS_NS::tagint &getParentId( ) const;
			const int &getParentType( ) const;
			PredefinedMonoatomicDesorption *getMonoDesTemplate( );
			
		protected:
			LAMMPS_NS::tagint parent_id = -1;
			int parent_type = -2;
			
		private:
			PredefinedMonoatomicDesorption *monodes_template = NULL;	
		
	};

}//end of PAPRECA namespace 


#endif
