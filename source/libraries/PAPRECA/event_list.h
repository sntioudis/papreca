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
/// @brief Declarations of PredefinedEvent classes.
///
/// DISCLAIMER: Predefined Event classes are not to be confused with PAPRECA::Event objects. For more information about Event objects please refer to the relevant header and/or the documentation (event.h and event.cpp).
/// Objects of Predefined Events classes are associated with the Predefined Catalog of kMC events. Objects of PredefinedEvent classes store information (e.g., rate of events, participating atom types, mol names, etc.) regarding ACCEPTABLE kMC events.
/// ALL MPI PROCS are initialized (by the input file) WITH THE SAME Predefined Event objects that do not change during the run (i.e., the Catalog of Events is PREDEFINED).
/// It should be straightforward to extend the software to include new PREDEFINED events. A new PREDEFINED event class can be developed to enable the management of new styles of predefined events by the PredefinedEventsCatalog class.
/// New Predefined Event style implementations should be placed in this header.

#ifndef EVENT_LIST_H
#define EVENT_LIST_H

//system headers
#include <utility>
#include <unordered_map>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <mpi.h>

//LAMMPS Headers
#include "lammps.h"
/// \cond
#include "atom.h"
#include "pointers.h"
#include "molecule.h"
/// \endcond

//kMC Headers
#include "utilities.h"


namespace PAPRECA{

	class PredefinedReaction{
		
		/// @class PAPRECA::PredefinedReaction
		/// @brief class associated with reaction events. This is a Predefined template for PAPRECA::BondBreak and PAPRECA::BondForm events.
		
		public:
			//Constructors/Destructors. The user can use the first constructor if no catalyzing_type(s) are defined. We know that a reaction is NOT catalyzed if the catalyzing_types vector is empty.
			PredefinedReaction( const int &atom1_type_in , const int &atom2_type_in , const int &bond_type_in , const double &rate_in  );
			PredefinedReaction( const int &atom1_type_in , const int &atom2_type_in , const int &bond_type_in , const double &rate_in , const std::vector< int > &catalyzing_types_in );
			~PredefinedReaction( );
			
			//Functions
			const int &getAtom1Type( ) const;
			const int &getAtom2Type( ) const;
			const int &getBondType( ) const;
			void setLengthEquil( const double &length_equil_in );
			const double &getLengthEquil( ) const;
			void setLimitLowSqr( const double &limit_low_sqr_in );
			const double &getLimitLowSqr( ) const;
			void setLimitHighSqr( const double &limit_high_sqr_in );
			const double &getLimitHighSqr( ) const;
			void setSqrLimits( const double &length_equil_in , const double &length_perc_in );
			const double &getRate( ) const;
			const std::vector< int > &getCatalyzingTypes( ) const;
			const bool &isForm( ) const;
		
		protected:
			int atom1_type = -1;
			int atom2_type = -2;
			int bond_type = -3;
			double length_equil = 0.0; //Optional:  Equilibrium bond length
			double limit_low_sqr = 0.0; //Optional: Acceptable event only if bond length>= length_limit
			double limit_high_sqr = 0.0;
			double rate = 0.0;
			std::vector< int > catalyzing_types; ///< these types of atoms have to be present in the neighborhood of atom1 or atom2. Otherwise, the reaction event is invalid. If the catalyzing_types vector is empty, the neighborhood of atom1 and atom2 is not scanned for catalyzing types.
			bool is_form = false; ///< this variable allows us to differentiate between the parent class PAPRECA::PredefinedReaction and the child class PAPRECA::PredefinedBondForm
			
			
	};
	
	class PredefinedBondForm : public PredefinedReaction{
	
		/// @class PAPRECA::PredefinedBondForm
		/// @brief child class of PAPRECA::PredefinedReaction associated with bond formation events. This is a Predefined template for PAPRECA::BondForm events.
		
		public:
			//Constructors/Destructors
			PredefinedBondForm( const int &atom1_type_in , const int &atom2_type_in , const int &bond_type_in , const double &rate_in , const double &bond_dist_sqr_in , const int &delete_atoms_in , const int &lone_candidates_in , const bool &same_mol_in );
			PredefinedBondForm( const int &atom1_type_in , const int &atom2_type_in , const int &bond_type_in , const double &rate_in , const double &bond_dist_sqr_in , const int &delete_atoms_in , const int &lone_candidates_in , const bool &same_mol_in , const std::vector< int > &catalyzing_types_in );
			~PredefinedBondForm( );
			
			//Functions
			const double &getBondDistSqr( ) const;
			const bool &isSameMol( ) const; 
			const int &isDeleteAtoms( ) const;
			const int &isLone( ) const;
			
		protected:
			double bond_dist_sqr = 0.0; ///< if the distance between atom1 and atom2 is less than or equal to bond_dist then the formation event is valid.
			bool same_mol = false; ///< if true, then bond formation can occur even if the molecule IDs of atom1 and atom2 are the same.
			int delete_atoms = 0;  ///< if delete_atoms=1, then atom1 and atom2 are deleted after the formation event. This variable is "technically boolean" but defined as an int to facilitate the communication of data among MPI processes (MPI_INT is supported by the MPI protocol).
			int lone_candidates = 0; ///< Lones are atoms without any bonds. If line_candidates=true, then both atom1 and atom2 have to be "lone". Otherwise, the formation event is invalid.
	};


	class PredefinedDiffusionHop{
		
		/// @class PAPRECA::PredefinedDiffusionHop
		/// @brief class associated with diffusion events. This is a Predefined template for PAPRECA::DiffusionHop events.
		
		public:
			//Constructor/Destructors
			PredefinedDiffusionHop( const int &parent_type_in , const double &insertion_vel_in , const double &diffusion_dist_in , const double &rate_in , const std::string &custom_style_in , const std::vector< int > &style_atomtypes_in ); //Initializing a PredefinedDiffusionHop object without a diffused type immediately sets the diffusion type to displacive (i.e., the parent atom moves).
			PredefinedDiffusionHop( const int &parent_type_in , const double &insertion_vel_in , const double &diffusion_dist_in , const double &rate_in , const std::string &custom_style_in , const std::vector< int > &style_atomtypes_in , const int &diffused_type_in , const bool &is_displacive_in ); //This constructor gives the freedom to set diffused type as well as type of diffusion hop (i.e., displacive==parent atom moves, non-displacive==new atom spawns on vacant site).
			~PredefinedDiffusionHop( );
			
			//Functions
			const int &getParentAtomType( ) const;
			const double &getInsertionVel( ) const;
			const int &getDiffusedAtomType( ) const;
			const bool &isDisplacive( ) const; ///< Displacive events "move" the parent atom, while non-displacive events, create a new atom at the vacancy.
			const double &getDiffusionDist( ) const;
			const double &getRate( ) const;
			const std::string &getCustomStyle( ) const;
			const std::vector< int > &getStyleAtomTypes( ) const;
			
		private:
			int parent_type = -1; ///< type of candidate parent atom.
			double insertion_vel = 0.0; ///< velocity of the displaced atom.
			int diffused_type = -2; ///< type of diffused atom. This can be the same as the parent type. Alternatively, it can be set to a different type if you want the atom type to change after diffusion.
			bool is_displacive = false; //You can have displacive diffusion and the parent type being different from the diffused type, if you want the atom type to change when it moves.
			double diffusion_dist = 0.0; /// displace atom by that much.
			std::string custom_style = "NONE"; ///< currently, only the Fe_4PO4neib style is an acceptable custom style for diffusion.
			std::vector< int > style_atomtypes; ///< vector that allows you to pass information about atom types from the PAPRECA input file.
			double rate = 0.0;
	
	};
	
	class PredefinedDeposition{
		
		/// @class PAPRECA::PredefinedDeposition
		/// @brief class associated with deposition events. This is a Predefined template for PAPRECA::Deposition events.
		public:
			//Constructor/Destructor
			PredefinedDeposition( LAMMPS_NS::LAMMPS *lmp , const int &parent_type_in , const double &rate_in , const double &depo_offset_in , const double &insertion_vel_in , const std::string &adsorbate_name_in );
			PredefinedDeposition( LAMMPS_NS::LAMMPS *lmp , const int &parent_type_in , const double &rate_in , const double &depo_offset_in , const double &insertion_vel_in , const std::string &adsorbate_name_in , const double &sticking_coeff_in );
			~PredefinedDeposition( );
			
			//Functions
			const int &getParentType( ) const;
			const double &getRate( ) const;
			const double &getDepoOffset( ) const;
			const double &getInsertionVel( ) const;
			const bool &hasVariableStickingCoeff( ) const;
			const double &getStickingCoeff( ) const;
			void setStickingCoeff( const double &sticking_coeff_in );
			void incrementDepositionTries( );
			const int &getDepositionTries( ) const;
			void resetDepositionTries( );
			void incrementDepositionSites( );
			const int &getDepositionSites( ) const;
			void resetDepositionSites( );
			void resetDepositionTriesAndSites( );
			const std::string &getAdsorbateName( ) const;
			const int &getAtomsNum( ) const;
			int *getAtomTypes( );
			double *getCenter( ) ;
			double **getCoords( );
			
		
		private:
			int parent_type = -1; ///< type of candidate parent atom.
			double rate = 0.0;
			double depo_offset = 0.0; ///< the molecule/particle is placed at a distance of depo_offset above the parent candidate atom.
			double insertion_vel = 0.0; ///< velocity of inserted molecule/particle.
			bool variable_sticking = true; ///< true or false for variable and fixed sticking coefficients, respectively.
			int deposition_tries = 0; ///< total sites (i.e., free + occupied) for that specific deposition event.
			int deposition_sites = 0; ///< available deposition sites for that specific deposition event.
			double sticking_coeff = -1.0; ///< sticking coefficient (i.e., free sites/ total sites). Does not change its value if variable_sticking=false.
			std::string adsorbate_name = "NONE"; ///< name of adsorbate. This has to be identical to the adsorbate name as initialized in the LAMMPS input file (e.g., if this command is used: "molecule mmmTCP ./TCP.xyz"), then your adsorbate name should be mmmTCP. 
			
			//Lammps template definitions
			int atoms_num = 0; ///< Number of molecule atoms in LAMMPS (and in the xyz molecule/particle input file).
			int *atom_types = NULL; ///< Template molecule types from lammps. This is an int[mol_natoms] array.
			double *center = NULL; //< Template molecule center coordinates array from lammps. This is a pointer to a double[3] array (defined in the LAMMPS header molecule.h).
			double **coords = NULL; ///< Template molecule xyz (coordinates) array from lammps. This is a pointer to double[mol_natoms][3].
	
	
	};
	
	class PredefinedMonoatomicDesorption{
		
		/// @class PAPRECA::PredefinedMonoatomicDesorption.
		/// @brief associated with monoatomic desorption events. This is a Predefined template for PAPRECA::MonoatomicDesorption events.
		public:
			//Constructors/Destructors
			PredefinedMonoatomicDesorption( const int &parent_type_in , const double &rate_in );
			~PredefinedMonoatomicDesorption( );
			
			//Functions
			const int &getParentAtomType( ) const;
			const double &getRate( ) const;
			
			private:
				int parent_type = -1;
				double rate = 0.0;
				
	};
	
	//Typedefs
	typedef std::unordered_map< int , PredefinedReaction* > TYPE2REACTION_MAP;
	typedef std::unordered_map< INT_PAIR , PredefinedBondForm* , PairHash > PAIR2BONDFORM_MAP;
						 
	typedef std::unordered_map< int , PredefinedDiffusionHop* > TYPE2DIFFUSION_MAP;
	typedef std::unordered_map< int , PredefinedDeposition* > TYPE2DEPOSITION_MAP;
	typedef std::unordered_map< int , PredefinedMonoatomicDesorption* > TYPE2MONODES_MAP;
	
	
	class PredefinedEventsCatalog{
		
		/// @class PAPRECA::PredefinedEventsCatalog
		/// @brief General class that stores ALL the predefined events in the system. You can consider this as the Predefined Catalog of Events.
		///
		/// Notice that in many cases we use a set and a map, even though the map/set keys are identical.
		/// We do that to search on sets but retrieve values from the map. This is faster than iterating a map and getting mapped values from the map.
		/// The drawback is that it requires a little bit more RAM (as you store the same key in a set AND a map).
			
		private:
		
			friend class PaprecaConfig; //PaprecaConfig is a friend class as it accesses some of those private members to return predefined events to the main function.
		
			//Functions used by the destructor
			void deletePredefinedDepositionsFromMap( );
			void deletePredefinedDiffusionsFromMap( );
			void deletePredefinedBondFormsFromMap( );
			void clearBondsMaxMap( );
			void deletePredefinedBondBreaksFromMap( );
			void deletePredefinedMonoatomicDesorptionsFromMap( );
			
			//Constructors/Destructors
			PredefinedEventsCatalog( );
			~PredefinedEventsCatalog( );
			
			//PredefinedReactions
			INT_SET bond_breaks_set; ///< Set of ints (representing bond types).
			TYPE2REACTION_MAP bond_breaks_map; ///< Mapping breaking bond types to corresponding predefined bond break event.
			
			//Predefined bond forms
			PAIR_SET bond_forms_set; ///< Set of bond formable pair of ints (representing types of atoms).
			PAIR2BONDFORM_MAP bond_forms_map; ///< Mapping pair of ints to corresponding predefined bond formation event.
			INT2INT_MAP bonds_max; ///< Limiting the maximum number of bonds per atom type.
			INT2INTSMAP_MAP bondtypes_max; ///< This can be used if you want to limit the number of specific bonds that an atom type can form. This is a map to a map of ints. You give this map an int (atom type) and it returns a map. You give the returned map another int (bond type) and it returns the max number of bond types.
			
			//Predefined diffusions
			INT_SET diffusions_set; ///< Set of diffusable types.
			TYPE2DIFFUSION_MAP diffusions_map; ///<Mapping diffusable type to corresponding diffusion event.
			
			
			//Predefined depositions
			INT_SET depositions_set; ///< Set of atom types (ints) that can be parents to deposition events.
			TYPE2DEPOSITION_MAP depositions_map; ///< Mapping atom type to corresponding deposition event.
			
			//Predefined Monoatomic Desorptions
			INT_SET monodes_set; ///< Set of atom types (ints) that can be parents to monoatomic desorption events.
			TYPE2MONODES_MAP monodes_map; ///< Mapping atom type to corresponding desorption event
			
	
	};
	
	
}//end of namespace PAPRECA


#endif
