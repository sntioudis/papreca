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
/// @brief Declarations for PAPRECA::Bond

#ifndef BOND_H
#define BOND_H

//System Headers
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <cstdio>


//LAMMPS headers
/// \cond
#include "pointers.h"
/// \endcond

//PAPRECA headers
#include "lammps_wrappers.h"


namespace PAPRECA{

	class Bond; //Forward definition of bond class
	typedef std::vector< Bond > BOND_VECTOR;
	typedef std::unordered_map< LAMMPS_NS::tagint , BOND_VECTOR > ATOM2BONDS_MAP; ///< maps atom IDs to their associated PAPRECA::BOND_VECTOR to allow easy access of bonds and bond types.
	
	class Bond{
		
		/// @class PAPRECA::Bond
		/// @brief Custom bond class (not to be confused with a LAMMPS bond).
		///
		/// Every atom in the system is associated with a vector of bond objects. The bonds vector of objects provides easy access to bonded atom IDs and bond types.
		/// PAPRECA::Bond objects enable bond formation, bond breaking, and custom diffusion events (e.g., the Fe_4PO4neib diffusion style deploys the bond vectors of the neighbors of the candidate Fe atom to determine whether or not there are 4 PO4 in the parent Fe neighborhood).
		/// The unordered_map PAPRECA::ATOM2BONDS_MAP can be used to retrieve the relevant vector of bond objects using a LAMMPS atom ID.
		

		private:
			LAMMPS_NS::tagint bond_atom = -1;
			int bond_type = -1;
			bool head_parent_atom = false;
			
		public:
			//Constructors/Destructors
			Bond( ); //Default
			Bond( const LAMMPS_NS::tagint &bond_atom_in , const int &bond_type_in , const bool &head_parent_atom_in );	
			~Bond( );
			
			//Member functions
			const LAMMPS_NS::tagint &getBondAtom( ) const;
			const int &getBondType( ) const;
			bool parentAtomIsHead( ) const; ///< For each bond in the system a tail and a head atom are (randomly) defined. This prevents the discovery of identical (i.e., defined twice) PAPRECA::BondBreak events.
			void assignBondAtom( const LAMMPS_NS::tagint &bond_atom_in );
			void assignBondType( const int &bond_type_in );
			
			//Static functions
			static const bool atomIDIsMapped( LAMMPS_NS::tagint &parent_atomID , ATOM2BONDS_MAP &atomID2bonds );
			static void addBond2BondVector( int &bond_type , LAMMPS_NS::tagint &parent_atomID , LAMMPS_NS::tagint &bond_atomID , const bool &head_atom_parent , ATOM2BONDS_MAP &atomID2bonds );
			static void fillAtomID2BondsContainer( ATOM2BONDS_MAP &atomID2bonds , LAMMPS_NS::tagint *bonds_list , const LAMMPS_NS::bigint &bonds_num );
			static void initAtomID2BondsMap( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , ATOM2BONDS_MAP &atomID2bonds );
			static bool atomHasBonds( const LAMMPS_NS::tagint &iatom_id , ATOM2BONDS_MAP &atomID2bonds );
			static void recursiveCollectBondedAtoms( LAMMPS_NS::tagint &atom_id , std::vector< LAMMPS_NS::tagint > &delids_local , TAGINT_SET &delids_set , ATOM2BONDS_MAP &atomID2bonds );
	};
	
} //End of PAPRECA namespace

#endif
