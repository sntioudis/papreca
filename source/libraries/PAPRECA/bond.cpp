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
/// @brief Definitions for PAPRECA::Bond

#include "bond.h"


namespace PAPRECA{

	//Constructors/Destructors
	Bond::Bond( ) : bond_atom( -1 ) , bond_type( -1 ){ }; //default
	Bond::Bond( const LAMMPS_NS::tagint &bond_atom_in , const int &bond_type_in , const bool &head_parent_atom_in ) : bond_atom( bond_atom_in ) , bond_type( bond_type_in ) , head_parent_atom( head_parent_atom_in ){ }
	Bond::~Bond( ){ };
	
	//Members Functions
	const LAMMPS_NS::tagint &Bond::getBondAtom( ) const{ return bond_atom; }
	const int &Bond::getBondType( ) const{ return bond_type; }
	bool Bond::parentAtomIsHead( ) const{ return head_parent_atom; }
	void Bond::assignBondAtom( const LAMMPS_NS::tagint &bond_atom_in ){ bond_atom = bond_atom_in; }
	void Bond::assignBondType( const int &bond_type_in ){ bond_type = bond_type_in; }
	
	//Static functions
	const bool Bond::atomIDIsMapped( LAMMPS_NS::tagint &parent_atomID , PAPRECA::ATOM2BONDS_MAP &atomID2bonds ){

		/// Checks if the parent_atomID is contained in the atomID2bonds map
		/// @param[in] parent_atomID ID of parent atom.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::Bond::fillAtomID2BondsContainer(), PAPRECA::Bond::addBond2BondVector(), PAPRECA::Bond::initAtomID2BondsMap(), PAPRECA::Bond::initAndGatherBondsList()
		
		return( ( atomID2bonds.find( parent_atomID ) == atomID2bonds.end( ) ) ? false : true );
	
	}
	
	void Bond::addBond2BondVector( int &bond_type , LAMMPS_NS::tagint &parent_atomID , LAMMPS_NS::tagint &bond_atomID , const bool &head_atom_parent , ATOM2BONDS_MAP &atomID2bonds ){
	
		/// Receives the bond type and atom IDs of a bond. Then, initializes a PAPRECA::Bond and inserts it in the relevant std::vector< PAPRECA::Bond > containers. The std::vector< PAPRECA::Bond > containers are stored in the atomID2bonds map.
		/// @param[in] bond_type type of bond.
		/// @param[in] parent_atomID ID of the first atom ID of the bond.
		/// @param[in] bond_atomID ID of the second atom ID of the bond.
		/// @param[in] head_atom_parent boolean dictating whether the parent_atomID is the head atom of the current bond
		/// @param[in,out] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::Bond::fillAtomID2BondsContainer(), PAPRECA::Bond::initAtomID2BondsMap(), PAPRECA::Bond::atomIDIsMapped(), PAPRECA::Bond::initAndGatherBondsList()
		
		if ( !atomIDIsMapped( parent_atomID , atomID2bonds ) ){
			BOND_VECTOR bonds;
			Bond bond = Bond( bond_atomID , bond_type , head_atom_parent ); //initialize bond object with bond_atomID
			bonds.push_back( bond ); //insert initialized bond object in BOND_VECTOR of parent_atom
			atomID2bonds[parent_atomID] = bonds; //map parent atom to BOND_VECTOR, since the BOND_VECTOR was created in here.
		}else{ //if the parent atom is already mapped
			BOND_VECTOR &bonds = atomID2bonds[ parent_atomID ]; //retrieve unordered set from map, if map already exists
			Bond bond = Bond( bond_atomID , bond_type , head_atom_parent ); //initialize bond object with bond_atomID
			bonds.push_back( bond );
		}
	
	}
	
	void Bond::fillAtomID2BondsContainer( ATOM2BONDS_MAP &atomID2bonds , LAMMPS_NS::tagint *bonds_list , const LAMMPS_NS::bigint &bonds_num ){
		
		/// Fills atomID2bonds container with bonds gathered from LAMMPS and stored in the bonds_list array.
		/// @param[in,out] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @param[in] bonds_list array containing the IDs of bonds.
		/// @param[in] bonds_num number of bonds in the bonds_list array.
		/// @see PAPRECA::Bond::initAtomID2BondsMap(), PAPRECA::Bond::addBond2BondVector(), PAPRECA::Bond::atomIDIsMapped(), PAPRECA::Bond::initAndGatherBondsList()

		for ( int i = 0; i < bonds_num; ++i ){
			int bond_type = bonds_list[3*i];
			LAMMPS_NS::tagint bond_atom1ID = bonds_list[3*i+1];
			LAMMPS_NS::tagint bond_atom2ID = bonds_list[3*i+2];
			addBond2BondVector( bond_type , bond_atom1ID , bond_atom2ID , true , atomID2bonds ); //We consider the first atom appearing in the list (atom1) to be the head atom of this bond
																								//We use the head atom to find which processor "owns" a bonding event and avoid scanning the same pair of bonds twice( this would happen if in the future the bond formation events are discovered using a full neib list).
			addBond2BondVector( bond_type , bond_atom2ID , bond_atom1ID , false , atomID2bonds ); //repeat inverse process so the other atom also knows about the current atom.

		}
	}
	
	void Bond::initAtomID2BondsMap( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , ATOM2BONDS_MAP &atomID2bonds ){
	
		/// Initializes/fills atomID2bonds map. This is done by 1) gathering all bonds from the LAMMPS instance, 2) creating a std::vector< PAPRECA::Bond > container for each atomID, 3) filling the std::vector< PAPRECA::Bond > container with the IDs of bonded atoms.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in,out] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @see PAPRECA::Bond::fillAtomID2BondsContainer(), PAPRECA::Bond::addBond2BondVector(), PAPRECA::Bond::atomIDIsMapped(), PAPRECA::Bond::initAndGatherBondsList()

		runLammps( lmp , 0 ); //force a run 0 to initialize/refresh bonds list and initialize/refresh neighbor list.
		
		LAMMPS_NS::tagint *bonds_list = NULL;
		LAMMPS_NS::bigint bonds_num = 0;
		initAndGatherBondsList( lmp , &bonds_list , bonds_num );
		fillAtomID2BondsContainer( atomID2bonds , bonds_list , bonds_num );
		
		delete[ ] bonds_list;

	
	}
	
	
	bool Bond::atomHasBonds( const LAMMPS_NS::tagint &iatom_id , ATOM2BONDS_MAP &atomID2bonds ){
			
		return !atomID2bonds[ iatom_id ].empty( ); //atomID2bonds[ iatom_id ] returns a BONDS_vector. We use empty( ) to see if the vector is empty and return false (use of !) if that's the case because it means that the atom has no bonds.
			
	}
	
	void Bond::recursiveCollectBondedAtoms( LAMMPS_NS::tagint &atom_id , std::vector< LAMMPS_NS::tagint > &delids_local , TAGINT_SET &delids_set , ATOM2BONDS_MAP &atomID2bonds ){
		
		/// Recursively collects all bonded atoms of an atom. The IDs of collected atoms are inserted in the std::vector< tagint > container (delids_local). The vector of PAPRECA::Bond objects is retrieved for the atom_id from the atomID2bonds map. The function uses an atom tagint to start adding bonded atoms and returns when all bonded atoms are recursively collected (uses an std::unordered_set to decide that).
		/// @param[in] atom_id ID of atom.
		/// @param[in,out] delids_local vector of collected atom IDs.
		/// @param[in,out] delids_set std::unordered_set< LAMMPS_NS::tagint > containing atom IDs marked for deletion. Used to avoid duplicate deletion of atoms.
		/// @param[in] atomID2bonds PAPRECA::ATOM2BONDS_MAP container (i.e., std::unordered_map< LAMMPS_NS::tagint parent_atomID , std::vector< PAPRECA::Bond > >). The atomID2bonds container provides direct access to all the bonds of the parent atom.
		/// @note This function works consistently because all MPI processes contain the same data in their atomID2bonds map.
		/// @see PAPRECA::Bond::initAtomID2BondsMap()
		
		BOND_VECTOR &bonds = atomID2bonds[ atom_id ];
		
		
		for( const auto &bond : bonds ){
			
			LAMMPS_NS::tagint bondatom_id = bond.getBondAtom( );
			if( !elementIsInUnorderedSet( delids_set , bondatom_id ) ){
				
				delids_set.insert( bondatom_id );
				delids_local.push_back( bondatom_id );
				
				recursiveCollectBondedAtoms( bondatom_id , delids_local , delids_set , atomID2bonds ); //Call the same function again on the bondatom_id. This will collect the bonded atoms of bondatom_id
																									   //The function will stop when there are no new bonded atoms to collect
				
				
			}
			
			
		}
		
	}
	
} //End of PAPRECA namespace
