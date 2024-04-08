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
///@brief Definitions for debug.h

//PAPRECA headers
#include "debug.h"



namespace PAPRECA{
	
	
	void debugPrintBondMapPairs( ATOM2BONDS_MAP const &bonds_map , const int &proc_id ){
		if ( proc_id == 0 ){
			for ( auto const &pair: bonds_map ) {
				int parent_id = pair.first;
				
				std::cout << "This is the bond list of atom with id: " << parent_id << " on proc " << proc_id << std::endl;
				std::cout << "~~~this bond list has " << pair.second.size( ) << " member(s) \n";
				for ( auto const &bond : pair.second ){
					std::cout <<  "			Atom " << pair.first << " HEADATOM1(?)= " << bond.parentAtomIsHead( ) << " is paired with atom " << bond.getBondAtom( ) << std::endl;
				}
			}
		}
	}
	
	void debugPrintBasicAtomInfo( LAMMPS_NS::LAMMPS *lmp , const int &proc_id ){
		
		int atoms_num = lmp->atom->nlocal;
		LAMMPS_NS::tagint *atom_ids = lmp->atom->tag;
		int *atom_types = lmp->atom->type;
		double **atom_xyz = ( double **)lammps_extract_atom( lmp , "x" );//extract atom positions
		double *atom_mass = ( double *)lammps_extract_atom( lmp , "mass" ); //extract atom mass
		LAMMPS_NS::tagint *molid = ( LAMMPS_NS::tagint *)lammps_extract_atom( lmp , "molecule" ); //Extract molecule id

		for ( int i = 0 ; i < atoms_num ; ++i ){
			
			if ( molid[i] == 1 ){
				printf( "This is atom with id %d of type %d and mass %f on proc %d at pos (%f,%f,%f) \n" , atom_ids[i] , atom_types[i] , atom_mass[atom_types[i]] , proc_id , atom_xyz[i][0] , atom_xyz[i][1] , atom_xyz[i][2] );
				std::cout << "The present atom belongs to molecule " << molid[i] << std::endl;
				printf( "\n \n" );
			}
			
		}
		
		
	}
	
	void debugPrintNeighborLists( LAMMPS_NS::LAMMPS *lmp , const int &proc_id ){
		
		int *atom_ids = lmp->atom->tag;
		int neiblist_id = lammps_find_pair_neighlist( lmp , "zero" , 1 , 0 , 0 );
		int inum = lammps_neighlist_num_elements( lmp , neiblist_id );
		int iatom , jneib;
		int numneigh;
		int *neighbors = NULL;
		for ( int i = 0 ; i < inum ; ++i ){ 
		
			lammps_neighlist_element_neighbors( lmp , neiblist_id , i , &iatom , &numneigh , &neighbors );
			
			if ( atom_ids[iatom] == 1087 || atom_ids[i] == 1087 ){ //Uncomment this if you wanna check a specific pair
				std::cout << "ATOM ID: " << atom_ids[iatom] << " on proc " << proc_id << std::endl;
				
				for ( int j = 0 ; j < numneigh ; ++j ){
				
					jneib = neighbors[j];
					jneib &= NEIGHMASK;
					
					std::cout << "		is neighbors with Atom id: " << atom_ids[jneib] << " on proc " << proc_id << std::endl;
					
				}
			}
			
		}

		
		
	}
	
	void debugCheckBondsInNeibLists( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , ATOM2BONDS_MAP &atomID2bonds){
		
		int *atom_ids = lmp->atom->tag;
		LAMMPS_NS::tagint **bond_ids = lmp->atom->bond_atom;
		int neiblist_id = lammps_find_pair_neighlist( lmp , "zero" , 1 , 0 , 0 );
		int inum = lammps_neighlist_num_elements( lmp , neiblist_id );
		int iatom , jneib;
		int numneigh;
		int *neighbors = NULL;
		
		std::cout << inum << std::endl;
		
		for ( int i = 0 ; i < inum ; ++i ){ 
		
			
			lammps_neighlist_element_neighbors( lmp , neiblist_id , i , &iatom , &numneigh , &neighbors );
				
		}
		std::cout << std::endl;
		
		
	}
	
	void debugPrintBondsList( LAMMPS_NS::tagint *bonds_list , LAMMPS_NS::bigint &bonds_num , const int &proc_id ){
		
		if ( proc_id == 0 ){
			for ( int i = 0; i < bonds_num; ++i ){
			
				printf( "bond % 4d: type= %d, atoms: % 4d % 4d\n" , i , bonds_list[3*i] , bonds_list[3*i+1] , bonds_list[3*i+2] );
			
			}
		}
	
	}
	
	void debugPrintType2SigmaMap( INTPAIR2DOUBLE_MAP &types2sigma ){
		
		for( const auto &element : types2sigma ){
			
			INT_PAIR pair = element.first;
			INT_PAIR pair_inverse( pair.second , pair.first );
			double sigma = element.second;
			printf( "types %d and %d have a sigma of %f \n" , pair.first , pair.second , sigma );
			
		}
		
	}
	
	void debugPrintEventInfo( Event *event , const int &proc_id ){
		
		
		if( event->getType( ) == "RXN-BREAK" ){
			BondBreak *bond_break = dynamic_cast<BondBreak*>( event );
			std::cout << "This is a bond breaking event between " << bond_break->getAtom1ID( ) << " and " << bond_break->getAtom2ID( ) << " of bond type " << bond_break->getBondType( ) << " with rate " << bond_break->getRate( ) << " on proc " << proc_id << std::endl;
		}else if( event->getType( ) == "RXN-FORM" ){
			BondForm *bond_form = dynamic_cast<BondForm*>( event );
			std::cout << "This is a bond forming event between " << bond_form->getAtom1ID( ) << " and " << bond_form->getAtom2ID( ) << " of bond type " << bond_form->getBondType( ) << " with rate " << bond_form->getRate( ) << " on proc " << proc_id << std::endl;
		}else if( event->getType( ) == "DEPO" ){
			Deposition *deposition = dynamic_cast<Deposition*>( event );
		}else if( event->getType( ) == "DIFF" ){
			Diffusion *diffusion = dynamic_cast<Diffusion*>( event );
		
		}
		
	
	}
	
} //End of PAPRECA namespace
