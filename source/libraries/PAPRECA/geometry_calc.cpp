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
///@brief Definitions for geometry_calc.h.

#include "geometry_calc.h"

namespace PAPRECA{
	
	//Film Height Calculation
	void calcLocalMassAndFillMassProfile( LAMMPS_NS::LAMMPS *lmp , double **mass_profiles , double &local_mass , const int &atom_type , double *atom_xyz , const double &atom_mass , const double &bin_width , const int &bins_num ){
		
		/// Appends atom_mass to the relevant bin of mass_profiles array.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in,out] mass_profiles 2D (i.e., mass_profiles[bins_num][atomtypes_num]) array containing the total mass of a specific atom type for each x-y slice (bin) in the system.
		/// @param[in,out] local_mass total mass in current MPI process.
		/// @param[in] atom_type type of the current atom.
		/// @param[in] atom_xyz array containing the coordinates of the current atom.
		/// @param[in] atom_mass mass of current atom.
		/// @param[in] bin_width width of current x-y slice (bin).
		/// @param[in] bins_num total number of x-y slices (bins).
		/// @see PAPRECA::initMassProfilesArr(), PAPRECA::deleteMassProfilesArr(), PAPRECA::fillMassProfilesTotalArrFromMassProfilesLocal(), PAPRECA::getFilmHeightFromMassBinsMethod(), PAPRECA::calcFilmHeight()
		
		int bin_num = round( ( atom_xyz[2] - lmp->domain->boxlo[2] ) / bin_width );
		mass_profiles[bin_num][atom_type] += atom_mass; //Store mass of current atom in relevant bin
		local_mass += atom_mass;		
	}
	
	double **initMassProfilesArr( const int &types_num , const int &bins_num ){
		
		/// Initializes a 2D array of mass profiles (i.e., arr[bins_num][types_num]).
		/// @param[in] types_num total number of atom types.
		/// @param[in] bins_num total number of x-y slices (bins).
		/// @return pointer to initialized 2D array of mass profiles.
		/// @see PAPRECA::calcLocalMassAndFillMassProfile(), PAPRECA::deleteMassProfilesArr(), PAPRECA::fillMassProfilesTotalArrFromMassProfilesLocal(), PAPRECA::getFilmHeightFromMassBinsMethod(), PAPRECA::calcFilmHeight()
		/// @note Arrays initialized with this function have to be deleted manually using PAPRECA::deleteMassProfilesArr().
		
		double **mass_profiles = new double*[bins_num];
		for( int i = 0; i < bins_num; ++i ){
			mass_profiles[i] = new double[types_num+1]; //types_num + 1 because in lammps 0 is not mapped to anything and we need to be able to access mass_profiles[i][types_num] (in the TCP example this is type 8);
			for( int j = 0; j < types_num + 1; ++j ){
				mass_profiles[i][j] = 0; //Initialize all elements to zero to avoid segmentation faults			
			}
				
		}	
		return mass_profiles;
	
	}
	
	void deleteMassProfilesArr( double **mass_profiles , const int &bins_num  ){
	
		/// Deletes a previously-initialized 2D array of mass profiles.
		/// @param[in,out] mass_profiles 2D array of mass profiles
		/// @param[in] bins_num number of x-y slices (bins).
		/// @see PAPRECA::initMassProfilesArr(), PAPRECA::calcLocalMassAndFillMassProfile(), PAPRECA::fillMassProfilesTotalArrFromMassProfilesLocal(), PAPRECA::getFilmHeightFromMassBinsMethod(), PAPRECA::calcFilmHeight()
		
		for( int i = 0; i < bins_num; ++i ){
				
			delete[] mass_profiles[i]; //Delete internal arrays
				
		}
			
		delete[] mass_profiles;
		
	}
	
	void fillMassProfilesTotalArrFromMassProfilesLocal( const int &bins_num , const int &types_num , double **mass_profiles_total , double **mass_profiles_local ){
	
		/// Fills global array of mass profiles on master process (i.e., proc_id==0). This is done by invoking MPI_Reduce on the local mass profile arrays.
		/// @param[in] bins_num total number of x-y slices (bins).
		/// @param[in] types_num total number of atom types.
		/// @param[in,out] mass_profiles_total global 2D mass profiles array initialized and filled ONLY for the master process (i.e., proc_id==0).
		/// @param[in] mass_profiles_local local (per MPI process) 2D array of mass profiles.
		/// @see PAPRECA::initMassProfilesArr(), PAPRECA::deleteMassProfilesArr(), PAPRECA::calcLocalMassAndFillMassProfile(), PAPRECA::getFilmHeightFromMassBinsMethod(), PAPRECA::calcFilmHeight()
		
		//Gather all mass profiles into master proc 0
		for( int i = 0; i < bins_num; ++i ){
					
			for( int j = 0; j < types_num + 1; ++j ){
					
				MPI_Reduce( &mass_profiles_local[i][j] , &mass_profiles_total[i][j] , 1 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD ); //Reduce all values in proc 0
						
			}
					
					
		}
		
		
		
	}
	
	void getFilmHeightFromMassBinsMethod( PaprecaConfig &papreca_config , LAMMPS_NS::LAMMPS *lmp , const int &proc_id , double &film_height , double **mass_profiles_total , const double &local_mass , double *atom_mass , const int &bins_num , const int &types_num , const double &bin_width ){
	
		/// Uses the global 2D mass profiles array (initialized/filled in the master process) to calculate the film height in the current step.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] proc_id id of current MPI process.
		/// @param[in,out] film_height height of current step.
		/// @param[in] mass_profiles_total global 2D array of mass_profiles (initialized/filled in master MPI process).
		/// @param[in] local_mass total mass in the current MPI process.
		/// @param[in] atom_mass LAMMPS array storing the masses of atom_types.
		/// @param[in] bins_num total number of x-y slices (bins).
		/// @param[in] types_num total number of atom types.
		/// @param[in] bin_width width of each x-y slice (bin).
		/// @see PAPRECA::initMassProfilesArr(), PAPRECA::deleteMassProfilesArr(), PAPRECA::fillMassProfilesTotalArrFromMassProfilesLocal(), PAPRECA::calcLocalMassAndFillMassProfile(), PAPRECA::calcFilmHeight()
		
		double total_mass = 0.0;
			
		//Get total mass on proc 0
		MPI_Reduce( &local_mass , &total_mass , 1 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD );

		double mass_cutoff = papreca_config.getHeightPercentage( ) * total_mass; //This is the percentage over the total mass (arbitrarily set, but you can try changing it in APRECAT config and see if this affects your system).
		double mass_cur = 0.0;
		film_height = 0.0;
			
		//Select film height on proc 0 using the user Defined
		if( proc_id == 0 ){
					
			for( int i = 0; i < bins_num; ++i ){
						
				double bin_mass = doubleArrSum( mass_profiles_total[i] , types_num + 1 );//The mass of the current bin is the sum of the mass of all the types. Here, note that we do an unnecessary loop (on types_num = 0);
				mass_cur += bin_mass;
				if( mass_cur >= mass_cutoff ){ //Mass total is non-zero ONLY on proc 0!
					film_height = lmp->domain->boxlo[2] + double( i ) * bin_width;
					break;
				}
			}			
		}
		
		//Broadcast value from master proc 0 to all other procs
		MPI_Bcast( &film_height , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
	
	}
	
	void calcFilmHeight( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const int &KMC_loopid , PaprecaConfig &papreca_config , double &film_height ){
	
		/// Calculates the film height based on the user defined film calculation method.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] KMC_loopid PAPRECA (kMC) step number.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in,out] film_height height at current PAPRECA (kMC) step.
		/// @see PAPRECA::initMassProfilesArr(), PAPRECA::deleteMassProfilesArr(), PAPRECA::fillMassProfilesTotalArrFromMassProfilesLocal(), PAPRECA::getFilmHeightFromMassBinsMethod(), PAPRECA::calcLocalMassAndFillMassProfile()
		/// @note This function also dumps ElementalDistributions files.
		
		if( papreca_config.getHeightMethod( ) != "mass_bins" && !papreca_config.getElementalDistributionsFile( ).isActive( ) ){ return; }//For now, we bin the types ONLY if we dump ElementalDistribution files OR if we calculate the film height using the mass bins method.
		
		const int natoms = *( ( int *)lammps_extract_global( lmp , "nlocal" ) );
		double **atom_xyz = ( double **)lammps_extract_atom( lmp , "x" );//extract atom positions
		double *atom_mass = ( double *)lammps_extract_atom( lmp , "mass" ); //extract atom mass: Careful, this is per-type, so you need to input the atom type number (goes from 0 to Ntypes+1)
		int *atom_types = ( int *)lammps_extract_atom( lmp , "type" );//extract atom types
		const int types_num = *( ( int * )lammps_extract_global( lmp , const_cast<char*>( "ntypes" ) ) ); //Extract total number of types

		int bins_num = round( ( lmp->domain->boxhi[2] - lmp->domain->boxlo[2] ) / papreca_config.getBinWidth( ) ) + 1; //Calculate number of bins using the dimensions of the simulation box. Add +1 to avoid segmentation faults when calculating the mass profile of the highest atom.
		double **mass_profiles = initMassProfilesArr( types_num , bins_num );
		double local_mass = 0; //Sum to hold local (proc) mass.
		
		
		for ( int i = 0; i < natoms; ++i ){
			calcLocalMassAndFillMassProfile( lmp , mass_profiles , local_mass , atom_types[i] , atom_xyz[i] , atom_mass[atom_types[i]] , papreca_config.getBinWidth( ) , bins_num );
		}
		
		
		//From the local mass bins we need to get the "full" mass bins. The full mass bins is the global mass bins which is the sum of all the local bins (on proc 0).
		double **mass_profiles_total = initMassProfilesArr( types_num , bins_num );
		fillMassProfilesTotalArrFromMassProfilesLocal( bins_num , types_num , mass_profiles_total , mass_profiles );
		
		//Calculate film height only if a method is defined (currently, only mass_bins is supported).
		if( papreca_config.getHeightMethod( ) == "mass_bins" ){ getFilmHeightFromMassBinsMethod( papreca_config , lmp , proc_id , film_height , mass_profiles_total , local_mass , atom_mass , bins_num , types_num , papreca_config.getBinWidth( ) ); }
		if( papreca_config.getElementalDistributionsFile( ).isActive( ) ){ papreca_config.dumpElementalDistributionFile( lmp , proc_id , KMC_loopid , mass_profiles_total , atom_mass , bins_num , types_num ); }
			
		deleteMassProfilesArr( mass_profiles , bins_num ); 
		deleteMassProfilesArr( mass_profiles_total , bins_num );
	}
	
	
	//Interference between atoms
	const bool atomsCollide( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config , double *atom1_xyz , const int &atom1_type , double *atom2_xyz , const int &atom2_type ){
	
		/// Checks for collisions between two atoms. It is assumed that two atoms of types i and j "collide" if the distance between them is smaller than sigma_ij.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] papreca_config object of the PAPRECA::PaprecaConfig class that stores global variables and settings for the current PAPRECA run.
		/// @param[in] atom1_xyz coordinates of the first atom (x,y, and z).
		/// @param[in] atom1_type atom type of the first atom.
		/// @param[in] atom2_xyz coordinates of the second atom (x,y, and z).
		/// @param[in] atom2_type atom type of the second atom.
		/// @return true or false if atoms collide or don't collide, respectively.
		/// @see PAPRECA::getDepoEventsFromAtom(), PAPRECA::atomHasCollisionWithMolAtoms(), PAPRECA::candidateDepoHasCollisions(), PAPRECA::getDiffEventsFromAtom()
		
		double sigma = papreca_config.getSigmaFromAtomTypes( atom1_type , atom2_type );
		if( sigma == 0 ){ 
			const std::string warn_message = "Sigma between types " + std::to_string( atom1_type ) + " and " + std::to_string( atom2_type ) + " is zero! Collisions might not be checked correctly! Please ensure that all sigmas are initialized properly. ";
			warnAll( MPI_COMM_WORLD , warn_message.c_str( ) );
		}
		double dist_sqr = get3DSqrDistWithPBC( lmp , atom1_xyz , atom2_xyz );
		if( dist_sqr < sigma * sigma ){
			return true;
			
		}
		
		return false;
		
	}
	
}//end of namespace PAPRECA
