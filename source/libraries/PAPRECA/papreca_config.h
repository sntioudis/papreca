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
/// @brief Declarations for PAPRECA::PaprecaConfig class storing settings and global variables.

#ifndef PAPRECA_CONFIG_H
#define PAPRECA_CONFIG_H

//LAMMPS Headers
#include "lammps.h"
#include "library.h"
#include "random_mars.h"


//System Headers
#include <mpi.h>
#include <limits>
#include <fstream>
#include <iostream>


//kMC Headers
#include "event_list.h"
#include "lammps_wrappers.h"
#include "utilities.h"
#include "export_files.h"

namespace PAPRECA{
	
	/*This function stores variables that are necessary for the PAPRECA run.
	You can consider this class as a container storing the simulation settings and global variables.
	Most variables in this class are initialized by the PAPRECA input file.*/
	
	class PaprecaConfig{
		
		/// @class PAPRECA::PaprecaConfig
		/// @brief Class storing settings and global variables for the PAPRECA run.
		public:
		
			//Constructors/Destructors
			PaprecaConfig( );
			~PaprecaConfig( );
			
			//Essential parameters
			void setKMCsteps( const unsigned long int &KMC_steps_in );
			const unsigned long int &getKMCsteps( ) const;
			void setKMCperMD( const unsigned long int &KMC_per_MD_in );
			const unsigned long int &getKMCperMD( ) const;
			void setTimeEnd( const double &time_end_in );
			const double &getTimeEnd( );
			
			//Random numbers
			void initRanNumGenerator( LAMMPS_NS::LAMMPS *lmp , const int &seed );
			const double getUniformRanNum( ) const;
			const bool ranNumGeneratorIsInitialized( ) const;
			
			//Atom groups
			void setFluidAtomTypes( const std::vector< int > &fluid_atomtypes_in );
			const std::vector< int > &getFluidAtomTypes( ) const;
			void setFrozenAtomTypes( const std::vector< int > &frozen_atomtypes_in );
			const std::vector< int > &getFrozenAtomTypes( ) const;
			
			//Predefined events
			PredefinedReaction *getReactionFromBondType( const int &bond_type );
			PredefinedBondForm *getBondFormFromAtomTypesPair( const INT_PAIR &types_pair );
			int getMaxBondsFromSpecies( const int &atom_type );
			int getMaxBondTypesOfSpecies( const int &atom_type , const int &bond_type );
			PredefinedDiffusionHop *getDiffusionHopFromAtomType( const int &atom_type );
			PredefinedDeposition *getDepositionFromParentAtomType( const int &atom_type );
			PredefinedMonoatomicDesorption *getMonoatomicDesorptionFromAtomType( const int &atom_type );
			void initPredefinedReaction( const int &atom1_type , const int &atom2_type , const int &bond_type , const double &rate , const std::vector< int > &catalyzing_types );
			void initPredefinedBondForm( const int &atom1_type , const int &atom2_type , const int &bond_type , const double &bond_dist , const int &delete_atoms , const int &lone_candidates , const bool &same_mol , const double &rate , const std::vector< int > &catalyzing_types );
			void initPredefinedDiffusionHop( const int &parent_type , const double &insertion_vel , const double &diff_dist , const bool &is_displacive , const int &diffused_type , const double &rate , const std::string &custom_style , const std::vector< int > &custom_atomtypes );
			void initPredefinedDeposition( LAMMPS_NS::LAMMPS *lmp , const int &parent_type , const double &depo_offset , const double &insertion_vel , const std::string &adsorbate_name , const double &rate , const bool &variable_sticking , const double &sticking_coeff );
			void initPredefinedMonoatomicDesorption( const int &parent_type , const double &rate );
			void setSpeciesMaxBonds( const int &species , const int &bonds_max );
			void setSpeciesMaxBondTypes( const int &species , const int &bond_type , const int &bonds_max );
			void calcStickingCoeffs( );
			const bool predefinedCatalogHasBondBreakEvents( ) const;
			const bool predefinedCatalogHasBondFormEvents( ) const;
			const bool predefinedCatalogHasDiffusionHopEvents( ) const;
			const bool predefinedCatalogHasDepositionEvents( ) const;
			const bool predefinedCatalogHasMonoDesEvents( ) const;
			const bool predefinedCatalogIsEmpty( ) const;
			
			//Random Deposition and Diffusion Vectors
			void setRandomDepoVecs( const bool &random_depovecs_in );
			const bool &depoVecsAreRandom( ) const;
			void setRandomDiffVecs( const bool &random_diffvecs_in );
			const bool &diffVecsAreRandom( ) const;
			void setRandomDiffVecsStyle( const std::string &diffvecs_style_in );
			const std::string &getRandomDiffVecsStyle( ) const;

			
			//Deposition height settings
			void setDepoHeights( const double &height_deposcan_in , const double &height_deporeject_in );
			const double &getHeightDepoScan( ) const;
			const double &getHeightDepoReject( ) const;
			
			//Desorption settings
			void setDesorptionHeight( const double &desorb_cut_in );
			const double &getDesorptionHeight( ) const;
			void setDesorbDelMax( const int &desorb_delmax_in );
			const int &getDesorbDelMax( ) const;
			void setDesorptionStyle( const std::string &desorb_style_in );
			const std::string &getDesorptionStyle( ) const;
			
			//Film height calculation settings
			void setHeightMethod( const std::string &height_method_in );
			const std::string &getHeightMethod( ) const;
			void setHeightPercentage( const double &height_percentage_in ); //Only for method mass_bins!
			const double &getHeightPercentage( ) const; //Only for film height percentage method!
			void setBinWidth( const double &bin_width_in );
			const double &getBinWidth( ) const;
			
			//type2sigma
			void initSigmasFromLammps( LAMMPS_NS::LAMMPS *lmp );
			void setSpeciesPair2Sigma( const int &species1 , const int &species2 , const double &sigma );
			void setSigmaStyle( const std::string &sigmastyle_in );
			const std::string &getSigmaStyle( )const;
			void setSigmaMix( const std::string &mixstyle_in );
			const std::string &getSigmaMixStyle( )const;
			void mixSigmas( LAMMPS_NS::LAMMPS *lmp );
			const double getSigmaFromAtomTypes( const int &atom1_type , const int &atom2_type );
			const bool type2SigmaMapIsEmpty( ) const;
			
			//Equilibration LAMMPS
			void setMinimize1( const std::string &minimize1_in );
			const std::string &getMinimize1( ) const;
			void setMinimize2( const std::string &minimize2_in );
			const std::string &getMinimize2( ) const;
			void setTrajDuration( const int &traj_duration_in );
			const int &getTrajDuration( ) const;
			void setCtimeConvert( const double &c_convert_in );
			const double &getCtimeConvert( );
			
			//Neighbor lists
			void setNeibLists( const std::string &neiblist_half_in , const std::string &neiblist_full_in );
			const std::string &getHalfNeibListName( ) const;
			const std::string &getFullNeibListName( ) const;
			
			//Output files
			Log &getLogFile( );
			HeightVtime &getHeightVtimeFile( );
			SurfaceCoverage &getSurfaceCoverageFile( );
			void calcSurfaceCoverage( );
			ElementalDistribution &getElementalDistributionsFile( );
			ExecTime &getExecTimeFile( );
			void setupExportFiles( const int &proc_id );
			void setHybridStartTimeStamp4ExecTimeFile( const int &KMC_loopid );
			void calcHybridAndKMCTimes4ExecTimeFile( const int &nprocs , const int &KMC_loopid );
			void setMDTimeStamp4ExecTimeFile( const int &KMC_loopid );
			void calcMDTime4ExecTimeFile( const int &nprocs , const int &KMC_loopid );
			void appendExportFiles( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const double &time , const char *event_type , const double &film_height , const int &KMC_loopid );
			void dumpElementalDistributionFile( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const int &KMC_loopid , double **mass_profiles_total , double *atom_mass , const int &bins_num , const int &types_num );
			void closeExportFiles( const int &proc_id );
			void setRestartDumpFreq( const int &restart_dumpfreq_in );
			const int &getRestartDumpFreq( )const;
			void dumpLAMMPSRestart( LAMMPS_NS::LAMMPS *lmp , const int KMC_loopid );
			
		protected:
		
		
			//Random number generator.
			LAMMPS_NS::RanMars *rnum_gen = NULL; ///< Points to RanMars object as implemented in LAMMPS (see random_mars.h). The LAMMPS RanMars object is initialized from a user defined random seed in the PAPRECA input file (seed > 0 && seed < 900000000).
			
			//Essential parameters
			unsigned long int KMC_steps = 0; ///< perform that many PAPRECA steps.
			unsigned long int KMC_per_MD = 0; ///< Perform that many KMC(PAPRECA) steps for every MD(LAMMPS) step.
			double time_end = std::numeric_limits< double >::max( ); ///< Optional parameter. The user can define an ending simulation time. If not set it stays at limits of max, so you effectively never go beyond time_end;
			
			//Atom groups
			std::vector< int > fluid_atomtypes; ///< stores fluid atom types (as defined in the PAPRECA input file).
			std::vector< int > frozen_atomtypes; ///< stores frozen atom types (as defined in the PAPRECA input file).
			
			//Predefined events
			PredefinedEventsCatalog predefined_catalog; ///< stores a PAPRECA::PredefinedEventsCatalog.
			bool random_depovecs = false; ///< Controls deposition sites. If true, the deposition sites are not directly above the parent atom, but on the surface of a sphere of radius depo_offset.
			bool random_diffvecs = false; ///< Controls diffusion sites. If true, the diffusion sites are not directly above the parent atom, but on the surface of a sphere centered at the parent atom.
			std::string diffvecs_style = "3D"; ///< Type 2D means that the random diffusion vector is always ABOVE the parent type. Type 3D means that the random diffusion vector can be anywhere in 3D space. Default is 3D as it is the more general case. Note that: we do not use 2D/3D random deposition vectors as it does not really make sense to get a random vector below the parent atom type.
			double height_deposcan = -1;  ///< Scan for deposition events only +- above/below the current film height. Default at -1 which means scan everywhere.
			double height_deporeject = -1; ///< Reject deposition event above height_current + height_deporeject. Default at -1 which means do not reject anything.
			
			//Desorption settings. For thin-film growth you might wanna delete atoms flying above a certain height, after the equilibration step.
			double desorb_cut = -1; ///< Atoms above film_height + desorb_cut are deleted. The default value is -1 which means that desorption is disabled.
			int desorb_delmax = std::numeric_limits< int >::max( ); ///< Maximum number of atoms that can be deleted at once. Initialized at max limits of int so if the user does not set that, the maximum number of deleted atoms will be unlimited
			std::string desorb_style = "";	///< User defined desorption algorithm. Currently, can be gather_local or gather_all (defined in the PAPRECA input file).
			
			//Film Height Calculation Settings
			std::string height_method = ""; ///< Algorithm to calculate height. Currently, only the mass_bins method is supported.
			double height_percentage = 0.0; ///< Only used for height calculation method: mass_bins. Defines the mass percentage cutoff. The film height is calculated as the first bin whose cumulative mass is above the percentage cutoff.
			double bin_width = 1.0; ///< Bin width for height calculation or for ElementaDistribution files. Initialized at 1.0, so even if the user does not use the height_caulcation or heightbin_width command, you would be able to calculate mass bins with a width of 1.0.
			
			//Map that returns the sigma distance (i.e., equilibrium distance of species as in the LJ potential). Currently initialized from the LAMMPS input.
			INTPAIR2DOUBLE_MAP type2sigma; ///< maps types of atom types to their corresponding sigma.
			std::string sigma_style = ""; ///< method for initialization of sigma values (can be manual/LAMMPS).
			std::string sigma_mix =""; //< Two types of sigma_mix are currently supported: geom/arithm. See: https://docs.lammps.org/pair_modify.html. This variable is initialized as NONE so we know there is no mixing even when the mix keyword is not used.
			
			//Equilibrations - LAMMPS
			//The user can pass valid LAMMPS commands to minimize the configuration before the LAMMPS trajectory (minimize1) and after (minimize2)
			std::string minimize1 = ""; ///< stores a valid LAMMPS minimization command to be executed before the LAMMPS trajectory. See: https://docs.lammps.org/minimize.html.
			std::string minimize2 = ""; ///< stores a valid LAMMPS minimization command to be executed after the LAMMPS trajectory. See: https://docs.lammps.org/minimize.html.
			int traj_duration = 0; ///< This control the duration of the LAMMPS minimization. The specific LAMMPS fixes have to be defined in the LAMMPS input file.
			double c_time_convert = -1; ///< This variable is initialized from lmp->update->unit_style to allow conversion of time units to second. The constant is also pre-multiplied with the LAMMPS timestep to allow easy conversion of timesteps (i.e., traj_duration) to time interval (in seconds).
			
			//Neighbor lists
			std::string neiblist_half = ""; ///< Name of LAMMPS half neighbors list.
			std::string neiblist_full = ""; ///< Name of LAMMPS full neighbors list.
			
			//Output files
			Log log_file; ///< stores a PAPRECA::Log file The log file is always active (no need to be set to active from the input file).
			HeightVtime heightVtime_file; ///< stores a PAPRECA::HeightVtime file.
			SurfaceCoverage surfcoverage_file; ///< stores a PAPRECA::SurfaceCoverage file.
			double surface_coverage = 0.0; ///< stores a surface coverage for easier printing (if necessary).
			ElementalDistribution elementalDistribution_files; ///< stores the PAPRECA::ElementalDistribution files generated in the simulation.
			ExecTime execTime_file; ///< stores a PAPRECA::ExecTime file.
			int restart_dumpfreq = std::numeric_limits< int >::max( ); ///< dump a restart every restart_dumpfreq PAPRECA steps. Initialized at int limits, so if it is not set you virtually never dump restarts (see how restarts are dumped in lammps_wrappers.h of papreca lib).
			
	};
			
	
}//end of namespace PAPRECA


#endif
