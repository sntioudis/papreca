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
/// @brief Declarations of export File classes.

#ifndef EXPORT_FILES_H
#define EXPORT_FILES_H

//System Headers
#include <cstdio>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>

//LAMMPS headers
#include "lammps.h"
#include "domain.h"
namespace PAPRECA{

	/* The classes below manage the result files exported by a run. New export files should be added here and a relevant implementation in the papreca_config.h and input_file.h headers should be included to read/manage/setup the files. */

	class File{
		
		/// @class PAPRECA::File
		/// @brief General parent file class. Any PAPRECA export file should be a child to this class.
		
		friend class PaprecaConfig; //PaprecaConfig reads the protected/private variables of file classes so it is declared as a friend class.
		
		public:
			//Constructors/Destructors
			File( );
			File( const int &print_freq_in );
			virtual ~File( );
			
			//Functions
			virtual void close( );
			void setActive( );
			void setIncative( );
			const bool isActive( ) const;
			void setPrintFreq( const int &print_freq_in );
			const int getPrintFreq( ) const;
			
		protected:
			std::ofstream file;
			bool is_active = false;
			int print_freq = 0;
		
	};


	class Log : public File{
	
		/// @class PAPRECA::Log
		/// @brief Child class of File, manages papreca.log files.
		
		friend class PaprecaConfig;
		
		public:
			//Constructors/Destructors
			Log( );
			~Log( );
			
			//functions
			void init( );
			void appendDeposition( const int &KMC_loopid , const double &time , const double *site_pos , const double *rot_pos , const double &rot_theta , const double &insertion_vel , const char *mol_name );
			void appendBondForm( const int &KMC_loopid , const double &time , const LAMMPS_NS::tagint &atom1_id , const LAMMPS_NS::tagint &atom2_id , const int &bond_type);
			void appendBondBreak( const int &KMC_loopid , const double &time , const LAMMPS_NS::tagint &atom1_id , const LAMMPS_NS::tagint &atom2_id , const int &bond_type);
			void appendDiffusion( const int &KMC_loopid , const double *vac_pos , const LAMMPS_NS::tagint &parent_id , const int &parent_type , const double &insertion_vel , const int &is_displacive , const int &diffused_type );
			void appendMonoatomicDesorption( const int &KMC_loopid , const double &time , const int &parent_id , const int &parent_type );
	};
	
	
	
	
	class HeightVtime : public File{
		
		/// @class PAPRECA::HeightVtime
		/// @brief Child class of File, manages heightVtime.log files.
		
		friend class PaprecaConfig;
		
		public:
			//Constructors/Destructors
			HeightVtime( );
			HeightVtime( const int &print_freq_in );
			~HeightVtime( );
			
			//functions
			void init( );
			void append( const double &time , const double &film_height );
	
	};
	
	class SurfaceCoverage : public File{
		
		/// @class PAPRECA::SurfaceCoverage
		/// @brief Child class of File, manages surface_coverage.log files.
		
		friend class PaprecaConfig;
		
		public:
			//Constructors/Destructors
			SurfaceCoverage( );
			SurfaceCoverage( const int &print_freq_in );
			~SurfaceCoverage( );
			
			//Functions
			void init( );
			void append( const double &time , const double &surface_coverage );	
		
	};
	
	class ElementalDistribution : public File{
		
		/// @class PAPRECA::ElementalDistribution
		/// @brief Child class of File, manages distribution.log files.
		
		friend class PaprecaConfig;
		
		public:
			//Constructors/Destructors
			ElementalDistribution( );
			ElementalDistribution( const int &print_freq_in );
			~ElementalDistribution( );
			
			//functions
			void init( const int &KMC_loopid , const int &types_num );
			void append( LAMMPS_NS::LAMMPS *lmp , double **mass_profiles , const int &types_num , const int &bins_num , const double &bin_width , double *atom_mass );
			
	};
	
	
	class ExecTime : public File{
		
		/// @class PAPRECA::ExecTime
		/// @brief Child class of File, manages execTimes.log files.
		
		friend class PaprecaConfig;
		
		public:
			//Constructors/Destructors
			ExecTime( );
			ExecTime( const int &print_freq_in );
			~ExecTime( );
			
			//functions
			void init( );
			void append( const int &step_num , const int &atoms_num );
			void close( ) override;
			//Total Time Calculation
			void setHybridStartTimeStamp( );
			void calcHybridTime( const int &nprocs );
			void resetHybridTimeVariables( );
			//MD time calculation
			void setMDstartTimeStamp( );
			void calcMDtime( const int &nprocs );
			void resetMDtimeVariables( );	
			//KMC time calculation
			void calcKMCtime( const int &nprocs );
			void resetKMCtimeVariables( );
			//General functions
			void calcTimes( const int &nprocs );
			void resetTimeVariables( );
			
		
		protected:
			
			double t_hybrid = 0.0 , t1_hybrid = 0.0 , t2_hybrid = 0.0 , t_md = 0.0 , t1_md = 0.0 , t2_md = 0.0 , t_kmc = 0.0;
			double thybrid_min = 0.0 , thybrid_avg = 0.0 , thybrid_max = 0.0;
			double tkmc_min = 0.0 , tkmc_avg = 0.0 , tkmc_max = 0.0;
			double tmd_min = 0.0 , tmd_avg = 0.0 , tmd_max = 0.0;
			double thybrid_total = 0.0 , tkmc_total = 0.0 , tmd_total = 0.0;
			
		
	};


}//end of namespace Papreca


#endif
