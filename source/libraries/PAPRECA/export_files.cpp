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
/// @brief Definitions for export_files.f

#include "export_files.h"

namespace PAPRECA{
	
	//---------------------------------------------General parent file Class---------------------------------------------
	//Constructors/Destructors
	File::File( ){ }
	File::File( const int &print_freq_in ) : print_freq( print_freq_in ) , is_active( true ){ };
	File::~File( ){ };
	
	void File::close( ){ file.close( ); }
	void File::setActive( ){ is_active = true; }
	void File::setIncative( ){ is_active = false; }
	const bool File::isActive( ) const{ return is_active; }
	void File::setPrintFreq( const int &print_freq_in ){ print_freq = print_freq_in; }
	const int File::getPrintFreq( ) const{ return print_freq; }
	//-------------------------------------------------End of file Class-------------------------------------------------

	//--------------------------------------------------Log File--------------------------------------------------
	//Constructors/Destructors
	Log::Log( ) : File( ){ }
	Log::~Log( ){ }
	
	//Functions
	void Log::init( ){ 
	
		file.open( "./papreca.log" );
		auto start_time = std::chrono::system_clock::now();
		time_t start_time_t = std::chrono::system_clock::to_time_t( start_time );
		
		file << "LOG FILE. PAPRECA kMC/MD Run started on " << ctime(  &start_time_t ) << " (MACHINE TIME) \n"; //Date/time
		file << "PLEASE CITE: https://doi.org/10.1016/j.commatsci.2023.112421 \n \n"; //Citations
		file << "Information about output data... \n";
		file << "For Deposition events events: site_pos (x,y,z) , rot_pos(x,y,z) , rot_theta , insertion_vel , mol_name \n";
		file << "For Bond-formation events: atom1_id , atom2_id , bond_type \n";
		file << "For Bond-breaking events: atom1_id , atom2_id , bond_type \n";
		file << "For Diffusion events: vac_pos (x,y,z) , parent_id , parent_type , insertion_vel , is_displacive , diffused_type \n";
		file << "For Monoatomic desorption events: parent_id , parent_type \n \n";
		
		
		file << std::fixed << "Step"
			<< std::setw( 14 ) << std::setprecision( 8 ) << std::fixed << "Event"
			<< std::setw( 22 ) << std::setprecision( 8 ) << std::fixed << "Time (s)" << std::endl;

	};
	
	void Log::appendDeposition( const int &KMC_loopid , const double &time , const double *site_pos , const double *rot_pos , const double &rot_theta , const double &insertion_vel , const char *mol_name ){
		
		file << std::setprecision( 8 ) << std::fixed << KMC_loopid
			<< std::setw( 20 ) << std::setprecision( 8 ) << std::fixed << "Deposition"
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << time
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << std::scientific << site_pos[0] << std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << site_pos[1] << std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << site_pos[2]
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << std::scientific << rot_pos[0] << std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << rot_pos[1] << std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << rot_pos[2]
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << std::scientific << rot_theta
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << std::scientific << insertion_vel
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << mol_name << std::endl;
	}

	void Log::appendBondForm( const int &KMC_loopid , const double &time , const LAMMPS_NS::tagint &atom1_id , const LAMMPS_NS::tagint &atom2_id , const int &bond_type){
		
		file << std::setprecision( 8 ) << std::fixed << KMC_loopid
			<< std::setw( 20 ) << std::setprecision( 8 ) << std::fixed << "Bond-form"
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << time
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << atom1_id
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << atom2_id
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << bond_type << std::endl;
	}

	void Log::appendBondBreak( const int &KMC_loopid , const double &time , const LAMMPS_NS::tagint &atom1_id , const LAMMPS_NS::tagint &atom2_id , const int &bond_type){
		
		file << std::setprecision( 8 ) << std::fixed << KMC_loopid
			<< std::setw( 20 ) << std::setprecision( 8 ) << std::fixed << "Bond-break"
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << time
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << atom1_id
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << atom2_id
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << bond_type << std::endl;
	}

	void Log::appendDiffusion( const int &KMC_loopid , const double &time , const double *vac_pos , const LAMMPS_NS::tagint &parent_id , const int &parent_type , const double &insertion_vel , const int &is_displacive , const int &diffused_type ){
		
		file << std::setprecision( 8 ) << std::fixed << KMC_loopid
			<< std::setw( 20 ) << std::setprecision( 8 ) << std::fixed << "Bond-break"
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << time
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << std::scientific << vac_pos[0] << std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << vac_pos[1] << std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << vac_pos[2]
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << parent_id
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << parent_type
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << std::scientific << insertion_vel
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << is_displacive
			<< std::setw( 5 ) << std::setprecision( 4 ) << std::fixed << diffused_type << std::endl;


	}

	//file << "For Monoatomic desorption events: parent_id , parent_type \n \n";
	void Log::appendMonoatomicDesorption( const int &KMC_loopid , const double &time , const int &parent_id , const int &parent_type ){
		
		file << std::setprecision( 8 ) << std::fixed << KMC_loopid
			<< std::setw( 20 ) << std::setprecision( 8 ) << std::fixed << "Bond-break"
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << std::scientific << time
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << parent_id
			<< std::setw( 20 ) << std::setprecision( 4 ) << std::fixed << parent_type << std::endl;
	}


	//--------------------------------------------------End of Log File --------------------------------------------------	
	
	//--------------------------------------------------HeightVTime Files--------------------------------------------------
	//Constructors/Destructors
	HeightVtime::HeightVtime( ) : File( ){ };
	HeightVtime::HeightVtime( const int &print_freq_in ) : File( print_freq_in ){ }
	HeightVtime::~HeightVtime( ){ } //In the current implementation we don't close the file in the destructor. HeightVtime files are managed by the PaprecaConfig class.
	
	//Functions
	void HeightVtime::init( ){
		
		file.open( "./heightVtime.log" );
		
		auto start_time = std::chrono::system_clock::now();
		time_t start_time_t = std::chrono::system_clock::to_time_t( start_time );
		
		file << "Height versus Time file. PAPRECA kMC/MD Run started on " << ctime(  &start_time_t ) << " (MACHINE TIME) \n"; //Date/time
		file << "PLEASE CITE: https://doi.org/10.1016/j.commatsci.2023.112421 \n\n"; //Citations
		
		file << "Time(sec)           Film Height (LAMMPS distance units) \n"; 
		
	}
	
	void HeightVtime::append( const double &time , const double &film_height ){
		
		file << std::setw( 10 ) << std::setprecision( 8 ) << std::fixed << std::scientific << time
								 << std::setw( 16 ) << std::setprecision( 8 ) << std::fixed << film_height << std::endl; 
		
		
	}
	//--------------------------------------------------End of HeightVTime Files--------------------------------------------------
	
	
	//--------------------------------------------------SurfaceCoverage Files--------------------------------------------------
	//Constructors/Destructors
	SurfaceCoverage::SurfaceCoverage( ) : File( ){ };
	SurfaceCoverage::SurfaceCoverage( const int &print_freq_in ) : File( print_freq_in ){ }
	SurfaceCoverage::~SurfaceCoverage( ){ } //In the current implementation we don't close the file in the destructor. HeightVtime files are managed by the PaprecaConfig class.
	
	//Functions
	void SurfaceCoverage::init( ){
		
		file.open( "./surface_coverage.log" );
		
		auto start_time = std::chrono::system_clock::now();
		time_t start_time_t = std::chrono::system_clock::to_time_t( start_time );
		
		file << "Surface coverage versus Time file. PAPRECA kMC/MD Run started on " << ctime(  &start_time_t ) << " (MACHINE TIME) \n"; //Date/time
		file << "PLEASE CITE: https://doi.org/10.1016/j.commatsci.2023.112421 \n\n"; //Citations
		
		file << "Time(sec)     Surface Coverage(-) \n"; 
		
	}
	
	void SurfaceCoverage::append( const double &time , const double &surface_coverage ){
		
		file << std::setw( 10 ) << std::setprecision( 8 ) << std::fixed << std::scientific << time
								 << std::setw( 12 ) << std::setprecision( 8 ) << std::fixed << surface_coverage << std::endl; 
		
		
	}
	//--------------------------------------------------End of SurfaceCoverage Files--------------------------------------------------
	
	
	
	
	
	
	//-----------------------------------------------ElementalDistribution files------------------------------------------------
	//Constructors/Destructors
	ElementalDistribution::ElementalDistribution( ) : File( ){ }	
	ElementalDistribution::ElementalDistribution( const int &print_freq_in ) : File( print_freq_in ){ }
	ElementalDistribution::~ElementalDistribution( ){	} //Again, destructor has no file.close( ), see above for explanation and for improvements
	
	//Functions
	void ElementalDistribution::init( const int &KMC_loopid , const int &types_num ){
		
		std::string file_name = "./distribution " + std::to_string( KMC_loopid ) + ".log" ;
		file.open( file_name.c_str( ) );
		
		auto start_time = std::chrono::system_clock::now();
		time_t start_time_t = std::chrono::system_clock::to_time_t( start_time );
		
		file << "PAPRECA kMC/MD Elemental Distribution file generated on " << ctime(  &start_time_t ) << " (MACHINE TIME) \n"; //Date/time
		file << "PLEASE CITE: https://doi.org/10.1016/j.commatsci.2023.112421 \n\n"; //Citations
		
		file << "Height(LAMMPS Distance Units)";
		
		for( int i = 1; i <= types_num; ++i ){ //Assumes that you use all atom types from 1 to types_num. Change accordingly for different setups
			
			file << std::setw( 12 ) << "TYPE_" << std::fixed << std::to_string( i ) << "      ";
			
			
		}
		
		file << std::setw( 12 ) << " TOTAL \n";
		
	}
	
	void ElementalDistribution::append( LAMMPS_NS::LAMMPS *lmp , double **mass_profiles , const int &types_num , const int &bins_num , const double &bin_width , double *atom_mass ){
		
		
		for( int i = 0; i < bins_num; ++i ){
			
			double height = lmp->domain->boxlo[2] + double( i ) * bin_width;
			
			int bin_total = 0;
			
			file << std::setprecision( 2 ) << std::fixed << height;
			
			for( int j = 1; j < types_num + 1; ++j ){
				
				
				int atoms_num = std::round( mass_profiles[i][j] / atom_mass[j] );
				bin_total += atoms_num;
				
				file << std::setw( 32 ) << std::fixed << atoms_num << "   ";
				
				
			}
			
			file << std::setw( 16 ) << std::fixed << bin_total << std::endl;
			
			
		}

		
		
	}
	//-----------------------------------------------End of ElementalDistribution files------------------------------------------------
	
	
	
	
	
	//---------------------------------------------------------ExecTime files----------------------------------------------------------
	//Constructors/Destructors
	ExecTime::ExecTime( ) : File( ){ };
	ExecTime::ExecTime( const int &print_freq_in ) : File( print_freq_in ){ };
	ExecTime::~ExecTime( ){ };
	
	
	//Functions
	void ExecTime::init( ){
		
		file.open( "./execTimes.log" );
		
		auto start_time = std::chrono::system_clock::now();
		time_t start_time_t = std::chrono::system_clock::to_time_t( start_time );
		
		file << "Execution Times file. PAPRECA kMC/MD Run started on " << ctime(  &start_time_t ) << " (MACHINE TIME) \n"; //Date/time
		file << "PLEASE CITE: https://doi.org/10.1016/j.commatsci.2023.112421 \n\n"; //Citations
		
		file << "Step   No.Atoms \t kMC runtime min/avg/max (sec) \t \t \t MD runtime min/avg/max (sec) \t \t \t TOTAL runtime min/avg/max (sec) \n \n"; 
		
		
	}
	
	void ExecTime::setHybridStartTimeStamp( ){
		
		t1_hybrid = MPI_Wtime( );
		
		
	}
	
	void ExecTime::calcHybridTime( const int &nprocs ){
		
		t2_hybrid = MPI_Wtime( );
		
		
		//This calculates the total elapsed time ON EACH INDIVIDUAL PROC.
		t_hybrid = t2_hybrid - t1_hybrid;
		
		MPI_Allreduce( &t_hybrid , &thybrid_min , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );
		MPI_Allreduce( &t_hybrid , &thybrid_max , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
		MPI_Allreduce( &t_hybrid , &thybrid_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
		thybrid_avg /= nprocs;
		
		thybrid_total += thybrid_avg;
		
		
	}
	
	void ExecTime::resetHybridTimeVariables( ){
		
		t_hybrid = 0.0;
		t1_hybrid = 0.0;
		t2_hybrid = 0.0;
		
		thybrid_min = 0.0;
		thybrid_max = 0.0;
		thybrid_avg = 0.0;
		
		
	}
	
	void ExecTime::setMDstartTimeStamp( ){
		
		t1_md = MPI_Wtime( );
		
		
	}
	
	void ExecTime::calcMDtime( const int &nprocs ){
		
		
		t2_md = MPI_Wtime( );
		
		t_md = t2_md - t1_md;
		
		MPI_Allreduce( &t_md , &tmd_min , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );
		MPI_Allreduce( &t_md , &tmd_max , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
		MPI_Allreduce( &t_md , &tmd_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
		tmd_avg /= nprocs;
		
		tmd_total += tmd_avg;
		
		
		
	}
	
	void ExecTime::resetMDtimeVariables( ){
		
		t_md = 0.0;
		t1_md = 0.0;
		t2_md = 0.0;
		
		tmd_min = 0.0;
		tmd_max = 0.0;
		tmd_avg = 0.0;
		
		
	}
	
	
	void ExecTime::calcKMCtime( const int &nprocs ){
		
		t_kmc = t_hybrid - t_md;
		
		MPI_Allreduce( &t_kmc , &tkmc_min , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );
		MPI_Allreduce( &t_kmc , &tkmc_max , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
		MPI_Allreduce( &t_kmc , &tkmc_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
		tkmc_avg /= nprocs;
		
		//The total KMC time is considered to be the sum of all average values. The same goes for the other time variables (hybrid, MD)
		tkmc_total += tkmc_avg;
		
	}
	
	void ExecTime::resetKMCtimeVariables( ){
		
		t_kmc = 0.0;
		
		tkmc_min = 0.0;
		tkmc_max = 0.0;
		tkmc_avg = 0.0;
		
		
		
	}
	
	void ExecTime::calcTimes( const int &nprocs ){
		
		calcMDtime( nprocs );
		calcHybridTime( nprocs );
		calcKMCtime( nprocs );
		
		
	}
	
	void ExecTime::resetTimeVariables( ){
		
		
		resetHybridTimeVariables( );
		resetMDtimeVariables( );
		resetKMCtimeVariables( );
		
		
	}
	
	void ExecTime::append( const int &step_num , const int &atoms_num ){
		
		
		file << std::setw( 0 ) << std::setprecision( 6 ) << step_num
			 << std::setw( 12 ) << std::setprecision( 8 ) << std::fixed << atoms_num
			 << std::setw( 15 ) << std::setprecision( 8 ) << std::fixed << tkmc_min << "  " << std::setw( 0 ) << std::setprecision( 8 ) << std::fixed << tkmc_avg << "  " << std::setw( 0 ) << std::setprecision( 8 ) << std::fixed << tkmc_max
			 << std::setw( 16 ) << std::setprecision( 8 ) << std::fixed << tmd_min << "  " << std::setw( 0 ) << std::setprecision( 8 ) << std::fixed << tmd_avg << "  " << std::setw( 0 ) << std::setprecision( 8 ) << std::fixed << tmd_max 
			 << std::setw( 16 ) << std::setprecision( 8 ) << std::fixed << thybrid_min << "  " << std::setw( 0 ) << std::setprecision( 8 ) << std::fixed << thybrid_avg << "  " << std::setw( 0 ) << std::setprecision( 8 ) << std::fixed << thybrid_max
			 << std::endl; 
		
		
		//Reset variables after appending
		resetTimeVariables( );
		
	}
	
	
	void ExecTime::close( ){
		
		//Write Final times before closing the file
		file << "\n \n";
		file << "Final stats: Total KMC walltime= \t \t" << tkmc_total << " sec (" << 100.0 * tkmc_total/thybrid_total << "%) \n"
		<< "\t \t \t Total MD walltime= \t \t" << tmd_total << " sec (" << 100.0 * tmd_total/thybrid_total << "%) \n"
		<< "\t \t \t Total HYBRID KMC/MD time=  " << thybrid_total << " sec" << std::endl; 
		
		
		file.close( );
		
	}
	//----------------------------------------------------End of ExecTime files-----------------------------------------------------
	
	
} //END OF NAMESPACE PAPRECA
