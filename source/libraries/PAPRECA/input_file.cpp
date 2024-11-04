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
/// @brief  Definitions for input_file.h

#include "input_file.h"

namespace PAPRECA{
	
	//SUPPLEMENTARY SETTER FUNCTIONS FOR APRECAT CONFIG
	void setTimeUnitsConversionConstant( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config ){
		
		/// Retrieves the LAMMPS units style from the provided LAMMPS object (lmp) and sets the conversion time variable of the PAPRECA::PaprecaConfig object.
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @note This operation is necessary, since PAPRECA will always report time in seconds (SI units), regardless of the LAMMPS time units. Any other variable in PAPRECA (e.g., length) has units consistent with LAMMPS units.
		
		double c_time_convert = -1; //Initialized at -1 to allow illegal run detection (e.g., if c_time_convert is negative).
		
		const char *style = lmp->update->unit_style;
		const double dt = lmp->update->dt;
		
		if( dt == 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "LAMMPS timestep was set to 0 in the LAMMPS input file." ); }
		
		if( strcmp( style , "lj" ) == 0 ){
			warnAll( MPI_COMM_WORLD , "Using LJ units does not allow direct conversion of time to seconds (or its derivatives). The time reported by PAPRECA will not account for the MD time intervals (i.e., PAPRECA time will only be KMC time).");
			c_time_convert = 0; //In that case we set the conversion constant to 0. When the conversion constant is multiplied by the number of timesteps, the resulting time interval will be zero.
			//We do that because LJ units cannot be directly converted to seconds (or its derived units).
			//In case the user decides to use LJ units, the time reported by PAPRECA will be pure KMC time (i.e., MD time intervals will be neglected).
		}else if( strcmp( style , "real" ) == 0 ){
			
			c_time_convert = 1.0e-15 * dt;
			
		}else if( strcmp( style , "metal" ) == 0 ){
			
			c_time_convert = 1.0e-12 * dt;
			
		}else if( strcmp( style , "si" ) == 0 ){
			
			c_time_convert = dt;
			
		}else if( strcmp( style , "cgs" ) == 0 ){
			
			c_time_convert = dt;
			
		}else if( strcmp( style , "electron" ) == 0 ){
			
			c_time_convert = 1.0e-15 * dt;
			
		}else if( strcmp( style , "micro" ) == 0 ){
			
			c_time_convert = 1.0e-6 * dt;
			
		}else if( strcmp( style , "nano" ) == 0 ){
			
			c_time_convert = 1.0e-9 * dt;
			
		}
		
		if( c_time_convert < 0.0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Configuration parameter c_time_convert was not initialized properly from LAMMPS unit style (input_file.cpp)." ); }
		
		papreca_config.setCtimeConvert( c_time_convert );	
		
	}
	
	//ACCEPTABLE OPTIONAL COMMAND KEYWORDS and KEYWORD IDENTIFICATION FUNCTIONS.
	void checkForAcceptableKeywordsUsedMultipleTimes( std::vector< std::string > &commands , const std::string &keyword ){
		
		/// Receives the commands vector of strings and aborts if the provided "keyword" is used multiple times.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in] keyword specific keyword to be checked.
		
		const int pos1 = getStringPosInStringVec( keyword , commands );
		
		if( pos1 != -1 ){
			
			std::vector< std::string > sub_commands = getSubVectorFromVector( commands , pos1 + 1 , commands.size( ) ); //Get a smaller vector from the position that we first encountered the command until the end of the command line.
			const int pos2 = getStringPosInStringVec( keyword , sub_commands );
			if( pos2 != -1 ){ //If the same keyword is encountered again, abort!
				allAbortWithMessage( MPI_COMM_WORLD , "Keyword: " + keyword + " used multiple times in command: " + commands[0] + "." );
			}
		}
		
	}
	
	void check4AcceptableKeywords( std::vector< std::string > &commands , const int &start , std::unordered_set< std::string > &acceptable_keywords , const bool &accept_bool ){
		
		/// Checks if a specific command line includes acceptable keyword. The "acceptable_keyword" unordered set has to be initialized correctly by the user, before the call to this function. The present function aborts if it detected a non-acceptable keyword (i.e, a keyword not included in the acceptable_keywords set).
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in] start current index position in the commands vector. Only check keywords for indices that are equal or greater than start.
		/// @param[in] acceptable_keywords unordered_set containing the acceptable (for this command) keywords. Has to be initialized by the user before the call to this function.
		/// @param[in] accept_bool if true, then yes and no keywords are also accepted in the command line.
		
		//If the accept_bool option is active (i.e., true), yes/no keywords are also accepted in the command line
		if( acceptable_keywords.size( ) == 0 ){
			
			allAbortWithMessage( MPI_COMM_WORLD , "Acceptable keywords unordered set was not initialized for command " + commands[0] + " with optional keywords in input_file.cpp." );
			
		}
		
		if( start < commands.size( ) ){
			
			//When checking for optional commands, the first encountered string HAS to be a keyword (and not a number or a bool).
			if( stringIsNumber( commands[start] ) || stringIsBool( commands[start] ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Illegal command line:" + getConcatenatedStringWithSpacesFromStringsVector( commands , 0 , commands.size( ) ) ); }
			
		}
		
		for( int i = start; i < commands.size( ); ++i ){
		
			if( ( accept_bool && !stringIsNumber( commands[i] ) && !stringIsBool( commands[i] ) ) || ( !accept_bool && !stringIsNumber( commands[i] ) ) ){

					
					if( !elementIsInUnorderedSet( acceptable_keywords , commands[i] ) ){
					
						allAbortWithMessage( MPI_COMM_WORLD , "Non-acceptable keyword: " + commands[i] + " for command " + commands[0] + "." );
					}else{
						
						checkForAcceptableKeywordsUsedMultipleTimes( commands , commands[i] ); //If the keyword is acceptable we still have to check if it is defined many times in the same command line
						
					}
					
					
			}		
			

		}
		
		
	}
	
	
	void processCatalyzedOption( std::vector< std::string > &commands , int &current_pos , std::vector< int > &catalyzing_types ){
		
		/// Checks if the optional "catalyzed" keyword is present in a command. If yes then the catalyzing types are inserted in the catalyzing_types std::vector< int > container.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] current_pos position on the commands vector. This variable gets updated after the call to the function to denote the exiting position in the commands vector.
		/// @param[in,out] catalyzing_types vector containing the catalyzing types for a specific command (typically create_BondBreak or create_BondForm).
		
		checkForAcceptableKeywordsUsedMultipleTimes( commands, "catalyzed" ); //This has to be done for all optional commands to ensure that no multiple-defined commands exist.
			
		if( commands.size( ) == current_pos + 1 ){ //Means that the catalyzed keyword is present but there are no other inputs after the keyword.
				
			allAbortWithMessage( MPI_COMM_WORLD , "Illegal catalyzed keyword in " + commands[0] + " command. Has to be catalyzed N type1 type2 ... typeN." );
				
		}else{
			const int catalyzed_num = string2Int( commands[current_pos + 1] );
			if( catalyzed_num <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Catalyzed num cannot be 0 or negative for catalyzed optional keyword." ); }
			
			const int catalyzed_start = current_pos + 2; //Because catalyzed is in option_pos and catalyzed_num is in option_pos+1
			const int catalyzed_end = catalyzed_start + catalyzed_num;
			if( commands.size( ) < catalyzed_end ){
				allAbortWithMessage( MPI_COMM_WORLD , "Illegal catalyzed keyword in " + commands[0] + " command. Has to be catalyzed N type1 type2 ... typeN." );
			}else{
				for( int i = catalyzed_start; i < catalyzed_end; ++i ){
					
					if( string2Int( commands[i] ) < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Negative atom type detected (" + commands[i] + ") during the initialization of catalyzing types for + " + commands[0] + " command." ); }
					
					if( elementIsInVector( catalyzing_types , string2Int( commands[i] ) ) ){
						warnAll( MPI_COMM_WORLD , "Catalyzing type : " + commands[i] + " defined multiple times for " + commands[0] + " command. The catalyzing type will only be inserted once in the catalyzing vector." );
					}else{
						catalyzing_types.push_back( string2Int( commands[i] ) );
					}
				}
				
				current_pos = catalyzed_end; //If you've reached this point there were no aborts and the function will normally update the current_pos variable
			}
				
		}
			
		
	}
	
	void processBondLimitOption( std::vector< std::string > &commands , int &current_pos , double &length_equil , double &length_perc ){
		
		length_equil = string2Double( commands[current_pos+1] );
		length_perc = string2Double( commands[current_pos+2] );
		
		if( length_perc >= 1.0 || length_perc <= 0.0 ){ allAbortWithMessage( MPI_COMM_WORLD , "length percentage for bonding events has to be between 0.0 and 1.0 (exclusive on both ends). Check input file!"); }
		
		current_pos += 3;
	}
	
	void processSigmaMixOptions( std::vector< std::string > &commands , PaprecaConfig &papreca_config , int &current_pos ){
		
		/// Sets the sigma mix variable of the PAPRECA::PaprecaConfig object (based on the input command-line).
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @param[in,out] current_pos position on the commands vector. This variable gets updated after the call to the function to denote the exiting position in the commands vector.
		
		checkForAcceptableKeywordsUsedMultipleTimes( commands, "mix" ); //This has to be done for all optional commands to ensure that no multiple-defined commands exist.

		std::string style = commands[current_pos+1];
		
		if( commands.size( ) < current_pos+2 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid mix keyword for " + commands[0] + " command." ); }
		if( style == "geom" ){
			papreca_config.setSigmaMix( style );
		}else if( style == "arithm" ){
			papreca_config.setSigmaMix( style );
					
		}else if( style =="no" ){
			papreca_config.setSigmaMix( style );
		}else{
			allAbortWithMessage( MPI_COMM_WORLD , "Illegal mixing options for " + commands[0] + " command. Has to be mix geom or mix arithm or mix no" );
		}
			
		current_pos += 2;
		
	}
	
	void processBinWidthOptionForElementalDistributions( std::vector< std::string > &commands , PaprecaConfig &papreca_config , int &current_pos ){
		
		/// Sets the bind width variable of the PAPRECA::PaprecaConfig object (based on the input command-line).
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @param[in,out] current_pos position on the commands vector. This variable gets updated after the call to the function to denote the exiting position in the commands vector.
		
		checkForAcceptableKeywordsUsedMultipleTimes( commands , "bin_width" );
		
		if( commands.size( ) < current_pos + 2 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid option bin_width option in " + commands[0] + " command." ); }
		
		double bin_width = string2Double( commands[current_pos+1] );
		if( bin_width <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Attempted to set negative bin_width in " + commands[0] + " command." ); }
		
		if( papreca_config.getHeightMethod( ) == "mass_bins" ){ warnAll( MPI_COMM_WORLD , "Bin width was already set in a previous height_calculation command. The value you've entered for " + commands[0] + " command will overwrite the previous bin_width value." ); }
		papreca_config.setBinWidth( bin_width );
		
		current_pos += 2;
		
		
		
	}
	
	double getStickingCoeffFromDepositionEventOptions( std::vector< std::string > &commands , int &current_pos ){
		
		/// Extracts the sticking coefficient from the commands vector returns it to the caller function.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] current_pos position on the commands vector. This variable gets updated after the call to the function to denote the exiting position in the commands vector.
		/// @return sticking coefficient (returns -1 if the sticking coefficient is variable).
		
		std::string error_message = "Illegal sticking_coeff keyword in " + commands[0] + " command. Has to be sticking_coeff variable/constant C (where C is the sticking coeff to be in the command ONLY if the constant option is used).";
		
			
		if( commands.size( ) < current_pos + 2 ){
				
			allAbortWithMessage( MPI_COMM_WORLD , error_message );
		}else{
				
			std::string coeff_type = commands[current_pos+1];
			if( coeff_type == "constant" ){
				if( commands.size( ) < current_pos+3 ){
					allAbortWithMessage( MPI_COMM_WORLD , error_message );
				}
				double sticking_coeff = string2Double( commands[current_pos+2] );
				if( sticking_coeff > 1.0 || sticking_coeff <= 0.0 ){
					allAbortWithMessage( MPI_COMM_WORLD , "Illegal sticking_coeff option in " + commands[0] + " command. The sticking coefficient has to be a double number between 0 and 1" );
				}
				current_pos += 3;
				return sticking_coeff;
			}else if( coeff_type == "variable" ){
					
				if( commands.size( ) < current_pos + 2 ){
					allAbortWithMessage( MPI_COMM_WORLD , error_message );
				}
				current_pos += 2;
				return -1.0;
			}
				
		}

		
		return -1.0;
		
		
	}
	
	
	double getRateFromInputRateOptions( std::vector< std::string > &commands , int &current_pos ){
		
		/// Calculates the predefined event rate by calling the chosen (in the command-line) function in rates_calc.h and returns it to the caller function.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] current_pos position on the commands vector. This variable gets updated after the call to the function to denote the exiting position in the commands vector.
		/// @return rate of predefined event.
		
		//Function that checks if the calc_arrhenius option is active in many commands and returns a rate accordingly.
		double rate = -1;
		
		if( commands[current_pos] == "rate_arrhenius" ){
			
			if( commands.size( ) < current_pos + 4 ){
				allAbortWithMessage( MPI_COMM_WORLD , "Incorrect number of inputs after rate_arrehnius option, must be rate_arrhenius activation_energy attempt_freq temperature." );
			}else{
				double activation_energy = string2Double( commands[current_pos+1] );
				double attempt_freq = string2Double( commands[current_pos+2] );
				double temperature = string2Double( commands[current_pos+3] );
				if( activation_energy < 0 || attempt_freq < 0 || temperature < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "activation_energy, attempt_freq, and temperature have to be non-negative numbers in " + commands[0] + " command." ); }
				
				rate = getRateFromArrhenius( activation_energy , attempt_freq , temperature );
				current_pos += 4;
			}
			
		}else if( commands[current_pos] == "rate_manual" ){
			
			if( commands.size( ) < current_pos + 2 ){
				allAbortWithMessage( MPI_COMM_WORLD , "Incorrect number of inputs after rate_manual option, must be rate_manual RATE." ); 
			}else{
				rate = string2Double( commands[current_pos+1] );
				current_pos += 2;
			}
			
		}else if( commands[current_pos] == "rate_hertz" ){
			
			if( commands.size( ) < current_pos + 5 ){
				allAbortWithMessage( MPI_COMM_WORLD , "Incorrect number of inputs after rate_hertz option, must be rate_hertz pressure ads_area ads_mass temperature." ); 
			}else{
				double pressure = string2Double( commands[current_pos+1] );
				double ads_area = string2Double( commands[current_pos+2] );
				double ads_mass = string2Double( commands[current_pos+3] );
				double temperature = string2Double( commands[current_pos+4] );
				if( pressure < 0 || ads_area < 0 || ads_mass < 0 || temperature < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "pressure, ads_area, ads_mass, and temperature have to be non-negative in " + commands[0] + " command."); }
				
				rate = getDepoRateFromHertzKnudsen( pressure , ads_area , ads_mass , temperature );
				current_pos += 5;
			}
		}else{
			
			allAbortWithMessage( MPI_COMM_WORLD , "Unknown rate option: " + commands[current_pos] + "." );
		}
		
		if( rate <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Calculation of rate resulted in non-positive rate:" + std::to_string( rate ) + " for " + commands[0] + " command." ); }
		
		return rate; //You should never return a rate of -1 (the code aborts if a rate is not defined within the if blocks).
		

	}
	
	void processCustomDiffEventOptions( std::vector< std::string > &commands , int &current_pos , std::string &custom_style , std::vector< int > &style_atomtypes ){
		
		/// Initializes custom_style and style_atomtypes for custom diffusion events.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] current_pos position on the commands vector. This variable gets updated after the call to the function to denote the exiting position in the commands vector.
		/// @param[in,out] custom_style type of custom style.
		/// @param[in,out] style_atomtypes potentially useful auxiliary vector storing atom types and passing them to the main function of papreca.cpp.
		
		checkForAcceptableKeywordsUsedMultipleTimes( commands, "custom" ); //This has to be done for all optional commands to ensure that no multiple-defined commands exist.
		std::string error_message = "Illegal custom diffusion keyword. Has to be custom STYLE N type1 type2 ... typeN (where N is the style_types num) If there are no style types the option should be custom STYLE 0. Currently only Fe_4PO4neib is supported as a custom diffusion style.";
		if( commands.size( ) < current_pos + 3 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		std::string style = commands[current_pos+1];
		if( style == "Fe_4PO4neib" ){
			custom_style = style;
			if( string2Int( commands[current_pos+2] ) != 1 ){ allAbortWithMessage( MPI_COMM_WORLD , "Custom diffusion style Fe_4PO4 only works with 1 custom atomtype is defined (i.e., the P type)." ); }
		}else{
			allAbortWithMessage( MPI_COMM_WORLD , error_message );
		}
		
		int types_num = string2Int( commands[current_pos+2] );
		if( types_num <= 0 ){ return; } //Exit if there are no style types
		
		
		int types_start = current_pos + 3;
		int types_end = types_start + types_num;
		
		if( commands.size( ) < types_end ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
			
		for( int i = types_start; i < types_end; ++i ){
			
			if( string2Int( commands[i] ) < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Negative atom type:" + commands[i] + " detected in " + commands[0] + " command." ); }
			
			if( elementIsInVector( style_atomtypes , string2Int( commands[i] ) ) ){
				warnAll( MPI_COMM_WORLD , "Style atom type : " + commands[i] + " defined multiple times for " + commands[0] + " command. The style atom type will only be inserted once in the style_atomtypes vector." );
			}else{
				style_atomtypes.push_back( string2Int( commands[i] ) );
			}
				
		}
		
		current_pos = types_end;
		
	
		
	}
	
	
	
	
	//ACCEPTABLE INPUT COMMANDS
	void executeKMCstepsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets the total number of KMC steps in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid KMC_steps command in PAPRECA input file. Correct formatting: KMC_Steps N (where N is the number of KMC steps)." ); }

		unsigned long int steps = string2UnsignedLongInt( commands[1] );
		if( steps <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Cannot start a PAPRECA simulation with a non-positive number of KMC steps" ); }
		
		papreca_config.setKMCsteps( steps );
		
	}
	
	void executeKMCperMDCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets the frequency of MD (LAMMPS) steps (i.e., how many KMC stages are performed for every MD stages) in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid KMC_per_MD command in PAPRECA input file. Correct formatting: KMC_per_MD N (where N is the frequency: N KMC steps per 1 MD step )." ); }
		
		unsigned long int KMC_per_MD = string2UnsignedLongInt( commands[1] );
		if( KMC_per_MD <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "KMC per MD has to be non-negative." ); }
		
		papreca_config.setKMCperMD( KMC_per_MD );
		
		
		
	}
	
	void executeTimeEndCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets the ending target time in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
	
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid time_end command in PAPRECA input file. Correct formatting: time_end N (where N is the target end time of the simulation." ); }
		
		double time_end = string2Double( commands[1] );
		if( time_end <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "The target ending simulation time has to be a positive (double) number." ); }
		
		papreca_config.setTimeEnd( time_end );
	
	
	}
	
	void executeRandomSeedCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Initializes the LAMMPS random number generator (RanMars) of the PAPRECA::PaprecaConfig object using the input random seed.
		/// @param[in] lmp pointer to LAMMPS object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
	
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid random_seed command in PAPRECA input file. Correct formatting: random_seed N (where N is is a random seed seed > 0 && seed < 900000000)." ); }
		
		
		int random_seed = string2Int( commands[1] );
		
		//Check that obtained integer is within the valid range for the random number generator. The range is defined by the random number generator of LAMMPS ( see RanMars.h header).
		if( random_seed <= 0 || random_seed >= 900000000 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid random_seed command in PAPRECA input file. Correct formatting: random_seed N (where N is is a random seed seed > 0 && seed < 900000000)." ); }
		
		papreca_config.initRanNumGenerator( lmp , random_seed );
		
		
	}
	
	void executeFluidAtomTypesCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Initializes the fluid atom types in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		std::string error_message = "Invalid fluid_atomtypes command in PAPRECA input file. Correct formatting: fluid_atomtypes N type1 type2 ... typeN (where N is the number of atomtypes).";
			
		std::vector< int > fluid_atomtypes;
		int types_num = string2Int( commands[1] );
		if( types_num <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "types num cannot be 0 or negative for fluid_atomtypes command." ); }
			
		int types_start = 2;
		int types_end = types_start + types_num;
			
		if( commands.size( ) != types_end ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
			
		for( int i = types_start; i < types_end; ++i ){
			
			if( string2Int( commands[i] ) < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Negative atom type number: " + commands[i] + " detected in " + commands[0] + " command." ); }
			if( elementIsInVector( fluid_atomtypes , string2Int( commands[i] ) ) ){
				warnAll( MPI_COMM_WORLD , "Fluid atom type : " + commands[i] + " defined multiple times for " + commands[0] + " command. The fluid atom type will only be inserted once in the fluid_atomtypes vector." );
			}else{
				fluid_atomtypes.push_back( string2Int( commands[i] ) );
			}
				
		}
			
		papreca_config.setFluidAtomTypes( fluid_atomtypes );
		
	}
	
	
	void executeFrozenAtomTypesCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Initializes the frozen atom types in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid frozen_atomtypes command in PAPRECA input file. Correct formatting: frozen_atomtypes N type1 type2 ... typeN (where N is the number of atomtypes).";
			
		std::vector< int > frozen_atomtypes;
		int types_num = string2Int( commands[1] );
		if( types_num <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "types num cannot be 0 or negative for frozen_atomtypes command." ); }
			
		int types_start = 2;
		int types_end = types_start + types_num;
			
		if( commands.size( ) != types_end ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
			
		for( int i = types_start; i < types_end; ++i ){
				
			if( string2Int( commands[i] ) < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Negative atom type number: " + commands[i] + " detected in " + commands[0] + " command." ); }
				
			if( elementIsInVector( frozen_atomtypes , string2Int( commands[i] ) ) ){
				warnAll( MPI_COMM_WORLD , "Frozen atom type : " + commands[i] + " defined multiple times for " + commands[0] + " command. The frozen atom type will only be inserted once in the frozen_atomtypes vector." );
			}else{
				frozen_atomtypes.push_back( string2Int( commands[i] ) );
			}
				
		}
			
		papreca_config.setFrozenAtomTypes( frozen_atomtypes );
		
	}
	
	
	
	void executeDesorptionCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
	
		/// Sets the necessary parameters for the desorption command in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid desorption command. Must be desorption N style (where N is a double number denoting the desorption height). Style can be gather_all or gather_local or LAMMPS_region. Acceptable keyword ONLY for the gather_all and gather_local styles: max N (where N is the maximum number of atoms that can be deleted at once)";
		if( commands.size( ) < 3 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		
		if( papreca_config.getHeightMethod( ).empty( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Cannot setup desorption before defining a height calculation method." ); }
		
		double desorption_height = string2Double( commands[1] );
		if( desorption_height <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "The desorption height in " + commands[0] + " command has to be a positive (double) number." ); }
		
		std::string style = commands[2];
		
		if( style == "gather_all" && style == "gather_local" && style != "LAMMPS_region" ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); } //If you reach any point below you are 100% certain that the correct styles are used
		if( style == "LAMMPS_region" && commands.size( ) > 3 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); } //Exit immediately if LAMMPS region is called with additional keywords (because LAMMPS_region has no extra arguments
		
		papreca_config.setDesorptionHeight( desorption_height );
		papreca_config.setDesorptionStyle( style );
		
		//Check for the acceptable argument "max" for styles gather_local and gather_all
		if( ( style == "gather_all" || style == "gather_local" ) && commands.size( ) > 3 ){
			
			if( commands[3] != "max" || commands.size( ) != 5 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
			
			int desorb_delmax = string2Int( commands[4] );
			papreca_config.setDesorbDelMax( desorb_delmax );
			
		}
	
	
	}
	
	void executeHeightCalculationCommand( std::vector< std:: string > &commands , PaprecaConfig &papreca_config ){
	
		/// Sets all parameters relevant to height calculation for the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( commands.size( ) != 4 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid height_calculation command. Must be height_calculation METHOD settings. Currently only one method is supported (mass_bins). Acceptable command: height_calculation mass_bins cutoff_percentage bin_width." ); }
	
		if( commands[1] == "mass_bins" ){
			
			papreca_config.setHeightMethod( commands[1] );
			
			double cutoff_percentage = string2Double( commands[2] );
			if( cutoff_percentage <= 0.0 || cutoff_percentage > 1.0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Cutoff percentage of " + commands[0] + " command has to be between 0.0 and 1.0." ); }
			
			double bin_width = string2Double( commands[3] );
			if( bin_width <= 0.0 ){ allAbortWithMessage( MPI_COMM_WORLD , "bin_width for " + commands[0] + " command has to be a positive (double) number." ); }
			
			papreca_config.setHeightPercentage( cutoff_percentage );
			papreca_config.setBinWidth( bin_width );
			
			
		
		}else{
			allAbortWithMessage( MPI_COMM_WORLD , "Invalid height_calculation method: " + commands[1] + "." );
		}
		
	}
	
	void executeSpeciesMaxBondsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
	
		/// Sets a maximum number of bonds for a species in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( commands.size( ) != 3 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid species_maxbonds command. Must be species_maxbonds N M. (N is the species and M the maximum permissible number of bonds for that species." ); }
		
		int species = string2Int( commands[1] );
		if( species < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Found negative species type:" + commands[1] + " in " + commands[0] + " command." ); }
		
		int bonds_max = string2Int( commands[2] );
		if( bonds_max < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Found negative species max bonds:" + commands[2] + " in " + commands[0] + " command." ); }
		
		papreca_config.setSpeciesMaxBonds( species , bonds_max );

	}
	
	void executeSpeciesMaxBondTypesCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets a maximum number of bonds of specific bond type for an atom type in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( commands.size( ) != 4 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid species_maxbondtypes command. Must be species_maxbondtypes N M K. (N is the atom species M is the bond type K is the maximum number of permissible bonds of type M for species N." ); }
		
		int species = string2Int( commands[1] );
		if( species < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Found negative species type:" + commands[1] + " in " + commands[0] + " command." ); }
		
		int bond_type = string2Int( commands[2] );
		if( bond_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Found negative bond type:" + commands[2] + " in " + commands[0] + " command." ); }
		
		int bonds_max = string2Int( commands[3] );
		if( bonds_max < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Found negative species max bonds:" + commands[3] + " in " + commands[0] + " command." ); }
		
		papreca_config.setSpeciesMaxBondTypes( species , bond_type , bonds_max );
	
	}
	
	void executeMinimizePriorCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
	
		/// Sets a prior minimization command in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @note This function will note check if the set minimization command is valid (i.e., if it is an acceptable LAMMPS command). If the user inputs (in the PAPRECA input file) an invalid minimization command, LAMMPS will fail with an error when the illegal minimize command gets executed.
		
		std::string error_message = "Invalid minimize command. Must be minimize_prior no or minimize_prior yes VALID_MINIMIZATION_LAMMPS_COMMAND. PAPRECA will not check the validity of the lammps command. However, you will get an error (and a relevant error message) during runtime if the command is invalid. See here for more info:https://docs.lammps.org/minimize.html.";
		
		if( commands.size( ) < 2 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		if( commands[1] == "yes" ){
			//concatenate all strings from command[2] onwards to a new string (which is going to be our minimize prior lammps command
			std::string lammps_command = getConcatenatedStringWithSpacesFromStringsVector( commands , 2 , commands.size( ) );
			papreca_config.setMinimize1( lammps_command );
			
		}else if( commands[1] != "no" ){
			allAbortWithMessage( MPI_COMM_WORLD , error_message );
		}
		
	}
	
	void executeMinimizeAfterCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
	
		/// Sets the after minimization command in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @note This function will note check if the set minimization command is valid (i.e., if it is an acceptable LAMMPS command). If the user inputs (in the PAPRECA input file) an invalid minimization command, LAMMPS will fail with an error when the illegal minimize command gets executed.
		
		std::string error_message = "Invalid minimize command. Must be minimize_after no or minimize_after yes VALID_MINIMIZATION_LAMMPS_COMMAND. PAPRECA will not check the validity of the lammps command. However, you will get an error (and a relevant error message) during runtime if the command is invalid. See here for more info:https://docs.lammps.org/minimize.html.";
		if( commands.size( ) < 2 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		if( commands[1] == "yes" ){
			//concatenate all strings from command[2] onwards to a new string (which is going to be our minimize prior lammps command)
			std::string lammps_command = getConcatenatedStringWithSpacesFromStringsVector( commands , 2 , commands.size( ) );
			papreca_config.setMinimize2( lammps_command );
		}else if( commands[1] != "no" ){
			allAbortWithMessage( MPI_COMM_WORLD , error_message );
		
		}
		
	}
	
	void executeTrajectoryDurationCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
	
		/// Sets the LAMMPS trajectory duration in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid trajectory_duration command. Must be trajectory_duration N (where N is an integer denoting the trajectory duration)."); }
		
		int traj_duration = string2Int( commands[1] );
		if( traj_duration <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "The trajectory duration in " + commands[0] + " command has to be a positive integer number." ); }
		
		papreca_config.setTrajDuration( traj_duration );
	
	}
	
	void executeDepoheightsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets the c_height_scan and c_height_reject parameters in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		
		if( commands.size( ) != 3 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid depoheights command. Must be depoheights height_scan height_reject. Scan for deposition events between film_height - height_scan and film_height + height scan. Reject deposition candidates above film_height + height_reject"); }
		
		if( papreca_config.getHeightMethod( ).empty( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Cannot setup depoheights before defining a height calculation method." ); }
		
		const double c_height_scan = string2Double( commands[1] );
		if( c_height_scan <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "height_scan in " + commands[0] + " command has to be a positive (double) number." ); }
		
		const double c_height_reject = string2Double( commands[2] );
		if( c_height_scan <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "height_reject in " + commands[0] + " command has to be a positive (double) number." ); }

		papreca_config.setDepoHeights( c_height_scan , c_height_reject );
		
	}
	
	void executeRandomDepovecsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets the type of deposition vector (random or above) in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid random_depovecs command. Must be random_depovecs yes/no."); }
		
		const bool random_depovecs = string2Bool( commands[1] );
		papreca_config.setRandomDepoVecs( random_depovecs );
		
	}
	
	void executeRandomDiffvecsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets the type of diffusion vector in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( commands.size( ) != 2 && commands.size( ) != 3 ){ allAbortWithMessage( MPI_COMM_WORLD , "Invalid random_diffvecs command. Must be random_diffvecs yes/no. Optional keyword(s): diffvecs_style (2D/3D). Choose 2D for random diffvecs ONLY above the parent atom or 3D for random diffvecs anywhere in the 3D space."); }
	
		const bool random_diffvecs = string2Bool( commands[1] );
		papreca_config.setRandomDiffVecs( random_diffvecs );
		
		if( commands.size( ) == 3 ){ 
			
			if( commands[2] == "2D" || commands[2] == "3D" ){
				papreca_config.setRandomDiffVecsStyle( commands[2] );
			}else{
				allAbortWithMessage( MPI_COMM_WORLD , "Unknown random diffvecs style: " + commands[2] + " in command: " + commands[0] + " the only supported options are 2D and 3D." );
			}
		}
		
	}
	
	
	void executeCreateBondBreakCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
	
		/// Initializes a PAPRECA::PredefinedBondBreak event template and inserts it in the PAPRECA::PredefinedEventsCatalog of the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid create_BondBreak command. Must be create_BondBreak atom1_type atom2_type bond_type rate_(valid rate calc option). Optional argument(s): 1) catalyzed Ntypes types(1-Ntypes) (separate types by spaces:e.g., catalyzed 3 7 8 10), limit length_equil length_perc";
		
		
		if( commands.size( ) < 6 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		int atom1_type = string2Int( commands[1] );
		int atom2_type = string2Int( commands[2] );
		if( atom1_type < 0 || atom2_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Detected non-positive atom type in " + commands[0] + " command." ); }
		
		int bond_type = string2Int( commands[3] );
		if( bond_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Detected non-positive bond type in " + commands[0] + " command." ); }
		
		int current_pos = 4;
		double rate = getRateFromInputRateOptions( commands , current_pos );
		
		std::vector< int > catalyzing_types;
		double length_equil = 0.0;
		double length_perc = 0.0;
		
		//Optional Commands update the current_pos value. Exit when current_pos reached the end of the command line (or if an error occurs).
		if( current_pos != commands.size( ) ){ //Exit immediately if there are no optional keywords
		
			std::unordered_set< std::string > processed; //Initialize set to keep track of processed optional commands
			
			do{
				if( !elementIsInUnorderedSet( processed , std::string( "catalyzed" ) ) && commands[current_pos] == "catalyzed" ){	
					processCatalyzedOption( commands , current_pos , catalyzing_types ); //If no catalyzing types are detected the size of catalyzing_types is 0. initReaction (from PaprecaConfig) takes this into account when initializing the event
					processed.insert( std::string( "catalyzed" ) );
				}else if( !elementIsInUnorderedSet( processed , std::string( "limit" ) ) && commands[current_pos] == "limit" ){
					processBondLimitOption( commands , current_pos , length_equil , length_perc );
					processed.insert( std::string( "limit" ) );
				}else if( commands[current_pos] != "catalyzed" && commands[current_pos] != "limit" ){ allAbortWithMessage( MPI_COMM_WORLD , "Unknown option " + commands[current_pos] + " for command " + commands[0] + "." ); }
				
			}while( commands.size( ) > current_pos );
		}
		
		
		papreca_config.initPredefinedReaction( atom1_type , atom2_type , bond_type , rate , catalyzing_types , length_equil , length_perc );

	
	}
	
	void executeCreateBondFormCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Initializes a PAPRECA::PredefinedBondForm event template and inserts it in the PAPRECA::PredefinedEventsCatalog of the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @note In the current version the catalyzed option for bond formation events is not available. It can be passed from the input file to papreca_config, but this is left here as a placeholder. In the documentation we do not even mention that the catalyzed option (keyword) can be used with bond-formation events.
		
		std::string error_message = "Invalid create_BondForm command. Must be create_BondForm atom1_type atom2_type bond_type bond_dist delete_atoms(yes/no) lone_candidates(yes/no) same_mol(yes/no) rate_(valid rate calc option). Optional argument(s): 1) catalyzed Ntypes types(1-Ntypes) (separate types by spaces:e.g., catalyzed 3 7 8 10).";
		
		if( commands.size( ) < 10 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		int atom1_type = string2Int( commands[1] );
		int atom2_type = string2Int( commands[2] );
		if( atom1_type < 0 || atom2_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Detected non-positive atom type in " + commands[0] + " command." ); }
				
		int bond_type = string2Int( commands[3] );
		if( bond_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Detected non-positive bond type in " + commands[0] + " command." ); }
		
		double bond_dist = string2Double( commands[4] );
		if( bond_dist <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "Bond distance in " + commands[0] + " command has to be a positive (double) number." ); }
		
		int delete_atoms = boolString2Int( commands[5] );
		int lone_candidates = boolString2Int( commands[6] );
		const bool same_mol = string2Bool( commands[7] );
		int current_pos = 8;
		double rate = getRateFromInputRateOptions( commands , current_pos );
		
		std::vector< int > catalyzing_types;
		double length_equil = 0.0;
		double length_perc = 0.0;
		
		//Optional Commands update the current_pos value. Exit when current_pos reached the end of the command line (or if an error occurs).
		if( current_pos != commands.size( ) ){ //Exit immediately if there are no optional keywords
		
			std::unordered_set< std::string > processed; //Initialize set to keep track of processed optional commands
			
			do{
				if( !elementIsInUnorderedSet( processed , std::string( "catalyzed" ) ) && commands[current_pos] == "catalyzed" ){	
					processCatalyzedOption( commands , current_pos , catalyzing_types ); //If no catalyzing types are detected the size of catalyzing_types is 0. initReaction (from PaprecaConfig) takes this into account when initializing the event
					processed.insert( std::string( "catalyzed" ) );
				}else if( !elementIsInUnorderedSet( processed , std::string( "limit" ) ) && commands[current_pos] == "limit" ){
					processBondLimitOption( commands , current_pos , length_equil , length_perc );
					processed.insert( std::string( "limit" ) );
				}else if( commands[current_pos] != "catalyzed" && commands[current_pos] != "limit" ){ allAbortWithMessage( MPI_COMM_WORLD , "Unknown option " + commands[current_pos] + " for command " + commands[0] + "." ); }
				
			}while( commands.size( ) > current_pos );
		}
	
		papreca_config.initPredefinedBondForm( atom1_type , atom2_type , bond_type , bond_dist , delete_atoms , lone_candidates , same_mol , rate , catalyzing_types , length_equil , length_perc );

		
		
	}
	
	void executeCreateDiffusionHopCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Initializes a PAPRECA::PredefinedDiffusionHop event template and inserts it in the PAPRECA::PredefinedEventsCatalog of the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid create_DiffusionHop command. Must be create_DiffusionHop parent_type velocity diffusion_distance is_displacive(yes/no) diffused_type rate_(valid rate calculation command).";
		
		if( commands.size( ) < 8 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		int parent_type = string2Int( commands[1] );
		if( parent_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "parent_type in " + commands[0] + " command has to be a non-negative integer number."); }
		
		double velocity = string2Double( commands[2] );
		if( velocity < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "velocity in " + commands[0] + " command has to be non-negative."); }

		double distance = string2Double( commands[3] );
		bool is_displacive = string2Bool( commands[4] );
		
		int diffused_type = string2Int( commands[5] );
		if( diffused_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "diffused_type in " + commands[0] + " command has to be a non-negative integer number."); }
		
		int current_pos = 6;
		double rate = getRateFromInputRateOptions( commands , current_pos );
		std::string custom_style = "NONE";
		std::vector< int > style_atomtypes;
		
		//Optional Commands update the current_pos value. Exit when current_pos reached the end of the command line (or if an error occurs).
		if( current_pos != commands.size( ) ){ //Exit immediately if there are no optional keywords
			do{
				if( commands[current_pos] == "custom" ){	
					processCustomDiffEventOptions( commands , current_pos , custom_style , style_atomtypes );
					
				}else{
					allAbortWithMessage( MPI_COMM_WORLD , "Unknown option " + commands[current_pos] + " for command " + commands[0] + "." );
				}
			}while( current_pos < commands.size( ) );
		}
		
		papreca_config.initPredefinedDiffusionHop( parent_type , velocity , distance , is_displacive , diffused_type , rate , custom_style , style_atomtypes );
		
		
		
	}
	
	void executeCreateDepositionCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Initializes a PAPRECA::PredefinedDeposition event template and inserts it in the PAPRECA::PredefinedEventsCatalog of the PAPRECA::PaprecaConfig object.
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @note This function will not check if a molecule with the provided adsorbate_name was previously initialized in the LAMMPS input file. If the adsorbate nate was note correctly set by the user LAMMPS will fail with an error the moment it attempts to execute the deposition event.
		
		std::string error_message = "Invalid create_Deposition command. Must be create_Deposition parent_type deposition_offset insertion_velocity adsorbate_name rate_(valid rate calculation command). Optional keyword(s): sticking_coeff variable/constant C (where C is the sticking coefficient ONLY to be used if the sticking_coeff is variable).";
		
		if( commands.size( ) < 7 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		int parent_type = string2Int( commands[1] );
		if( parent_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "parent_type in " + commands[0] + " command has to be a non-negative integer number."); }
		
		double depo_offset = string2Double( commands[2] );
		if( depo_offset < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "depo_offset in " + commands[0] + " command has to be non-negative."); }
		
		double insertion_velocity = string2Double( commands[3] );
		if( insertion_velocity < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "insertion_velocity in " + commands[0] + " command has to be non-negative."); }
		
		std::string adsorbate_name = commands[4];
		int current_pos = 5;
		double rate = getRateFromInputRateOptions( commands , current_pos );
		double sticking_coeff = getStickingCoeffFromDepositionEventOptions( commands , current_pos ); //if sticking coeff returns -1.0 this means that the sticking coeff is variable
		bool variable_sticking;
		if( sticking_coeff == -1.0 ){
			variable_sticking = true;
		}else{
			variable_sticking = false;
		}
		
		if( commands.size( ) > current_pos ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		
		//VERY IMPORTANT to compute center
		computeMolCenter( lmp , adsorbate_name.c_str( ) ); //invoke a mol center calculation to update the mol center from (0,0,0) to the actual mol center, before the run
		

		papreca_config.initPredefinedDeposition( lmp , parent_type , depo_offset , insertion_velocity , adsorbate_name , rate , variable_sticking , sticking_coeff );
			
		
		
	}
	
	void executeCreateMonoatomicDesorptionCommand( std::vector< std:: string > &commands , PaprecaConfig &papreca_config ){
		
		/// Initializes a PAPRECA::PredefinedMonoatomicDesorption event template and inserts it in the PAPRECA::PredefinedEventsCatalog of the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid create_MonoatomicDesorption command. Must be create_MonoatomicDesorption parent_type rate_(valid rate calculation command).";
		
		if( commands.size( ) < 4 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		int parent_type = string2Int( commands[1] );
		if( parent_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "parent_atom in " + commands[0] + " command has to be a non-negative integer number." ); }
		
		int current_pos = 2;
		double rate = getRateFromInputRateOptions( commands , current_pos );
		
		if( commands.size( ) > current_pos ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		papreca_config.initPredefinedMonoatomicDesorption( parent_type , rate );
		
	}
	
	void executeExportHeightVtimeCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Enables the generation of PAPRECA::HeightVtime files with a give frequency in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid export_HeightVtime command. Must be export_HeightVtime N (where N is the export frequency:i.e., every N steps we write to the file).";
		
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }

		int print_freq = string2Int( commands[1] );
		if( print_freq <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , " print_freq in " + commands[0] + " command has to be a positive integer number." ); }
		
		HeightVtime &heightVtime_file = papreca_config.getHeightVtimeFile( );
		
		heightVtime_file.setActive( );
		heightVtime_file.setPrintFreq( print_freq );
		
		
	}
	
	void executeExportSurfaceCoverageCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Enables the generation of PAPRECA::SurfaceCoverage files with a give frequency in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid export_SurfCoverage command. Must be export_SurfCoverage N (where N is the export frequency:i.e., every N steps we write to the file).";
		
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }

		int print_freq = string2Int( commands[1] );
		if( print_freq <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , " print_freq in " + commands[0] + " command has to be a positive integer number." ); }
		
		SurfaceCoverage &surfcoverage_file = papreca_config.getSurfaceCoverageFile( );
		
		surfcoverage_file.setActive( );
		surfcoverage_file.setPrintFreq( print_freq );
		
	}
	
	void executeExportElementalDistributionsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Enables the generation of PAPRECA::ElementalDistribution files with a give frequency in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid export_ElementalDistributions command. Must be export_ElementalDistributions N (where N is the export frequency: i.e., every N steps we write to the file). Optional keyword(s): bin_width M (where M is a double number).";
		if( commands.size( ) < 2 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		int print_freq = string2Int( commands[1] );
		if( print_freq <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , " print_freq in " + commands[0] + " command has to be a positive integer number." ); }
		
		ElementalDistribution &elementalDistributions_file = papreca_config.getElementalDistributionsFile( );
		
		elementalDistributions_file.setActive( );
		elementalDistributions_file.setPrintFreq( print_freq );
		
		int current_pos = 2;
		
		
		//Optional Commands update the current_pos value. Exit when current_pos reached the end of the command line (or if an error occurs).
		if( commands.size( ) != current_pos ){
			do{
				if( commands[current_pos] == "bin_width" ){	
					processBinWidthOptionForElementalDistributions( commands , papreca_config , current_pos );
					
				}else{
					allAbortWithMessage( MPI_COMM_WORLD , "Unknown option " + commands[current_pos] + " for command " + commands[0] + "." );
				}
			}while( current_pos < commands.size( ) );
		}
		
		
	}
	
	void executeExportExecutionTimesCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Enables the generation of PAPRECA::ExecTime files with a give frequency in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid export_ExecTimes command. Must be export_ExecTimes N (where N is the export frequency:i.e., every N steps we write to the file).";
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		int print_freq = string2Int( commands[1] );
		if( print_freq <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , " print_freq in " + commands[0] + " command has to be a positive integer number." ); }
		
		ExecTime &execTime_file = papreca_config.getExecTimeFile( );
		
		execTime_file.setActive( );
		execTime_file.setPrintFreq( print_freq );
		
		
		
	}
	
	void executeRestartFreqCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Enables the generation of LAMMPS restart files with a give frequency in the PAPRECA::PaprecaConfig object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid restart_freq command. Must be restart_freq N (where N is the dump restart frequency:i.e., every N steps a restart file is dumped).";
		if( commands.size( ) != 2 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		int restart_freq = string2Int( commands[1] );
		if( restart_freq <= 0 ){ allAbortWithMessage( MPI_COMM_WORLD , " restart_freq in " + commands[0] + " command has to be a positive integer number." ); }
		
		papreca_config.setRestartDumpFreq( restart_freq );
		
		
	}
	
	void executeSigmasOptionsCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets the sigma options in the PAPRECA::PaprecaConfig object.
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		std::string error_message = "Invalid sigmas_options command. Must be sigmas_options LAMMPS/manual. Optional keyword(s):mix geom/arithm (for geometric and arithmetic mixing of sigmas).";
		if( commands.size( ) < 2 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		std::string style = commands[1];
		papreca_config.setSigmaStyle( style );
		
		if( style == "LAMMPS" ){
			papreca_config.initSigmasFromLammps( lmp );
		}else if( style != "manual" ){ 
			allAbortWithMessage( MPI_COMM_WORLD , error_message );
		}
		
		
		int current_pos = 2;
		
		
		if( commands.size( ) != current_pos ){
			//Optional Commands update the current_pos value. Exit when current_pos reached the end of the command line (or if an error occurs).
			do{
				if( commands[current_pos] == "mix" ){	
					processSigmaMixOptions( commands , papreca_config , current_pos );
					
				}else{
					allAbortWithMessage( MPI_COMM_WORLD , "Unknown option " + commands[current_pos] + " for command " + commands[0] + "." );
				}
			}while( current_pos < commands.size( ) );
		}
		
		
	}
	
	void executeInitSigmaCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// Sets the sigma value for a pair of atom types in the PAPRECA::PaprecaConfig object.
		/// @param[in] lmp pointer to previously initialized LAMMPS object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		if( papreca_config.getSigmaStyle( ) == "" ){ allAbortWithMessage( MPI_COMM_WORLD , "Use of init_sigma command without prior use of sigmas_options command. Please set sigmas_option first before initializing sigmas." ); }
		if( papreca_config.getSigmaStyle( ) != "manual" ){ warnAll( MPI_COMM_WORLD , "You are attempting to manually init a sigma in a init_sigma command but the set sigmas option WAS NOT MANUAL. You may be modify existing values..." ); }
		
		std::string error_message = "Invalid init_sigma command. Must be init_sigma N M sigma (N is atom1_type M is atom2_type)";
		if( commands.size( ) != 4 ){ allAbortWithMessage( MPI_COMM_WORLD , error_message ); }
		
		const int atom1_type = string2Int( commands[1] );
		if( atom1_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , " atom1_type in " + commands[0] + " command has to be a non-negative integer number." ); }
		
		const int atom2_type = string2Int( commands[2] );
		if( atom2_type < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , " atom1_type in " + commands[0] + " command has to be a non-negative integer number." ); }
		
		const double sigma = string2Double( commands[3] );
		if( sigma < 0 ){ allAbortWithMessage( MPI_COMM_WORLD , " sigma in " + commands[0] + " command has to be a non-negative (double) number as it represents a distance." ); }
		
		papreca_config.setSpeciesPair2Sigma( atom1_type , atom2_type , sigma );
	
		
	}
	
	
	void executePaprecaCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config ){
		
		/// General PAPRECA input function that redirects commands to the relevant command function.
		/// @param[in] lmp pointer to previously initialized LAMMPS object.
		/// @param[in] commands trimmed/processed vector of strings. This is effectively the entire command line with each vector element (i.e., std::string) being a single word/number.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @note The first word (string) in the commands vector of strings ( i.e., commands[0] ) indicates the type of command.
		
		//Use commands[0] to redirect the command to the appropriate function for execution.
		//When invalid commands are encountered the code Aborts with an error.
		
		std::string command_class = commands[0];
		
		
		//Acceptable commands
		if( command_class == "KMC_steps" ){
			executeKMCstepsCommand( commands , papreca_config );
		}else if( command_class == "KMC_per_MD" ){ 
			executeKMCperMDCommand( commands , papreca_config );
		}else if( command_class == "time_end" ){
			executeTimeEndCommand( commands , papreca_config );
		}else if( command_class == "random_seed" ){
			executeRandomSeedCommand( lmp , commands, papreca_config );
		}else if( command_class == "fluid_atomtypes" ){
			executeFluidAtomTypesCommand( commands , papreca_config );
		}else if( command_class == "frozen_atomtypes" ){
			executeFrozenAtomTypesCommand( commands , papreca_config );
		}else if( command_class == "desorption" ){
			executeDesorptionCommand( commands, papreca_config );
		}else if( command_class == "height_calculation" ){
			executeHeightCalculationCommand( commands , papreca_config );
		}else if( command_class == "species_maxbonds" ){
			executeSpeciesMaxBondsCommand( commands, papreca_config );
		}else if( command_class == "species_maxbondtypes" ){
			executeSpeciesMaxBondTypesCommand( commands , papreca_config );
		}else if( command_class == "minimize_prior" ){
			executeMinimizePriorCommand( commands , papreca_config );
		}else if( command_class == "minimize_after" ){
			executeMinimizeAfterCommand( commands , papreca_config );
		}else if( command_class == "trajectory_duration" ){
			executeTrajectoryDurationCommand( commands , papreca_config );
		}else if( command_class == "depoheights" ){
			executeDepoheightsCommand( commands , papreca_config );
		}else if( command_class == "random_depovecs" ){
			executeRandomDepovecsCommand( commands , papreca_config );
		}else if( command_class == "random_diffvecs" ){
			executeRandomDiffvecsCommand( commands, papreca_config );
		}else if( command_class == "create_BondBreak" ){
			executeCreateBondBreakCommand( commands , papreca_config );
		}else if( command_class == "create_BondForm" ){
			executeCreateBondFormCommand( commands, papreca_config );
		}else if( command_class == "create_DiffusionHop" ){
			executeCreateDiffusionHopCommand( commands, papreca_config );
		}else if( command_class == "create_Deposition" ){
			executeCreateDepositionCommand( lmp , commands , papreca_config );
		}else if( command_class == "create_MonoatomicDesorption" ){
			executeCreateMonoatomicDesorptionCommand( commands , papreca_config );
		}else if( command_class == "export_HeightVtime" ){
			executeExportHeightVtimeCommand( commands , papreca_config );
		}else if( command_class == "export_SurfCoverage" ){
			executeExportSurfaceCoverageCommand( commands , papreca_config );
		}else if( command_class == "export_ElementalDistributions" ){
			executeExportElementalDistributionsCommand( commands , papreca_config );
		}else if( command_class == "export_ExecutionTimes" ){
			executeExportExecutionTimesCommand( commands , papreca_config );
		}else if( command_class == "restart_freq" ){
			executeRestartFreqCommand( commands , papreca_config );
		}else if( command_class == "sigmas_options" ){
			executeSigmasOptionsCommand( lmp , commands , papreca_config );
		}else if( command_class == "init_sigma" ){
			executeInitSigmaCommand( lmp , commands , papreca_config );
		}else{ 
			allAbortWithMessage( MPI_COMM_WORLD , "Invalid PAPRECA command:" + command_class + " in PAPRECA input file." );
		}
		
	}
	
	
	std::vector< std::string > processLine( char *line ){
		
		/// Receives a line from the PAPRECA input file, removes blank characters and comments (i.e., text to the right of a '#' character). Then, converts the char* to an std::vector< std::string > container that stores each word/command/keyword on a separate std::string.
		/// @param[in] line raw (unprocessed command-line sent from the PAPRECA input file.
		/// @returns vector of string container with processed line (command).
		
		//First we convert the traditional c-style string to an std string for easier processing (We used a traditional string before for easier communication of chars through MPI).
		std::string line_str( line );
		
		//Before executing we need to process the line
		//Firstly, anything to the right of '#' characters is considered a commend and thus trimmed from the string.
		size_t comment_pos = line_str.find( '#' ); //Ignore everything to the right of # (i.e., comments)
		if( comment_pos != std::string::npos ){ line_str = line_str.substr( 0 , comment_pos ); }
				
		//Use stringstream to split the line into words (i.e., specific commands)
		std::stringstream string_stream( line_str );
		std::vector< std::string > commands; //Contains the full command (i.e., many words/strings)
		std::string command; //Contains the individual words/commands
		while( string_stream >> command ){ commands.push_back( command ); }
		
		//At this point the commands vector of strings contains the full command (i.e., the whole line).
		//Hashtags and blank spaces have been trimmed so we can proceed with the execution of the command
		return commands;
		
	}
	
	
	void warn4IllegalRuns( const int &proc_id , PaprecaConfig &papreca_config ){
	
		/// After processing all commands, this function is called to warn all MPI processes for potentially illegal PAPRECA runs.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
	
		if( papreca_config.predefinedCatalogIsEmpty( ) ){ warnAll( MPI_COMM_WORLD , "No predefined events were defined!" ); }
		if( papreca_config.getKMCperMD( ) !=0 && papreca_config.getTrajDuration( ) == 0 && papreca_config.getMinimize1( ).empty( ) && papreca_config.getMinimize2( ).empty( ) ){ warnAll( MPI_COMM_WORLD , "KMC per MD defined but not equilibration scheme set (i.e., trajectory duration is 0, and no prior or after minimization commands were set" ); }
	
	}
	
	void abortIllegalRun( const int &proc_id , PaprecaConfig &papreca_config ){
	
		/// After processing all commands, this function is called abort illegal PAPRECA runs.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		
		//Basic settings aborts
		if( papreca_config.getKMCsteps( ) == 0 ){ allAbortWithMessage( MPI_COMM_WORLD , "PAPRECA KMC steps were not set (or set to zero)." ); }
		
		//Random number generator aborts
		if( !papreca_config.ranNumGeneratorIsInitialized( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Random number generator was not initialized properly. This typically indicates that the random seed WAS NOT provided in the PAPRECA input file." ); }
		
		//Sigmas aborts
		if( ( papreca_config.predefinedCatalogHasBondBreakEvents( ) || papreca_config.predefinedCatalogHasBondFormEvents( ) ) && papreca_config.type2SigmaMapIsEmpty( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Bond form/break events were defined but no sigma_options and sigmas were detected! Please defined sigmas." ); }
		
		//Film height aborts
		if( papreca_config.getHeightMethod( ).empty( ) && papreca_config.predefinedCatalogHasDepositionEvents( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Cannot setup deposition events without setting up a film calculation method." ); }
		if( papreca_config.getHeightMethod( ).empty( ) && !papreca_config.getDesorptionStyle( ).empty( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Cannot setup desorption without setting up a film calculation method." ); }
		if( papreca_config.getHeightMethod( ).empty( ) && papreca_config.getHeightVtimeFile( ).isActive( ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Cannot dump a heightVtime file without setting up a film calculation method." ); }
		
	
	}
	
	void readInputAndInitPaprecaConfig( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const char *file_name , PaprecaConfig &papreca_config ){
		
		/// Reads PAPRECA input line-by-line on the master MPI process and casts each line to all other processes. Then, all MPI processes execute the relevant PAPRECA command in executePaprecaCommand(). In the end, all MPI processes initialize an identical copy of PAPRECA::PaprecaConfig (containing the settings and global variables for the PAPRECA run).
		/// @param[in] lmp pointer to previously instantiated LAMMPS object.
		/// @param[in] proc_id ID of current MPI process.
		/// @param[in] file_name name of PAPRECA input file.
		/// @param[in,out] papreca_config previously instantiated PAPRECA::PaprecaConfig object storing the settings and global variables for the PAPRECA simulation.
		/// @note After initializing the PAPRECA::PaprecaConfig object in all procs, this function warns for or aborts illegal runs.
		
		//Open input file
		FILE *input_file = NULL;
		if( proc_id == 0 ){
			
			input_file = fopen( file_name , "r" );
			if( input_file == NULL ){ allAbortWithMessage( MPI_COMM_WORLD , " Could not open PAPRECA input sciprt." ); }
		}
		
		//Master proc reads the commands, line by line and then Bcasts information to other procs. Each proc calls the command functions separately so PaprecaConfig is correctly initialized (with the same information) on all procs.
		int line_length;
		char line[1024];
		while( true ){
			
			//Read the line and get line length
			if( proc_id == 0 ){
				
				if( std::fgets( line , 1024 , input_file ) == NULL ){
						line_length = 0;
				}else{
					
					line_length = strlen( line ) + 1;
				}
				if( line_length == 0 ){ fclose( input_file ); }
				
			}
			//Bcast line to all other procs
			MPI_Bcast( &line_length , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
			if( line_length == 0 ){ break; } //Exit immediately if line_length is 0. This means that we have reached the end of file since fgets returns NULL.
			MPI_Bcast( line , line_length , MPI_CHAR , 0 , MPI_COMM_WORLD );
			//Now every proc know what's in the line. You can go ahead and process then line, and then execute PAPRECA commands on all procs to initialize PaprecaConfig
			std::vector< std::string > commands = processLine( line );
			if( !commands.empty( ) ){ executePaprecaCommand( lmp , commands , papreca_config ); } //Execute command (unless the line was empty due to commends or blank lines).
			
				
		}
			
		
		papreca_config.mixSigmas( lmp ); //Mix sigmas before the start of the simulation and AFTER reading all the commands
		setTimeUnitsConversionConstant( lmp , papreca_config );
		
		//After reading the whole input file, prevent simulations that lead to runtime errors (e.g., dereference of NULL ptr, etc.).
		//Also, warn if certain settings MIGHT cause runtime errors.
		
		warn4IllegalRuns( proc_id , papreca_config );
		abortIllegalRun( proc_id , papreca_config );
		
			
			
	}	
		
}//End of namespace PAPRECA
