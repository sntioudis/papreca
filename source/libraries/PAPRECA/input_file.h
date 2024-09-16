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
/// @brief Declarations of functions related to the reading of the PAPRECA input file, and the initialization of the PAPRECA::PaprecaConfig object from the PAPRECA input file.

#ifndef INPUT_FILE_H
#define INPUT_FILE_H

//System Headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <mpi.h>

//LAMMPS Headers
#include "lammps.h"
#include "random_mars.h"
#include "update.h"

//kMC Headers
#include "papreca_error.h"
#include "papreca_config.h"
#include "event_list.h"
#include "rates_calc.h"
#include "utilities.h"
#include "lammps_wrappers.h"

namespace PAPRECA{
	
	//SUPPLEMENTARY SETTER FUNCTIONS FOR APRECAT CONFIG
	void setTimeUnitsConversionConstant( LAMMPS_NS::LAMMPS *lmp , PaprecaConfig &papreca_config );
	
	//ACCEPTABLE OPTIONAL COMMAND KEYWORDS and KEYWORD IDENTIFICATION FUNCTIONS.
	void checkForAcceptableKeywordsUsedMultipleTimes( std::vector< std::string > &commands , const std::string &keyword );
	void check4AcceptableKeywords( std::vector< std::string > &commands , const int &start , std::unordered_set< std::string > &acceptable_keywords , const bool &accept_bool );
	void processCatalyzedOption( std::vector< std::string > &commands , int &current_pos , std::vector< int > &catalyzing_types );
	void processSigmaMixOptions( std::vector< std::string > &commands , int &current_pos );
	double getBinWidthFromElementalDistributions( std::vector< std::string > &commands , int &current_pos );
	double getStickingCoeffFromDepositionEventOptions( std::vector< std::string > &commands , int &current_pos );
	double getRateFromInputRateOptions( std::vector< std::string > &commands , int &current_pos );
	void processCustomDiffEventOptions( std::vector< std::string > &commands , int &current_pos , std::string &custom_style , std::vector< int > &style_atomtypes );
	
	//ACCEPTABLE COMMANDS
	void executeKMCstepsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeKMCperMDCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeTimeEndCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeRandomSeedCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeFluidAtomTypesCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeFrozenAtomTypesCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeDesorptionCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeHeightCalculationCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeSpeciesMaxBondsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeSpeciesMaxBondTypesCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeMinimizePriorCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeMinimizeAfterCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeTrajectoryDurationCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeDepoheightsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeRandomDepovecsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeRandomDiffvecsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeCreateBondBreakCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeCreateBondFormCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeCreateDiffusionHopCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeCreateDepositionCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeCreateMonoatomicDesorptionCommand( std::vector< std:: string > &commands , PaprecaConfig &papreca_config );
	void executeExportHeightVtimeCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeExportSurfaceCoverageCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeExportElementalDistributionsCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeExportExecutionTimesCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeRestartFreqCommand( std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeSigmasOptionsCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	void executeInitSigmaCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	
	//General command execution functions
	void executePaprecaCommand( LAMMPS_NS::LAMMPS *lmp , std::vector< std::string > &commands , PaprecaConfig &papreca_config );
	std::vector< std::string > processLine( char *line );
	void abortIllegalRun( const int &proc_id , PaprecaConfig &papreca_config );
	void readInputAndInitPaprecaConfig( LAMMPS_NS::LAMMPS *lmp , const int &proc_id , const char *file_name , PaprecaConfig &papreca_config ); //Main function reading the input file and passing information (i.e., commands) to the execute commands function

	
}//end of namespace PAPRECA

#endif

