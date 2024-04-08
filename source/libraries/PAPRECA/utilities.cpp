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
///@brief Definitions for utilities.h

#include "utilities.h"

namespace PAPRECA{
	
	//Array/Vector functions
	double doubleArrSum( double *arr , const int &size ){
		
		/// Receives a C-style array of doubles and returns the sum of its elements.
		/// @param[in] arr pointer to 1-D array of doubles.
		/// @param[in] size size (length) of arr.
		/// @return sum of all elements of arr.
		
		double sum = 0;
		
		for( int i = 0; i < size; ++i ){
			
			sum += arr[i];
			
		}

		return sum;	
	}
	
	
	void copyDoubleArray3D( double *copy , const double *source , const int start , const int end ){
		
		/// Receives an array of doubles and copies it to a different array of doubles.
		/// @param[in,out] copy copied array.
		/// @param[in] source original array.
		/// @param[in] start copy starting from that index.
		/// @param[in] end copy all elements whose index is less than this variable (i.e., copy[i] with i < end ).
		/// @note This function will result in a segmentation fault if the "start" or "end" indexes are out-of-bounds (i.e., not mapped for the source array).
		
		for( int i = start; i < end; ++i ){
			
			copy[i] = source[i];
			
		}
		
	}
	
	//String management
	const std::string getConcatenatedStringFromStringsVector( const std::vector< std::string > &strings , size_t start = -1 , size_t end = -1 ){
		
		/// Receives an std::vector of strings and returns a concatenated std::string (i.e., a string combining all elements of the std::vector< std::string > container)
		/// @param[in] strings std::vector< std::string > container
		/// @param[in] start concatenate strings starting for this index in the strings vector.
		/// @param[in] end concatenate strings whose index in the strings vector is less than end.
		/// @return std::string of concatenated strings from input std::vector< std::string > container.
		/// @see getConcatenatedStringWithSpacesFromStringsVector()
		
		//This means that we've entered the function using the default values. In such case we concatenate the whole vector of strings
		if( start == -1 && end == -1 ){
			start = 0;
			end = strings.size( );
		
		}else if( start == -1 ){
			start = 0;
		
		}else if( end == -1 ){
		
			end = strings.size( );
		}
		
		// Ensure start and end indices are within bounds.
		if ( start >= strings.size() || end > strings.size() || start > end || start < 0 ) {
			allAbortWithMessage( MPI_COMM_WORLD , "Attempted to concatenate string with incorrect start/end bounds.");
		}

		std::string result;

		for (size_t i = start; i < end; ++i) {
			result += strings[i];
		}

		return result;
	
	}
	
	const std::string getConcatenatedStringWithSpacesFromStringsVector( const std::vector< std::string > &strings , size_t start = -1 , size_t end = -1 ){
		
		/// Same as getConcatenatedStringFromStringsVector(). However, the present function also adds a space character between each concatenated string (received from std::vector< string > ).
		/// @param[in] strings std::vector< std::string > container
		/// @param[in] start concatenate strings starting for this index in the strings vector.
		/// @param[in] end concatenate strings whose index in the strings vector is less than end.
		/// @return std::string of concatenated strings from input std::vector< std::string > container.
		/// @see getConcatenatedStringFromStringsVector()
	
		//This means that we've entered the function using the default values. In such case we concatenate the whole vector of strings
		if( start == -1 && end == -1 ){
			start = 0;
			end = strings.size( );
		
		}else if( start == -1 ){
			start = 0;
		
		}else if( end == -1 ){
		
			end = strings.size( );
		}
		
		// Ensure start and end indices are within bounds.
		if ( start >= strings.size() || end > strings.size() || start > end || start < 0 ) {
			allAbortWithMessage( MPI_COMM_WORLD , "Attempted to concatenate string with incorrect start/end bounds.");
		}

		std::string result;

		for (size_t i = start; i < end; ++i) {
			result += strings[i];
			if( i != end-1 ){
				result += " "; //Add a space to all strings except for the last one
			}
		}

		return result;
	
	}
	
	const int getStringPosInStringVec( const std::string &string2check , const std::vector< std::string > &strings ){
		
		///Returns the position of a string in a vector of strings (std::vector< std::string > container).
		/// @param[in] string2check string to identify.
		/// @param[in] strings std::vector< std::string > container.
		/// @return position (index) of string in the std::vector< std::string > container or -1 if the string2check is not in the vector of strings.
		/// @note This function is useful to check for keywords in command lines (e.g., PAPRECA input file).
		
		int string_pos = -1;
		
		for( int i = 0; i < strings.size( ); ++i ){
			
			if( strings[i] == string2check ){ return i; }
			
		}
		
		return string_pos;
		
	}
	
	bool stringIsNumber( const std::string &string ){
		
		/// This function checks if the string is a number. We allow negative/positive signs and decimal points in there so we can potentially indetify double numbers or negative numbers or exponentials as well.
		/// @param[in] string input std::string.
		/// @return true if the string is a number or false if it is not.
		/// @note See this post in stackoverflow: https://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c#:~:text=The%20most%20efficient%20way%20would,the%20string%20not%20a%20number.&text=As%20pointed%20out%20in%20the,only%20works%20for%20positive%20integers.
		/// @see stringIsIntNumber()
		
		std::string::const_iterator it = string.begin( );
		while( it != string.end( ) && ( std::isdigit( *it ) || *it == '.' || *it == '+' || *it == '-' || *it == 'e' || *it == 'E' ) ) ++it;
		return !string.empty( ) && it == string.end( );
		
	}
	
	bool stringIsIntNumber( const std::string &string ){
		
		
		/// This function checks if the string is a number.
		/// @param[in] string std::string to be checked.
		/// @return false if a decimal or exponential is detected. true otherwise.
		/// @note Same as stringIsNumber() but now we allow positive/negative signs.
		/// @see stringIsNumber()
		
		std::string::const_iterator it = string.begin( );
		while( it != string.end( ) && ( std::isdigit( *it ) || *it == '+' || *it == '-' ) ) ++it;
		return !string.empty( ) && it == string.end( );
		
	}
	
	const bool stringIsBool( const std::string &string ){
		
		/// Checks if the string is a bool (i.e., yes or no).
		/// @param[in] string std::string to be checked.
		/// @return true if the the string is "yes" or "no". False otherwise.
		
		return ( string == "yes" || string == "no" ) ? true : false;
		
	}
	
	const unsigned long int string2UnsignedLongInt( std::string &string ){
		
		/// Converts input string to unsigned long int.
		/// @param[in] string std::string to be converted.
		/// @return converted unsigned long integer.
		
		//Check if string contains ONLY digits.
		if( !stringIsNumber( string ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Tried to convert string " + string + " containing invalid characters to integer. This usually indicates that a non-numeric character is mixed-up with numbers in the PAPRECA input file." );  }
		
		unsigned long int num;
		
		try{ 
			num = static_cast< unsigned long int >( std::stoul( string ) );
		}catch( const std::invalid_argument& e ){ 
			allAbortWithMessage( MPI_COMM_WORLD , "Invalid conversion of string: " + string + " to unsigned long int. This error probably indicates an invalid unsigned long int input in the PAPRECA input file." );
		}catch( const std::out_of_range& e){ 
			allAbortWithMessage( MPI_COMM_WORLD , "Conversion of string: " + string + " to unsigned long int led to out-of-range unsigned long integer. This error probably indicates an invalid unsigned long int input i the PAPRECA input file." );
		}
		
		return num;
		
	}
	
	const int string2Int( std::string &string ){
		
		///  Converts an std::string to an integer.
		/// @param[in] string std::string to be converted.
		/// @return converted integer.
		
		//Check if string contains ONLY digits.
		if( !stringIsIntNumber( string ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Tried to convert string " + string + " containing invalid characters to integer. This usually indicates that a non-numeric character OR A DECIMAL POINT is mixed-up with integer numbers in the PAPRECA input file." );  }
		
		int num;
		
		try{ 
			num = static_cast< int >( std::stoi( string ) );
		}catch( const std::invalid_argument& e ){ 
			allAbortWithMessage( MPI_COMM_WORLD , "Invalid conversion of string: " + string + " to int. This error probably indicates an invalid int input in the PAPRECA input file." );
		}catch( const std::out_of_range& e){ 
			allAbortWithMessage( MPI_COMM_WORLD , "Conversion of string: " + string + " to int led to out-of-range integer. This error probably indicates an invalid int input i the PAPRECA input file." );
		}
		
		return num;
			
	}
	
	const double string2Double( std::string &string ){
	
		/// Converts an input std::string to a double number.
		/// @param[in] string std::string to be converted.
		/// @return converted double number.
		
		//Check if string contains ONLY digits.
		if( !stringIsNumber( string ) ){ allAbortWithMessage( MPI_COMM_WORLD , "Tried to convert string " + string + " containing invalid characters to double. This usually indicates that a non-numeric character is mixed-up with numbers in the PAPRECA input file." );  }
			
		double num;
		
		try{ 
			num = static_cast< double >( std::stod( string ) );
		}catch( const std::invalid_argument& e ){ 
			allAbortWithMessage( MPI_COMM_WORLD , "Invalid conversion of string: " + string + " to double. This error probably indicates an invalid double input in the PAPRECA input file." );
		}catch( const std::out_of_range& e){ 
			allAbortWithMessage( MPI_COMM_WORLD , "Conversion of string: " + string + " to double led to out-of-range double. This error probably indicates an invalid double input i the PAPRECA input file." );
		}
		
		return num;
		
	}
	
	const bool string2Bool( std::string &string ){
	
		/// Converts an input std::string to a boolean.
		/// @param[in] string std::string to be converted.
		/// @return converted boolean.
		/// @note This function only works if the string is either yes (i.e., 1 ) or no (i.e., 0 );
		
		if( !stringIsBool( string ) ){
			allAbortWithMessage( MPI_COMM_WORLD , "Tried convert a string that is neither yes nor no to a bool. String: " + string  + "." );
		}else if( string == "yes" ){
			return true;
		}else if( string == "no" ){
			return false;
		}
		
		
		return false; //You should never reach this point. The function either aborts or returns a valid value
		
	}
	
	const int boolString2Int( std::string &string ){
		
		/// This function receives an std::string which is yes/no and converts it to an int (1/0, respectively).
		/// @param[in] string std::string to be converted.
		/// @return converted integer.
		
		if( !stringIsBool( string ) ){
			allAbortWithMessage( MPI_COMM_WORLD , "Tried convert a string that is neither yes nor no to a bool integer (0/1). String: " + string  + "." );
		}else if( string == "yes" ){
			return 1;
		}else if( string == "no" ){
			return 0;
		}
		
		
		return -1; //Indicates an error. You should never reach this point
		
	}
	
}//End of namespace PAPRECA
