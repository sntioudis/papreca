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
///@brief Utility functions (e.g., for arrays, strings etc.) and typedefs used in papreca.cpp main and in the other header files.
 
#ifndef UTILITIES_H
#define UTILITIES_H

//system Headers
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>
#include <mpi.h>

//LAMMPS headers
#include "pointers.h"

//kMC headers
#include "papreca_error.h"

namespace PAPRECA{

	//Functions with template are included ALONG WITH THEIR IMPLEMENTATION in the .h file to avoid undefined reference compiler errors. See: https://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor

	//Array/Vector functions
	double doubleArrSum( double *arr , const int &size );
	void copyDoubleArray3D( double *copy , const double *source , const int start = 0 , const int end = 3 );
	
	template< typename T > 
	T getSumOfVecElements( const std::vector< T > &vec ){
		
		/// Receives an std::vector and returns the sum of its elements.
		/// @param[in] vec std::vector
		/// @return sum of elements.
		/// @note This is a template function so it works with std::vector containers of any primitive type (i.e., integers, doubles, floats, etc.).
		
		T sum = T( ); //Initialize sum to zero for numeric types
		
		for( const T &element : vec ){
			sum += element;
			
		}
		
		return sum;
		
	}
	
	template < typename T , size_t N >
	void pushArrayToVector( const T *arr , std::vector<T> &vec ) {
		
		/// Receives an std::vector and an array and pushes all the array elements to the vector.
		/// @param[in] arr input array.
		/// @param[in,out] vec std::vector container.
		/// @note This is a template function so it works with std::vector containers of any primitive type (i.e., integers, doubles, floats, etc.).
		
		for( size_t i = 0; i < N; ++i ) {
			vec.push_back( arr[i] );
		}
		
	}
	
	template <typename T>
	std::vector< T > getSubVectorFromVector( const std::vector< T > &vec , size_t start , size_t end ) {
		
		/// Receives an std::vector and returns a new std::vector containing all the elements of the original vector whose indexes are equal to or greater than "start" and smaller than "end".
		/// @param[in] vec input std::vector.
		/// @param[in] start starting index (i.e., include all elements whose indexes are equal to or greater than start).
		/// @param[in] end ending index (i.e., include elements whose indexes are smaller than end).
		/// @note This is a template function so it works with std::vector containers of any primitive type (i.e., integers, doubles, floats, etc.).
		
		if( start > vec.size() || end > vec.size() || start > end ) {
			allAbortWithMessage( MPI_COMM_WORLD , "Invalid range used in getSubVectorFromVector in utilities.h." );
			return { };
		}

		std::vector< T > sub_vector( vec.begin( ) + start , vec.begin( ) + end );
		return sub_vector;
	}
	
	//Elements in unordered set
	template < typename T , typename Hash >
	bool elementIsInUnorderedSet( const std::unordered_set< T , Hash > &set , const T &element ){
		
		/// Searches if an element is included in an unordered set.
		/// @param[in] set unordered set to be searched.
		/// @param[in] element element to be searched in the unordered set.
		/// @return true if the element is in the unordered set or false otherwise.
		/// @note This is a template function that allows searches on  unordered sets of any data type. Custom hash functions can also be used.
		
		return set.count( element ) > 0;
		
	}
	
	//Elements in unordered map
	template< typename KeyType , typename ValueType , typename Hash >
		bool mappingExists( const std::unordered_map< KeyType , ValueType , Hash > &map, const KeyType &key ) {
			
		/// Searches if an element is included in an unordered map.
		/// @param[in] map unordered_map to be searched.
		/// @param[in] key element to be searched in the unordered map.
		/// @return true if the element is in the unordered map or false otherwise.
		/// @note This is a template function that allows searches on  unordered maps of any data type. Custom hash functions can also be used.
		
		return map.count( key ) > 0;
	}
	
	
	//Elements in vector.
	template < typename T >
	bool elementIsInVector( const std::vector<T> &vec , const T &element) {
		
		/// Searches if an element is included in an std::vector.
		/// @param[in] vec vector to be searched.
		/// @param[in] element element to be searched in the vector.
		/// @return true if the element is in the vector or false otherwise.
		/// @note This is a template function that allows searches on vector of any primitive data types.
		
		return std::find( vec.begin( ), vec.end( ), element ) != vec.end();
		
	}
	
	
	//String management
	const std::string getConcatenatedStringFromStringsVector( const std::vector< std::string > &strings , size_t start , size_t end );
	const std::string getConcatenatedStringWithSpacesFromStringsVector( const std::vector< std::string > &strings , size_t start , size_t end );
	const int getStringPosInStringVec( const std::string &string2check , const std::vector< std::string > &strings );
	
	bool stringIsNumber( const std::string &string );
	bool stringIsIntNumber( const std::string &string );
	const bool stringIsBool( const std::string &string );
	const unsigned long int string2UnsignedLongInt( std::string &string );
	const int string2Int( std::string &string );
	const double string2Double( std::string &string );
	const bool string2Bool( std::string &string );
	const int boolString2Int( std::string &string );
	
	struct PairHash{
		
		/// @class PAPRECA::PairHash
		/// @brief Utility hash function (struct) for custom typedefs.
		/// @note This struct is necessary because an std pair is not a primitive datatype. This means that we have to provide a custom hash function for the unordered containers (both map and set). See: https://stackoverflow.com/questions/15160889/how-can-i-make-an-unordered-set-of-pairs-of-integers-in-c. Could also use std functions as a solution.
	
		inline std::size_t operator()(const std::pair<int,int> & v) const {
		return v.first*31+v.second;
		}
		
	};
	

	//Note 1: It is faster to search for a pair in a set and then map it (if it exists) than searching directly on a map. This is, of course, assuming that no RAM bottleneck exists (of course, storing both a map and an unordered set occupies more memory).
	//Note 2:  The unordered map gives you a pointer to a (parent) PredefinedReaction object. Use it as it is for bond breaking events, however, cast it to a (child) PredefinedBondForm object if necessary.
	//typedefs for identifying predefined events easily.
	
	typedef std::pair< int , int > INT_PAIR;
	typedef std::unordered_set< INT_PAIR , PairHash > PAIR_SET;
	
	typedef std::unordered_map< INT_PAIR , double , PairHash > INTPAIR2DOUBLE_MAP;
	typedef std::unordered_map< INT_PAIR , int , PairHash > INTPAIR2INT_MAP;
	
	typedef std::unordered_map< int , int > INT2INT_MAP;
	typedef std::unordered_map< int , INT2INT_MAP > INT2INTSMAP_MAP;
	typedef std::unordered_set< int > INT_SET;
	
	typedef std::unordered_set< LAMMPS_NS::tagint > TAGINT_SET;
	typedef std::unordered_map< LAMMPS_NS::tagint , int > TAGINT2INT_MAP;
	typedef std::vector< LAMMPS_NS::tagint > TAGINT_VEC;
	
	//Typedef for modern c++ containers of 3xdouble arrays
	typedef std::array< double , 3 > ARRAY3D;
	
	//Typedef for proc event selection
	typedef std::pair< double , int > DOUBLE2INTPAIR;
	typedef std::vector< DOUBLE2INTPAIR > DOUBLE2INTPAIR_VEC;
	
}//End of namespace PAPRECA







#endif
