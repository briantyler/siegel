/*
 * ### autorights BEGIN ###
 * siegel / siegelcl is an application that computes Siegel sets and cusp
 * lists for the action of SU(n,1;O) on HCn. Its development was supported
 * by EPSRC and University College London.
 * 
 * Copyright (C) 2009  Brian Tyler
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ### autorights END ###
 */

/**
 *  @brief    An inline header file for the initialization functions.
 *  @note     Include this file to set up a standard testing scheme.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-08-28
 */


#ifndef _SG_INITIALIZATION_H
#define _SG_INITIALIZATION_H 1

// Global includes
#include <cstddef>
#include <iostream>
#include <exception>
#include <sstream>
#include <set>

// Local includes
#include "engine/siegel_engine.hpp"
#include "geometry/algorithms/make_siegel.hpp"
#include "utility/pause.hpp"


namespace sg
{

/**
 *  @brief  Print help information to the screen.
 */
void help_information( void ) {
  using std::cout;
  using std::endl;
  
  cout << "Usage: siegel [options]\n"
          "Options:\n"
          "  -d=<arg>\tThe complex hyperbolic dimension of the space.\n"
          "\t\tDefault is 2.\n"
          "  -g=<arg>\tGenerator of the number field. Default is -1.\n"
          "  -b=<arg>\tBalance: determines how many cusps to generate\n"
          "\t\tbefore starting the search. Must be in (0, 1.0]\n"
          "\t\tDefault is 0.9.\n"
          "  -m=<arg>\tMatch: determines how accurate the result should be.\n"
          "\t\tMust be in (0, 1.0]. Default is 0.9.\n"
          "  -r=<arg>\tResolution: the number of subdivisions to divide\n"
          "\t\tthe search into (see --get-resolution)\n"
          "\t\tDefault is 1.\n"
          "  -f=<arg>\tFirst subdivision to search (first index is 1).\n"
          "\t\tDefault is 1.\n"
          "  -c=<arg>\tMaximum number of subdivisions to search.\n"
          "\t\tDefault is the maximum possible.\n"
          "  -s=<arg>\tSample: The number of sample points used in\n"
          "\t\tin estimating the minimum and maximum heights.\n"
          "\t\tDefault is 100,000.\n"
          "  -o=<arg>\tName of the output file. Default is /dev/null\n"
          "  --sieve=<arg>\tDepth of the cusp sieve. Default is m.\n"
	  "\t\t\tOptions are n=1, t=2, s=4, m=8, l=12\n"
          "  --get-resolution\tDisplays the computed resolution when used in\n"
          "\t\t\tconjunction with the -d, -g and -r parameters.\n"

#ifdef __GMP_PLUSPLUS__
          "  --use-gmp\t\tUse the GMP big number library (very slow).\n"
#endif

          "  --help\t\tDisplays this information.\n";
  cout << endl;
}


/**
 *  @struct bad_sg_input
 *  @brief  A lazy exception object
 */
struct bad_sg_input: public std::exception
{
  std::string parameters;
  
  virtual const char* what() const throw()
  {
    return parameters.c_str();
  }
  
  virtual ~bad_sg_input() throw() { }
};


/**
 *  @brief  Determine if the input generator is valid.
 *  @param  generator The generator of the number field.
 *  @note   Throws if \c generator is not a Heegner number.
 */
template <class _Integer>
    void validate_generator( const _Integer __generator )
{
  typedef _Integer integer_type;
  
  std::set<integer_type> generators;
  generators.insert( -1 );
  generators.insert( -2 );
  generators.insert( -3 );
  generators.insert( -7 );
  generators.insert( -11 );
  generators.insert( -19 );
  generators.insert( -43 );
  generators.insert( -67 );
  generators.insert( -163 );
  
  if( generators.count( __generator ) == 0 ) {
    std::stringstream ss;
    ss << "Bad generator input, generator given: " << __generator
       << ". Generator must be a Heegner number (-1, -2, -3, -7, -11, -19,"
          " -43, -67, -163).";
    
    bad_sg_input ex;
    ex.parameters = ss.str();
    throw ex;
  }
}


/**
 *  @brief  Determine if the balance is valid.
 *  @param  balance The balance variable.
 *  @note   Throws if \c balance is not in the required bounds.
 */
template <class _Float>
    void validate_balance( const _Float __balance )
{
  typedef _Float float_type;
  
  static const float_type minimum = 1e4 * utility::precision::zero();
  
  if( __balance < minimum || __balance > 1.0 ) {
    std::stringstream ss;
    ss << "Bad balance input, balance given: " << __balance
       << ". Balance must be in the range [" << minimum << ", 1.0].";
    
    bad_sg_input ex;
    ex.parameters = ss.str();
    throw ex;
  }
}


/**
 *  @brief  Determine if the match is valid.
 *  @param  match The match variable.
 *  @note   Throws if \c match is not in the required bounds.
 */
template <class _Float>
    void validate_match( const _Float __match )
{
  typedef _Float float_type;
  
  static const float_type minimum = 1e4 * utility::precision::zero();
  
  if( __match < minimum || __match > 1.0 ) {
    std::stringstream ss;
    ss << "Bad match input, match given: " << __match
       << ". Match must be in the range [" << minimum << ", 1.0].";
    
    bad_sg_input ex;
    ex.parameters = ss.str();
    throw ex;
  }
}

/**
 *  @brief  Determine if the balance is valid.
 *  @param  resolution The computed resolution.
 *  @param  first The first subspace to search.
 *  @note   Throws if \c first >= \c resolution
 */
template <class _Integer>
    void validate_first( const _Integer __resolution, const _Integer __first )
{
  typedef _Integer integer_type;
  
  if( __first >= __resolution ) {
    std::stringstream ss;
    ss << "Bad first input, first given: " << __first
       << ". First must be in the range [0, " << ( __resolution - 1 ) << "].";
    
    bad_sg_input ex;
    ex.parameters = ss.str();
    throw ex;
  }
}


/**
 *  @brief  A helper function to find the computed resolution for an input
 *          resolution.
 *  @param  N The dimension of the space.
 *  @param  Integer The integer type.
 *  
 *  @param  generator The generator of the number field.
 *  @param  resolution The minimum number of subspaces to subdivide the main
 *          Siegel container into.
 *  @return Nothing. Outputs the computed resolution to the screen.
 *  @note   This is a helperfunctin used to set up an iteeffectn scheme.
 */
template <std::size_t N, class _Integer>
    void get_resolution
    ( const _Integer __generator, const std::size_t __resolution )
{
  using std::cout;
  using std::endl;
  using std::size_t;
  using structure::geometric::hyperbolic::heisenberg_slice;
  using structure::numerical::iq_field;
  
  typedef double float_type;
  typedef _Integer integer_type;
  
  typedef heisenberg_slice<N,float_type,std::size_t> slice_type;
  typedef typename slice_type::space_type space_type;
  typedef iq_field<double,integer_type> field_type;
  
  // Set the field generator
  field_type& field = field_type::instance();
  field.initialize( __generator );
  
  slice_type slice;
  geometry::algorithms::make_siegel<field_type,space_type>( slice.space() );
  slice.resolution() = __resolution;
  slice.initialize();
  
  cout << "Total Space:\n" << slice.space() << endl;
  cout << "Resolution: " << slice.resolution() << endl;
}


/**
 *  @brief  A helper function for creating Siegel Sets
 *  @param  Depth The depth of the sieve to use.
 *  @param  N The dimension of the space.
 *  @param  Float The floating point type.
 *  @param  Integer The integer type.
 *  
 *  @param  generator The generator of the number field.
 *  @param  balance The \c balance determines the effect of the maximum and
 *          minimum heights during initialization. It should be in the set
 *          \f$ \[1,0) \f$.
 *  @param  match The \c balance determines the effect of the maximum and
 *          minimum heights during computation. It should be in the set
 *          \f$ \[1,0) \f$.
 *  @param  resolution Determines the number of subspaces to subdivide the
 *          main Siegel container into. Note that since a heisenberg slice is
 *          used, the only guaranteed is that the final resolution will be at
 *          least as big as the input value (probably larger). Defaults to
 *          \c 1, in this case no slicing is done and the entire Siegel
 *          container is searched.
 *  @param  first The subspace to begin the search in, must be less than
 *          the computed resolution (but slightly confusingly not necessarily
 *          less than the input resolution) this allows a large search scheme
 *          to be run on multiple computers simultaneously.
 *  @param  count The number of iteeffectns to perform. If 0 then the maximum
 *          number of iteeffectns will be performed, otherwise min( first + count,
 *          [computed resolution] ) will be performed.
 *  @param  sample Determines how many random points to sample when
 *          estimating a minimum and maximum height, clearly the larger
 *          sample is, the more accurate the result, but the slower the
 *          algorithm during the sampling phase. Something like 10^6 gives
 *          a good balance between accuracy and speed.
 */
template <std::size_t _Depth, std::size_t N, class _Float, class _Integer>
    void compute_siegel
        ( std::ostream& __output, const _Integer __generator,
          const _Float __balance, const _Float __match,
          const std::size_t __resolution = 1, const std::size_t __first = 0,
          const std::size_t __count = 0, const std::size_t __sample = 100000ul
        )
{
  using std::cout;
  using std::endl;
  using std::size_t;
  using engine::siegel_engine;
  using structure::geometric::hyperbolic::heisenberg_slice;
  
  typedef _Float float_type;
  typedef _Integer integer_type;
  
  typedef siegel_engine<_Depth,N,float_type,integer_type,0> engine_type;
  typedef heisenberg_slice<N,float_type,std::size_t> slice_type;
  typedef typename engine_type::space_type space_type;
  typedef typename engine_type::field_type field_type;
  
  // Validate inputs
  validate_generator( __generator );
  validate_balance( __balance );
  validate_match( __match );
  
  // Set the field generator
  field_type& field = field_type::instance();
  field.initialize( __generator );
  
  // Initialize the slice in case the resolution is greater than 1.
  slice_type slice;
  if( __resolution > 1 ) {
    // Set the outer search space and resolution.
    slice.resolution() = __resolution;
    geometry::algorithms::make_siegel<field_type,space_type>( slice.space() );
    slice.initialize();
    
    // Ensure that the first subspace is in the valid range.
    validate_first( slice.resolution(), __first );
  }
  
  engine_type engine;
  engine.sample() = __sample;
  engine.balance() = __balance;
  engine.match() = __match;
  
  __output << "Siegel Set Generator v1.0 (written by Brian Tyler)" << endl;
  __output << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" << endl;
  
  // Output the setup information
  cout << "Dimension: " << engine_type::dimension << endl;
  cout << "Depth: " << engine_type::depth << endl;
  cout << "Generator: " << field.generator() << endl;
  cout << "Balance: " << engine.balance() << endl;
  cout << "Match: " << engine.match() << endl;
  cout << "Sample: " << engine.sample() << endl;
  
  // Log the setup information
  __output << "Dimension: " << engine_type::dimension << endl;
  __output << "Depth: " << engine_type::depth << endl;
  __output << "Generator: " << field.generator() << endl;
  __output << "Balance: " << engine.balance() << endl;
  __output << "Match: " << engine.match() << endl;
  __output << "Sample: " << engine.sample() << endl;
  
  // Begin timing the operation.
  boost::progress_timer t;
  
  // Set the output stream.
  engine.set_output( __output );
  
  
  if( __resolution <= 1 ) {
    // Set the search space
    geometry::algorithms::make_siegel<field_type,space_type>( engine.space() );
    cout << "Space:\n" << engine.space() << endl;
    __output << "Space:\n" << engine.space() << "\n" << endl;
    
    // Wait a few seconds to give the user a chance to read the details
    //utility::pause()( 5 );
    cout << endl;
    
    // Initialize the engine
    engine.initialize();
    
    // And begin the search.
    engine();
  }
  else {
    // Output the search space and the number of subdivisions
    cout << "Space:\n" << slice.space() << endl;
    cout << "Resolution: " << slice.resolution() << endl;
    
    // Log the search space and the number of subdivisions
    __output << "Space:\n" << slice.space() << endl;
    __output << "Resolution: " << slice.resolution() << endl;
    
    size_t last = (   __count == 0
                    ? slice.resolution()
                    : std::min( slice.resolution(), __first + __count )
                  );
    
    // Output the indices of the first and last subspaces.
    cout << "First: " << __first << endl;
    cout << "Last: " << last << endl;
    
    // Log the indices of the first and last subspaces.
    __output << "First: " << __first << endl;
    __output << "Last: " << last << endl;
    
    // Perform the search over all subspaces between First and Last.
    float_type minHeight = 2.0;
    for( size_t i = __first; i != last; ++i ) {
      // Set the engine space.
      engine.space() = slice.subspace_at(i);
      cout << "Subspace (" << (i+1) << " of " << slice.resolution() << ")\n"
           << engine.space() << "\n" << endl;
      
      __output << "Subspace (" << (i+1) << " of " << slice.resolution() << ")\n"
               << engine.space() << "\n" << endl;
      
      cout << endl;
      
      // Initialize the engine
      engine.initialize();
      
      // And compute the minimum for this loop
      float_type tmp = engine();
      
      // If this is a new minimum then update.
      if( tmp < minHeight ) { minHeight = tmp; }
      cout << endl;
      __output << endl;
    }
    
    // Output the minimum height found.
    cout << "\nMinimum Height: " << minHeight << endl;
    __output << "\nMinimum Height: " << minHeight << endl;
  }
}

}// sg
#endif
