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
 *  @file     cusp_candidates_engine.hpp
 *  @brief    An inline header file for the \c cusp_candidates_engine class.
 *  @note     An engine which constructs a complete list of candidate
 *            cusps which are effective on a given region.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2009-04-30
 */


#ifndef _SG_CUSP_CANDIDATES_ENGINE_H
#define _SG_CUSP_CANDIDATES_ENGINE_H 1

// Global includes
#include <iostream>
#include <sstream>
#include <vector>
#include <map>

#include <boost/timer.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "engine/cusp_engine.hpp"
#include "engine/cusp_validator.hpp"


namespace sg
{
namespace engine
{
/**
 *  @class    cusp_candidates_engine
 *  @note     An engine which constructs a complete list of candidate cusps
 *            or a given space and group.
 *  @param    N The dimension of the space
 *  @param    Float A floating point type; defaults to \c double
 *  @param    Integer An integral type; defaults to \c long
 *  @param    Id The identifier of the field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2009-04-30
 */
template < std::size_t N, class _Float = double, class _Integer = long,
           std::size_t _Id = 0 >
         class cusp_candidates_engine
{
  public:
    //! The object type
    typedef cusp_candidates_engine<N,_Float,_Integer,_Id> self_type;
    //! The float type
    typedef _Float float_type;
    //! The integral type
    typedef _Integer integer_type;
    //! The type of the cusp engine
    typedef cusp_engine<N,float_type,integer_type,_Id> cusp_engine_type;
    //! The cusp selector engine type
    typedef cusp_validator<N,float_type,integer_type,_Id> validator_type;
    //! The cusp type
    typedef typename cusp_engine_type::cusp_type cusp_type;
    //! The space type
    typedef typename cusp_engine_type::space_type space_type;
    //! The point type
    typedef typename space_type::point_type point_type;
    //! The field type
    typedef typename cusp_engine_type::field_type field_type;
    //! A vector to hold all of the generated cusps.
    typedef std::vector<cusp_type> cusp_vector_type;
    
    //! Template constants
    enum constants {
      dimension_size = N
    };
    
    
  private:
    //! Floating point param type
    typedef typename boost::call_traits<float_type>::param_type flt_param_type;
    //! Integer param type
    typedef typename boost::call_traits<integer_type>::param_type int_param_type;
    
    //! A map to hold the number of cusps of a given dilation factor.
    typedef std::map<integer_type,unsigned long> count_map_type;
    //! The map value type
    typedef typename count_map_type::value_type count_type;
    
    
  public:
    /**
     *  @brief  Default constructor
     */
    cusp_candidates_engine( )
  : space_(), engine_(), validator_()
    { engine_.bind_space(space_); }
    
    
    /**
     *  @brief  Get a constant reference to the underlying space.
     *  @return A constant reference to the underlying space.
     */
    const space_type& space( void ) const { return space_; }
    
    /**
     *  @brief  Get a reference to the underlying space.
     *  @return A reference to the underlying space.
     */
    space_type& space( void ) { return space_; }
    
    
    /**
     *  @brief  Set the imaginary quadratic field generator.
     */
    void set_field_generator( int_param_type __generator ) {
      field_type::instance().initialize( __generator );
    }
    
    
    /**
     *  @brief  Set the output stream.
     */
    void set_output( std::ostream& __output ) { output_ = &__output; }
    
    
    /**
     *  @brief  Get the vector of cusp candiates.
     *  @note   Only constant access is supplied because you shouldn't change
     *          the cusps outside the class, you should only observe them.
     */
    const cusp_vector_type& cusps() const { return cusps_; }
    
    
    /**
     *  @brief  Get the total number of cusps
     */
    unsigned long total_cusps() const {
      unsigned long count = 0;
      for( typename count_map_type::const_iterator it = counts_.begin();
           it != counts_.end(); ++it )
      { count += it->second; }
      
      return count;
    }
    
    
    /**
     *  @brief  Get the total number of cusps
     */
    unsigned long total_cusps_of_dilation( integer_type __dilation ) const {
      typename count_map_type::const_iterator it = counts_.find( __dilation );
      
      return it != counts_.end() ? it->second : 0;
    }
    
    
    /**
     *  @brief  Start the engine to find a candidate cusp list.
     */
    void operator()( bool __store = true ) {
      typedef utility::math::is_zero<float_type> _isz;
      store_ = __store;
      
      sys_initialize();
      
      if( _isz()( space_.height_ref().lower() ) ) {
        sys_output_message( "Aborted: Height bound must be set prior to use" );
        return;
      }
      
      sys_output_message( "Cusp Generation" );
      sys_output_message( "^^^^^^^^^^^^^^^" );
      
      sys_loop();
      sys_output_message( "Done!\n" );
      sys_clear_ss();
      ss_ << "Total number of cusp candidates = " << total_cusps() << "\n";
      sys_output_message( ss_.str() );
    }
    
    
    void tex() const {
      *output_ << "\n\n";
      *output_ << "LaTeX Output\n";
      *output_ << "^^^^^^^^^^^^\n\n";
      // Title
      *output_ << "\\noindent{\\raggedright Space = $"
              << space_.tex_complex_hyperbolic() 
              << "$, Group = $" << engine_.tex_group() << "$, "
              << "Total number of cusps $=" << total_cusps() << "$\n";
      
      // Output Cusps
      integer_type currentDilationFactor = 0;
      typedef typename cusp_vector_type::const_iterator _citer;
      for( _citer it = cusps_.begin(); it != cusps_.end(); ++it ) {
        if( it->dilation() != currentDilationFactor ) {
          *output_ << "\\ \\\\$\\Delta = " << it->dilation()
                   << "$, Number of cusps candidates$=" 
                   << total_cusps_of_dilation( it->dilation() )
                   << "$ :\\\\\n";
          
          currentDilationFactor = it->dilation();
        }
        *output_ << "$" << it->tex() << "$\\ \\ \n";
      }
      *output_ << "}\n";
      *output_ << std::flush;
    }
    
    
  private:
    /**
     *  @brief  Initialize the engine ready for use
     */
    void sys_initialize( void ) {
      using std::cout;
      using std::endl;
      
      sys_output_message( "Initialization" );
      sys_output_message( "^^^^^^^^^^^^^^" );
      
      // Clear all cusp data
      cusps_.clear();
      
      // Initialize the space in case it has not already been done.
      space_.initialize();
      sys_output_message( "Space initialized!" );
      
      // Initialize the engine and move to the first valid cusp.
      
      engine_.initialize(1);
      engine_();
      sys_output_message( "Engine initialized!" );
      
      // Prepare the cusp validator
      validator_.bind_space( space_ );
      validator_.bind_cusp( engine_.cusp() );
      sys_output_message( "Validator initialized!\n" );
       
      sys_short_precision();
      sys_clear_ss();
      ss_  << "Initialization complete.\n"
           << "--\n\n";
      sys_output_message( ss_.str() );
      sys_long_precision();
    }
    
    
    /**
     *  @brief  Begin the search loop.
     */
    void sys_loop( void ) {
      boost::timer timer;
      
      // The engine knows when it can't generate more cusp candidates for the
      // starting space, so it is left to the engine to decide when to stop
      // processing.
      integer_type currentDilationFactor = 0;
      unsigned long count = 0;
      while( !engine_.finished() ) {
        // If the dilation factor changes notify the user
        if( engine_.cusp().dilation() != currentDilationFactor ) {
          sys_clear_ss();
          ss_  << "Number of cusps with dilation factor " << currentDilationFactor
               << " = " << count << "\n"
               << "Time Elapsed: " << timer.elapsed() << "s\n"
               << "Dilation factor increased to: " << engine_.cusp().dilation();
          sys_output_message( ss_.str() );
          counts_.insert( count_type( currentDilationFactor, count) );
          currentDilationFactor = engine_.cusp().dilation();
          count = 0;
        }
        
        if( validator_() ) {
          ++count;
          if( store_ ) {
            cusps_.push_back( engine_.cusp() );
          }
        }
        
        // Build the next cusp
        engine_();
      }
      sys_clear_ss();
      ss_  << "Number of cusps with dilation factor " << currentDilationFactor
           << " = " << count << "\n"
           << "Total Time Elapsed: " << timer.elapsed() << "s\n";
      sys_output_message( ss_.str() );
      counts_.insert( count_type( currentDilationFactor, count) );
    }
    
    
    void sys_output_message( std::string __message ) {
      using std::cout;
      using std::endl;
      
      cout << __message << endl;
      *output_ << __message << endl;
    }
    
    void sys_clear_ss( void ) { ss_.str( std::string() ); }
    
    void sys_short_precision( void ) {
      std::cout.precision( 2 );
      ss_.precision( 2 );
      output_->precision( 2 );
    }
    
    void sys_long_precision( void ) {
      std::cout.precision( utility::precision::stream() );
      ss_.precision( utility::precision::stream() );
      output_->precision( utility::precision::stream() );
    }
    
    
    //! The space that is to become the Siegel set
    space_type space_;
    //! The engine which builds the cusps
    cusp_engine_type engine_;
    //! The cusp validator which validates all candiate cusps
    validator_type validator_;
    
    //! List of candidate cusps
    cusp_vector_type cusps_;
    //! The total number of cusps
    count_map_type counts_;
    
    //! The output stream to output to
    std::ostream* output_;
    //! A stringstream to write things to
    std::stringstream ss_;
    //! Should the engine store the cusps: true, or count them:false?
    bool store_;
};

}// engine
}// sg
#endif
