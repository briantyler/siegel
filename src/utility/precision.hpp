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
 *  @file     precision.hpp
 *  @brief    This is a header implementation file for precision.
 *  @note     Include this file to set the precision of streams and other
 *            objects.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-05
 */


#ifndef _SG_PRECISION_H
#define _SG_PRECISION_H 1

// Global includes
#include <sstream>


namespace sg
{
namespace utility
{
/**
 *  @struct   precision
 *  @brief    A singleton struct containing information about precision
 *  @note     Contains the precision rules for all type conversion.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-05
 *  @note     Probably want to set the stream precision higher than the
 *            zero precision, since that should keep things fairly accurate
 *            on a restart.
 */
struct precision {
  //! Type of the stream precision variable
  typedef std::streamsize stream_prec_type;
  //! Type of the floating point precision variable
  typedef double float_prec_type;
  
  /**
   *  @brief  Get the stream precision
   */
  static std::streamsize stream( void ){ return instance().stream_; }
  
  /**
   *  @brief  Set the stream precision
   */
  static void set_stream( stream_prec_type __stream ) {
    instance().stream_ = __stream;
  }
  
  
  /**
   *  @brief  Get the computational zero
   */
  static double zero( void ){ return instance().zero_; }
  
  /**
   *  @brief  Set the computational zero
   */
  static void set_zero( float_prec_type __zero ) {
    instance().zero_ = __zero;
  }
   
  private:
    /**
     *  @brief  Singleton constructor
     */
    static precision& instance( void ){ static precision p; return p; }
    
    /* private construct, copy, assign and delete ensure that the object
     * is a singleton.
     */
    precision( ) : stream_(12), zero_(1e-10) { }
    precision( const precision& that )
    : stream_( that.stream_ ), zero_( that.zero_ ) { }
    precision& operator=( const precision& that ) {
      stream_ = that.stream_; zero_ = that.zero_; return *this;
    }
    ~precision( ) { }
    
    //! Precision of stream operators
    stream_prec_type stream_;
    //! The value which is taken to be the computational zero
    float_prec_type zero_;
};

}// utility
}// sg
#endif
