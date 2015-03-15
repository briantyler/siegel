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
 *  @file     truncate.hpp
 *  @brief    This is a header implementation file for truncate
 *  @note     Include this file to truncate floating point numbers to integers.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-14
 */


#ifndef _SG_TRUNCATE_H
#define _SG_TRUNCATE_H 1

// Global includes
#include <cmath>

#include <boost/call_traits.hpp>

// Local includes
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
   *  @struct  truncate
   *  @brief   Truncate a floating point number and convert to an integer type.
   *  @author  Brian Tyler
   *  @version 1.0
   *  @date    2008-03-14
   *  @param   To The integer type to convert to.
   *  @param   From The floating point type to convert from.
   *  @note    Specialised for native and gmp types.
 */
template <class _To, class _From> struct truncate
  : public std::unary_function
      <typename boost::call_traits<_From>::param_type,_To>
{
  //! Object type
  typedef truncate<_To,_From> self_type;
  //! Type to convert to
  typedef _To to_type;
  //!Type to convert from
  typedef _From from_type;
  
  /**
   *  @brief  Truncate a float and convert to an integral type.
   *  @param  value The floating point to truncate.
   *  @return The truncated value cast to an integral type.
   */
  _To operator( ) ( typename self_type::argument_type __value ) const {
    _To truncated;
    functors::stream_cast<_From,_To>( )( __value, truncated );
    return truncated;
  }
};


template <> struct truncate<double,double>
  : public std::unary_function <double,double>
{
  //! Object type
  typedef truncate<double,double> self_type;
  //! Type to convert to
  typedef double to_type;
  //!Type to convert from
  typedef double from_type;
  
  double operator( ) ( double __value ) const {
    double truncated;
    std::modf( __value, &truncated );
    return truncated;
  }
};


template <> struct truncate<long double,long double>
  : public std::unary_function <long double,long double>
{
  //! Object type
  typedef truncate<long double,long double> self_type;
  //! Type to convert to
  typedef long double to_type;
  //!Type to convert from
  typedef long double from_type;
  
  long double operator( ) ( long double __value ) const {
    long double truncated;
    std::modf( __value, &truncated );
    return truncated;
  }
};


#ifdef __GMP_PLUSPLUS__
template <> struct truncate<mpf_class,mpf_class>
  : public std::unary_function <const mpf_class&,mpf_class>
{
  //! Object type
  typedef truncate<mpf_class,mpf_class> self_type;
  //! Type to convert to
  typedef mpf_class to_type;
  //!Type to convert from
  typedef mpf_class from_type;
  
  mpf_class operator( ) ( const mpf_class& __value ) const {
    return mpf_class( ::trunc( __value ) );
  }
};


template <> struct truncate<mpq_class,mpq_class>
  : public std::unary_function <const mpq_class&,mpq_class>
{
  //! Object type
  typedef truncate<mpq_class,mpq_class> self_type;
  //! Type to convert to
  typedef mpq_class to_type;
  //!Type to convert from
  typedef mpq_class from_type;
  
  mpq_class operator( ) ( const mpq_class& __value ) const {
    return mpq_class( ::trunc( mpf_class(__value) ) );
  }
};
#endif

}// math
}// utility
}// sg
#endif
